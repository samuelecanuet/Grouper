#include "Merger.hh"
#include <ctime>

int main(int argc, char *argv[])
{

    if (argc == 2)
    {
        if (argv[1] == "Al")
        {
            Catcher = "Al";
        }
        else
        {
            Catcher = "Mylar";
        }
    }
    else
    {
        Catcher = "Al";
    }

    clock_t start = clock(), Current;

    InitFiles();
    string dirNameMerged = "./Merged/";
    string dirNameCleaned = "./Cleaned/";
    string dirNameMatched = "./Matched/";

    Merged_File = new TFile((dirNameMerged + "runs_32Ar_merged.root").c_str(), "RECREATE");
    Tree_Merged = new TTree("Tree", "Tree");
    Tree_Merged->Branch("Tree_Silicon", &Tree_Merged_Silicon);
    Tree_Merged->Branch("Tree_SiPM", &Tree_Merged_SiPM);

    InitDetectors("../Analysis/ReadFaster/Data/detectors.dat");
    InitPeakWindow(); /////// find window range for peaks in config file
    InitSelection();

    for (auto &run : runs)
    {
        Matched_File = new TFile((dirNameMatched + "run_0" + to_string(run) + "_32Ar_matched.root").c_str(), "READ");
        SetTime(Matched_File);
    }

    ///////////////////// CREATE TIME RUNs FILES & BUILD PROTON HISTOGRAM FOR CALIBRATION ///////////////
    InitHistograms();
    InitCalib();       /////// find detector calib from simulation
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////SILICON CALIBRATION////////
    for (auto &run : runs)
    {
        cout << "<SAM> Silicon calibration on run : " << run << " | " << endl;
        Matched_File = new TFile((dirNameMatched + "run_0" + to_string(run) + "_32Ar_matched.root").c_str(), "READ");
        Tree_Read = (TTree *)Matched_File->Get("Tree");
        TTreeReader *Reader = new TTreeReader(Tree_Read);
        TTreeReaderArray<Signal> Tree_Silicon(*Reader, "Tree_Silicon");

        pair<string, string> string_time = GetTime(Matched_File);
        start_time = Convert_DatetoTime(string_time.first, file_string_Time[0].first);
        stop_time = Convert_DatetoTime(string_time.second, file_string_Time[0].first);

  
        InitMatching();

        bool CurrentDetectorSelection[SIGNAL_MAX];
        copy(begin(RunDetectorSelection[run]), end(RunDetectorSelection[run]), begin(CurrentDetectorSelection));

        int TotalEntries = Tree_Read->GetEntries();
        while (Reader->Next())
        {
            
            ULong64_t cEntry = Reader->GetCurrentEntry();
            if (cEntry % 100 == 0 && cEntry > 0)
            {
                ProgressBar(cEntry, TotalEntries, start, Current);
            }

            if (CurrentDetectorSelection[Tree_Silicon[0].Label])
            {
                
                HSilicon_Channel_Unmatched[Tree_Silicon[0].Label]->Fill(1 / detectorMatching[Tree_Silicon[0].Label] * Tree_Silicon[0].Channel);
                HSilicon_Channel_Matched[Tree_Silicon[0].Label]->Fill(Tree_Silicon[0].Channel);
                HSilicon_Channel_Matched_CURRENTRUN[Tree_Silicon[0].Label]->Fill(Tree_Silicon[0].Channel);

                HSiliconRun_Channel_Unmatched[Tree_Silicon[0].Label]->Fill(start_time + (double)Tree_Silicon[0].Time * 1e-9 / 3600, 1 / detectorMatching[Tree_Silicon[0].Label] * Tree_Silicon[0].Channel);
                HSiliconRun_Channel_Matched[Tree_Silicon[0].Label]->Fill(start_time + (double)Tree_Silicon[0].Time * 1e-9 / 3600, Tree_Silicon[0].Channel);
            }
        }
        Matched_File->Close();
        MakeSiliconCalibration(run);
        cout<<endl;
    }

    MakeSiliconCalibration(); /////// compute calibration for merged coeficient

    ////////SiPM CALIBRATION////////
    // SiPM_Runs_Tree = new TTree("SiPM_Runs_Tree", "SiPM_Runs_Tree");
    // TList *list = new TList;
    for (auto &run : runs)
    {
        run = 30;
        cout << "<SAM> SiPM calibration on run : " << run << " | ";
        Matched_File = new TFile((dirNameMatched + "run_0" + to_string(run) + "_32Ar_matched.root").c_str(), "READ");
        Tree_Read = (TTree *)Matched_File->Get("Tree");
        TTreeReader *Reader = new TTreeReader(Tree_Read);
        TTreeReaderArray<Signal> Tree_Silicon(*Reader, "Tree_Silicon");
        TTreeReaderArray<Signal> Tree_SiPM(*Reader, "Tree_SiPM");
        Ref_Hist = new TH1D("Ref_Hist", "Ref_Hist", eLowN, 0, eLowMax/1000);
        while (Reader->Next())
        {
            if (Is_F(1 / SiliconCalib[99][Tree_Silicon[0].Label].first * Tree_Silicon[0].Channel, Tree_Silicon[0].Label))
            {
                for (auto &SiPM : Tree_SiPM)
                {
                    Ref_Hist->Fill(SiPM.Channel/1000);
                }
            }
        }
        
        // if (run == runs[0])
        // {
        //     SiPM_Runs_Tree = Tree_Read->CloneTree(0);
        // }
        // else
        // {
        //     list->Add(Tree_Read);
        // }
        pair<string, string> string_time = GetTime(Matched_File);
        start_time = Convert_DatetoTime(string_time.first, file_string_Time[0].first);
        stop_time = Convert_DatetoTime(string_time.second, file_string_Time[0].first);

        bool CurrentDetectorSelection[SIGNAL_MAX];
        copy(begin(RunDetectorSelection[run]), end(RunDetectorSelection[run]), begin(CurrentDetectorSelection));

        MakeSiPMCalibration(run);
        Matched_File->Close();
        cout << endl;
    }

    // SiPM_Runs_Tree->Merge(list);
    // MakeSiPMCalibration(); /////// compute calibration for merged coeficient

    ////////MERGING////////

    for (auto &run : runs)
    {
        cout << "<SAM> Merging on run : " << run << endl;
        Matched_File = new TFile((dirNameMatched + "run_0" + to_string(run) + "_32Ar_matched.root").c_str(), "READ");
        Tree_Read = (TTree *)Matched_File->Get("Tree_LowHighMatched");
        TTreeReader *Reader = new TTreeReader(Tree_Read);
        TTreeReaderArray<Signal> Tree_Silicon(*Reader, "Tree_Silicon");
        TTreeReaderArray<Signal> Tree_SiPM(*Reader, "Tree_SiPM");

        pair<string, string> string_time = GetTime(Matched_File);
        start_time = Convert_DatetoTime(string_time.first, file_string_Time[0].first);
        stop_time = Convert_DatetoTime(string_time.second, file_string_Time[0].first);

        InitMatching();
        bool CurrentDetectorSelection[SIGNAL_MAX];
        copy(begin(RunDetectorSelection[run]), end(RunDetectorSelection[run]), begin(CurrentDetectorSelection));

        int TotalEntries = Tree_Read->GetEntries();
        while (Reader->Next())
        {
            // ULong64_t cEntry = Reader->GetCurrentEntry();
            // if (cEntry % 10000 == 0 && cEntry > 10000)
            // {
            //     ProgressBar(cEntry, TotalEntries, start, Current);
            // }

            if (CurrentDetectorSelection[Tree_Silicon[0].Label])
            {
                double Silicon_Energy = 1/SiliconCalib[99][Tree_Silicon[0].Label].first * Tree_Silicon[0].Channel;


                HSiliconRun[Tree_Silicon[0].Label]->Fill(start_time + (double)Tree_Silicon[0].Time * 1e-9 / 3600, Silicon_Energy);
                HSilicon[Tree_Silicon[0].Label]->Fill(Silicon_Energy);

                if (Tree_SiPM.GetSize() > 3)
                {
                    HSilicon_coinc[Tree_Silicon[0].Label]->Fill(Silicon_Energy);
                }
                else
                {
                    HSilicon_no_coinc[Tree_Silicon[0].Label]->Fill(Silicon_Energy);
                }

                if (Tree_SiPM.GetSize() > 0)
                {
                    double sum = 0;
                    for (int j = 0; j < Tree_SiPM.GetSize(); j++)
                    {
                        sum += Tree_SiPM[j].Channel;
                    }

                    for (int i = 1; i < Tree_SiPM.GetSize() + 1; i++)
                    {
                        HSiPM[i]->Fill(sum);
                    }

                    if (Is_F(Silicon_Energy, Tree_Silicon[0].Label))
                    {
                        for (int i = 1; i < Tree_SiPM.GetSize() + 1; i++)
                        {
                            HSiPM_F[i]->Fill(sum);
                        }
                    }

                    else if (Is_GT(Silicon_Energy, Tree_Silicon[0].Label))
                    {
                        for (int i = 1; i < Tree_SiPM.GetSize() + 1; i++)
                        {
                            HSiPM_GT[i]->Fill(sum);
                        }
                    }
                }
            }
        }
        Matched_File->Close();
    }

    WriteHistograms();
    Merged_File->Close();

    return 0;
}