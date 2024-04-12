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

    dark_file = new TFile("../../SiPM_model/test.root", "UPDATE");
    dark_file->cd();
    dark = (TF1 *)dark_file->Get("f1");
    
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
    InitCalib(); /////// find detector calib from simulation
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////SILICON CALIBRATION////////
    for (auto &run : runs)
    {
        string Prefix = Form("<SAM> Silicon calibration on run : %d | ", run);
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
            if ((cEntry % 100 == 0 && cEntry > 0) || (cEntry == TotalEntries-1))
            {
                ProgressBar(cEntry, TotalEntries, start, Current, Prefix);
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
    }

    MakeSiliconCalibration(); /////// compute calibration for merged coeficient
    InitCalibSiPM();          /////// find SiPM calib from simulation

    ////////SiPM CALIBRATION////////
    // SiPM_Runs_Tree = new TTree("SiPM_Runs_Tree", "SiPM_Runs_Tree");
    // TList *list = new TList;
    for (int i = 0; i < BETA_SIZE+1; i++)
    {
        for (int j = 0; j < BETA_SIZE+1; j++)
        {
            Ref_Hist_F[j][i] = new TH1D(("Fermi_M"+ to_string(j)+"_SiPM" + to_string(i)).c_str(), ("Fermi_M"+ to_string(j)+"_SiPM" + to_string(i)).c_str(), 800, 0, 8000);
        }
        Ref_Hist_Multi_F[i] = new TH1D(("Fermi_M" + to_string(i)).c_str(), ("Fermi_M" + to_string(i)).c_str(), 2000, 0, 8000);
        Ref_Hist_Surface_F[i] = new TH1D(("Fermi_Surface" + to_string(i)).c_str(), ("Fermi_Surface" + to_string(i)).c_str(), 2000, 0, 8000);
    }

    double Channel;
    int Label;
    for (int i = 0; i < BETA_SIZE+1; i++)
    {
        // Ref_Tree_F[i] = new TTree(("Fermi_Tree_SiPM" + to_string(i)).c_str(), ("Fermi_Tree_SiPM" + to_string(i)).c_str());
        // Ref_Tree_F[i]->Branch("Channel", &Channel, "Channel/D");
        for (int j = 0; j < BETA_SIZE+1; j++)
        {
            Ref_Tree_F[j][i] = new TTree(("Fermi_Tree_M"+ to_string(j)+"_SiPM" + to_string(i)).c_str(), ("Fermi_Tree_M"+ to_string(j)+"_SiPM" + to_string(i)).c_str());
            Ref_Tree_F[j][i]->Branch("Channel", &Channel, "Channel/D");
        }
        Ref_Tree_Multi_F[i] = new TTree(("Fermi_Tree_M" + to_string(i)).c_str(), ("Fermi_Tree_M" + to_string(i)).c_str());
        Ref_Tree_Multi_F[i]->Branch("Channel", &Channel, "Channel/D");
        Ref_Tree_Multi_F[i]->Branch("Label", &Label, "Label/I");

        Ref_Tree_Surface_F[i] = new TTree(("Fermi_Tree_Surface" + to_string(i)).c_str(), ("Fermi_Tree_Surface" + to_string(i)).c_str());
        Ref_Tree_Surface_F[i]->Branch("Channel", &Channel, "Channel/D");
        Ref_Tree_Surface_F[i]->Branch("Label", &Label, "Label/I");
    }

    for (auto &run : runs)
    {
        string Prefix = Form("<SAM> SiPM calibration on run : %d | ", run);
        Matched_File = new TFile((dirNameMatched + "run_0" + to_string(run) + "_32Ar_matched.root").c_str(), "READ");
        Tree_Read = (TTree *)Matched_File->Get("Tree");
        TTreeReader *Reader = new TTreeReader(Tree_Read);
        TTreeReaderArray<Signal> Tree_Silicon(*Reader, "Tree_Silicon");
        TTreeReaderArray<Signal> Tree_SiPM(*Reader, "Tree_SiPM");

        vector<pair<int, double>> vec;
        vector<int> result;
        int TotalEntries = Tree_Read->GetEntries();

        while (Reader->Next())
        {
            ULong64_t cEntry = Reader->GetCurrentEntry();
            if ((cEntry % 100 == 0 && cEntry > 0) || (cEntry == TotalEntries-1))
            {
                ProgressBar(cEntry, TotalEntries, start, Current, Prefix);
            }

            if (Is_F(1 / SiliconCalib[99][Tree_Silicon[0].Label].first * Tree_Silicon[0].Channel, Tree_Silicon[0].Label))
            {
                for (int i = 1; i <= Tree_SiPM.GetSize(); i++)
                {
                    for (auto &SiPM : Tree_SiPM)
                    {
                        Channel = SiPM.Channel / 1000;
                        Ref_Tree_F[i][GetDetectorChannel(SiPM.Label)]->Fill();
                    }
                }

                for (int i = 1; i <= Tree_SiPM.GetSize(); i++)
                {
                    for (auto &SiPM : Tree_SiPM)
                    {
                        Channel = SiPM.Channel / 1000;
                        Label = GetDetectorChannel(SiPM.Label);
                        Ref_Tree_Multi_F[i]->Fill();
                    }
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

        Matched_File->Close();
    }
    Merged_File->cd();
    MakeSiPMCalibration(0);
    MakeSiPM_SiPMPlots();
    MakeSiPM_MultiplicityPlots();

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
                double Silicon_Energy = 1 / SiliconCalib[99][Tree_Silicon[0].Label].first * Tree_Silicon[0].Channel;

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