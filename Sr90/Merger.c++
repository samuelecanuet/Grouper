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
    string dirNameMerged = "./../../../../../../../mnt/hgfs/shared-2/";
    string dirNameCleaned = "./../../../../../../../mnt/hgfs/shared-2/";
    string dirNameMatched = "./../../../../../../../mnt/hgfs/shared-2/";

    string filename = "run_069_90Sr";
    Merged_File = new TFile((dirNameMerged + filename + "_merged.root").c_str(), "RECREATE");
    Tree_Merged = new TTree("Tree", "Tree");
    Tree_Merged->Branch("Tree_Silicon", &Tree_Merged_Silicon);
    Tree_Merged->Branch("Tree_SiPM", &Tree_Merged_SiPM);

    InitDetectors("../../Analysis/ReadFaster/Data/detectors.dat");

    Matched_File = new TFile((dirNameMatched + filename + "_matched.root").c_str(), "READ");
    SetTime(Matched_File);
    

    ///////////////////// CREATE TIME RUNs FILES & BUILD PROTON HISTOGRAM FOR CALIBRATION ///////////////
    InitHistograms();
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    InitCalibSiPM();          /////// find SiPM calib from simulation

    ////////SiPM CALIBRATION////////
    for (int i = 0; i < BETA_SIZE+1; i++)
    {
        for (int j = 0; j < BETA_SIZE+1; j++)
        {
            Ref_Hist_F[j][i] = new TH1D(("Fermi_M"+ to_string(j)+"_SiPM" + to_string(i)).c_str(), ("Fermi_M"+ to_string(j)+"_SiPM" + to_string(i)).c_str(), 1200, 0, 12000);
        }
        Ref_Hist_Multi_F[i] = new TH1D(("Fermi_M" + to_string(i)).c_str(), ("Fermi_M" + to_string(i)).c_str(), 1200, 0, 12000);
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
    }


        Matched_File = new TFile((dirNameMatched + filename + "_matched.root").c_str(), "READ");
        Tree_Read = (TTree *)Matched_File->Get("Tree");
        TTreeReader *Reader = new TTreeReader(Tree_Read);
        TTreeReaderArray<Signal> Tree_SiPM(*Reader, "Tree_SiPM");

        vector<pair<int, double>> vec;
        vector<int> result;
        int TotalEntries = Tree_Read->GetEntries();

        while (Reader->Next() && Reader->GetCurrentEntry() < 1000000)
        {
            ULong64_t cEntry = Reader->GetCurrentEntry();
            ProgressBar(cEntry, TotalEntries, start, Current, "");
            

            // if (Is_F(1 / SiliconCalib[99][Tree_Silicon[0].Label].first * Tree_Silicon[0].Channel, Tree_Silicon[0].Label))
            // {
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
            // }

            if (Tree_SiPM.GetSize() == 1)
            {
                histGT->Fill(Tree_SiPM[0].Channel*4.9/1000);
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
        // pair<string, string> string_time = GetTime(Matched_File);
        // start_time = Convert_DatetoTime(string_time.first, file_string_Time[0].first);
        // stop_time = Convert_DatetoTime(string_time.second, file_string_Time[0].first);

        // bool CurrentDetectorSelection[SIGNAL_MAX];
        // copy(begin(RunDetectorSelection[run]), end(RunDetectorSelection[run]), begin(CurrentDetectorSelection));

        Matched_File->Close();
    
    Merged_File->cd();
    // MakeSiPMCalibration(0); // convoluted + th but no multiplicity convoluting
    // MakeSiPM_SiPMPlots();
    MakeSiPM_MultiplicityPlots();
    Merged_File->cd();
    histGT->Write();
    //MakeSiPM_ParameterPlots();

    ////////MERGING////////
    // for (auto &run : runs)
    // {
    //     cout << "<SAM> Merging on run : " << run << endl;
    //     Matched_File = new TFile((dirNameMatched + "run_0" + to_string(run) + "_32Ar_matched.root").c_str(), "READ");
    //     Tree_Read = (TTree *)Matched_File->Get("Tree_LowHighMatched");
    //     TTreeReader *Reader = new TTreeReader(Tree_Read);
    //     TTreeReaderArray<Signal> Tree_Silicon(*Reader, "Tree_Silicon");
    //     TTreeReaderArray<Signal> Tree_SiPM(*Reader, "Tree_SiPM");

    //     pair<string, string> string_time = GetTime(Matched_File);
    //     start_time = Convert_DatetoTime(string_time.first, file_string_Time[0].first);
    //     stop_time = Convert_DatetoTime(string_time.second, file_string_Time[0].first);

    //     InitMatching();
    //     bool CurrentDetectorSelection[SIGNAL_MAX];
    //     copy(begin(RunDetectorSelection[run]), end(RunDetectorSelection[run]), begin(CurrentDetectorSelection));

    //     int TotalEntries = Tree_Read->GetEntries();
    //     while (Reader->Next())
    //     {
    //         // ULong64_t cEntry = Reader->GetCurrentEntry();
    //         // if (cEntry % 10000 == 0 && cEntry > 10000)
    //         // {
    //         //     ProgressBar(cEntry, TotalEntries, start, Current);
    //         // }

    //         if (CurrentDetectorSelection[Tree_Silicon[0].Label])
    //         {
    //             double Silicon_Energy = 1 / SiliconCalib[99][Tree_Silicon[0].Label].first * Tree_Silicon[0].Channel;

    //             HSiliconRun[Tree_Silicon[0].Label]->Fill(start_time + (double)Tree_Silicon[0].Time * 1e-9 / 3600, Silicon_Energy);
    //             HSilicon[Tree_Silicon[0].Label]->Fill(Silicon_Energy);

    //             if (Tree_SiPM.GetSize() != 0)
    //             {
    //                 HSilicon_coinc[Tree_Silicon[0].Label]->Fill(Silicon_Energy);
    //             }
    //             else
    //             {
    //                 HSilicon_no_coinc[Tree_Silicon[0].Label]->Fill(Silicon_Energy);
    //             }
    //         }
    //     }
    //     Matched_File->Close();
    // }

    // WriteHistograms();
    Merged_File->Close();

    return 0;
}