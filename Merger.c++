#include "Merger.hh"
#include <ctime>

int main()
{

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
    
    for (auto &files : fileMap)
    {
        Matched_File = new TFile((dirNameMatched + files.first+"_matched.root").c_str(), "READ");
        SetTime(Matched_File);
    }

    InitHistograms();
    InitRunHistograms();

    for (auto &files : fileMap)
    {
        cout<<"<SAM> Calibration : "<<files.first<<endl;
        Matched_File = new TFile((dirNameMatched + files.first+"_matched.root").c_str(), "READ");
        Tree_Read = (TTree *)Matched_File->Get("Tree_LowHighMatched");
        TTreeReader *Reader = new TTreeReader(Tree_Read);
        TTreeReaderArray<Signal> Tree_Silicon(*Reader, "Tree_Silicon");
        TTreeReaderArray<Signal> Tree_SiPM(*Reader, "Tree_SiPM");

        InitMatching();
        
        int TotalEntries = Tree_Read->GetEntries();
        while (Reader->Next())
        {
            ULong64_t cEntry = Reader->GetCurrentEntry();
            if (cEntry % 100 == 0 && cEntry > 100)
            {
                ProgressBar(cEntry, TotalEntries, start, Current);
            }

            HSilicon_Channel_Unmatched[Tree_Silicon[0].Label]->Fill(Tree_Silicon[0].Channel);
            HSilicon_Channel_Matched[Tree_Silicon[0].Label]->Fill(detectorMatching[Tree_Silicon[0].Label] * Tree_Silicon[0].Channel);

            WriteRunHistograms(Matched_File, Tree_Silicon[0].Label, Tree_Silicon[0].Channel, detectorMatching[Tree_Silicon[0].Label] * Tree_Silicon[0].Channel);
        }

        
        Matched_File->Close();
    }
    
    InitCalib();         /////// find detector calib from simulation
    MakeCalibration();   /////// compute calibratio coeficient 

    for (auto &files : fileMap)
    {
        cout<<"<SAM> Merging : "<<files.first<<endl;
        Matched_File = new TFile((dirNameMatched + files.first + "_matched.root").c_str(), "READ");
        Tree_Read = (TTree *)Matched_File->Get("Tree_LowHighMatched");
        TTreeReader *Reader = new TTreeReader(Tree_Read);
        TTreeReaderArray<Signal> Tree_Silicon(*Reader, "Tree_Silicon");
        TTreeReaderArray<Signal> Tree_SiPM(*Reader, "Tree_SiPM");

        InitMatching();
        int TotalEntries = Tree_Read->GetEntries();
        while (Reader->Next())
        {
            ULong64_t cEntry = Reader->GetCurrentEntry();
            if (cEntry % 10000 == 0 && cEntry > 10000)
            {
                ProgressBar(cEntry, TotalEntries, start, Current);
            }

            double Silicon_Energy = detectorCalib[Tree_Silicon[0].Label] * detectorMatching[Tree_Silicon[0].Label] * Tree_Silicon[0].Channel;

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
                for (int j=0; j< Tree_SiPM.GetSize(); j++)
                {
                    sum += Tree_SiPM[j].Channel;
                }

                for (int i=1; i<Tree_SiPM.GetSize()+1; i++)
                {   
                    HSiPM[i]->Fill(sum);
                }

                if (Is_F(Silicon_Energy, Tree_Silicon[0].Label) )
                {
                    for (int i=1; i<Tree_SiPM.GetSize()+1; i++)
                    {   
                        HSiPM_F[i]->Fill(sum);
                    }
                }

                else if (Is_GT(Silicon_Energy, Tree_Silicon[0].Label) )
                {
                    for (int i=1; i<Tree_SiPM.GetSize()+1; i++)
                    {   
                        HSiPM_GT[i]->Fill(sum);
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