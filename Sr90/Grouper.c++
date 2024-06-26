#include "Grouper.hh"
#include <ctime>


int main(int argc, char *argv[])
{
    clock_t start = clock(), Current;

    ///////////////////////////////////  INPUT ///////////////////////////////////
    baseFileName = "run_069_90Sr";
    dirNameGrouped = "../../../../../../../mnt/hgfs/shared-2/";
    string ROOTFileName = "../../../../../../../mnt/hgfs/shared-2/" + baseFileName + ".root";
    TFile *Data_File = new TFile(ROOTFileName.c_str(), "READ");
    TTree *Tree = (TTree *)Data_File->Get("Tree");
    TTreeReader *Reader = new TTreeReader(Tree);

    TTreeReaderValue<double> Channel(*Reader, "Channel");
    TTreeReaderValue<int> Label(*Reader, "Label");
    TTreeReaderValue<double> Time(*Reader, "Time");

    ///////////////////////////////////  DELETE OLD FILE ///////////////////////////////////
    int status = remove((dirNameGrouped+baseFileName+"_grouped.root").c_str());
    status = remove((dirNameCleaned+baseFileName+"_cleaned.root").c_str());

    ///////////////////////////////////  Grouped ///////////////////////////////////
    File_Grouped = new TFile((dirNameGrouped+baseFileName+"_grouped.root").c_str(), "RECREATE");

    ///////////////////////////////////  INITIALISATION ///////////////////////////////////
    InitDetectors("../../Analysis/ReadFaster/Data/detectors.dat");

    InitHistograms_Grouped();
    InitTree_Grouped();

    ///////////////////////////////////  PROCESSING TREE ///////////////////////////////////
    

    Verbose = 0;
    ULong64_t TotalEntries = Tree->GetEntries();
    GLogMessage("<SAM> Creating Groups : ");
    while (Reader->Next())
    {
        ULong64_t cEntry = Reader->GetCurrentEntry();
        if (cEntry%100000 == 0 && cEntry > 2000000)
        {
            ProgressBar(cEntry, TotalEntries, start, Current, "");
        }

        if (Verbose > 1)
            cout << "DATA n°" << Reader->GetCurrentEntry() << " ---> "
                 << "\t" << detectorName[*Label] << "\t" << setprecision(15) << *Time << "\t" << *Channel << endl;


        SearchForCoincidence(Reader, Time, Label, Channel);
        
    
    }
    cout << "\n" <<"\e[0m" << flush ;
    cout << "Writing Data" << endl;
    WriteHistograms_Grouped();
    WriteTree_Grouped();
    WriteTime(Data_File, File_Grouped);
    File_Grouped->Close();
    Data_File->Close();


    // ///////////////////////////////////  CLEANING ///////////////////////////////////
    // File_Grouped = new TFile((dirNameGrouped+baseFileName+"_grouped.root").c_str(), "READ");
    // File_Cleaned = new TFile((dirNameCleaned+baseFileName+"_cleaned.root").c_str(), "RECREATE");
    // WriteTime(File_Grouped, File_Cleaned);
    // Tree_Grouped = (TTree *)File_Grouped->Get("Tree");
    // GLogMessage("<SAM> Cleaning Proton Spectrum : ");
    // InitCleaning();
    // TTreeReader *Grouped_Reader = new TTreeReader(Tree_Grouped);
    // TTreeReaderArray<Signal> Grouped_Silicon(*Grouped_Reader, "Tree_Silicon");
    // TTreeReaderArray<Signal> Grouped_SiPMHigh(*Grouped_Reader, "Tree_SiPMHigh");
    // TTreeReaderArray<Signal> Grouped_SiPMLow(*Grouped_Reader, "Tree_SiPMLow");

    // InitHistograms_Cleaned();
    // InitTree_Cleaned();

    // Tree_Grouped->SetBranchAddress("Tree_Silicon", &Tree_Cleaned_Silicon);
    // Tree_Grouped->SetBranchAddress("Tree_SiPMHigh", &Tree_Cleaned_SiPMHigh);
    // Tree_Grouped->SetBranchAddress("Tree_SiPMLow", &Tree_Cleaned_SiPMLow);

    // TotalEntries = Tree_Grouped->GetEntries();
    // Event = 0;
    // double sigma_acceptance = 3;
    // while(Grouped_Reader->Next())
    // {
    //     ULong64_t cEntry = Grouped_Reader->GetCurrentEntry();
    //     if (cEntry%10000 == 0 && cEntry > 10000)
    //     {
    //         ProgressBar(cEntry, TotalEntries, start, Current, "");
    //     }
        
    //     double min = detectorCleaning[Grouped_Silicon[0].Label].first - sigma_acceptance*detectorCleaning[Grouped_Silicon[0].Label].second;
    //     double max = detectorCleaning[Grouped_Silicon[0].Label].first + sigma_acceptance*detectorCleaning[Grouped_Silicon[0].Label].second;
    //     double frac = static_cast<double>(Grouped_Silicon[0].Channel)/Grouped_Silicon[1].Channel;
    //     if ( frac > min && frac < max)
    //     {   
    //         Tree_Grouped->GetEntry(cEntry);
    //         SavingData_Cleaned(*Tree_Cleaned_Silicon, *Tree_Cleaned_SiPMHigh, *Tree_Cleaned_SiPMLow);


    //         ////FOR STRIP
    //         Tree_Channel = Grouped_Silicon[0].Channel;
    //         Tree_Silicons[Grouped_Silicon[0].Label]->Fill();
    //         ////FOR REAR
    //         Tree_Channel = Grouped_Silicon[1].Channel;
    //         Tree_Silicons[Grouped_Silicon[1].Label]->Fill();

    //         Tree_Cleaned->Fill();
    //     }
    // }
    // cout << "\n" <<"\e[0m" << flush ;

    // WriteHistograms_Cleaned();
    // WriteTree_Cleaned();
    // File_Cleaned->Close();
    // File_Grouped->Close();
}
