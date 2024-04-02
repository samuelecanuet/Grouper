#ifndef MATCHER_HH
#define MATCHER_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <algorithm> 

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRandom.h"
#include  "TKey.h"

#include "Signal.h"
#include "Detectors.hh"
#include "/home/local1/Documents/lib/GTools1.0/include/GString.hh"
using namespace std;

int RunRef_1234 = 16;
int RunRef_5678 = 30;
int RunCal = 34;   
TFile *Ref_File;
TH1D *Ref_Hist;

TTreeReader *Reader;
TTree *Tree;
TGraph* graph;
int counter;
double chi2;
int range_low;
int range_high;
int guess_high = 0.9;
int guess_low = 1.1;

int Verbose = 0;

bool Calibration = false;
bool Peaks = false;
int current_SiPM;


double SiliconMatching[SIGNAL_MAX];
double SiPMLowHighMatching[SIGNAL_MAX];
double SiPMsMatching[SIGNAL_MAX];

TTree* Tree_Read;
TTree* Tree_Write; 
TTree* Tree_Final;

TH1D* HStrip_Channel[SIGNAL_MAX];
TH1D* HRear_Channel[SIGNAL_MAX];
TH1D* HSiPMHigh_Channel[BETA_SIZE+1];
TH1D* HSiPMHigh_Channel_all[BETA_SIZE+1];
TH1D* HSiPMLow_Channel[BETA_SIZE+1];
TH1D* HSiPMLow_Channel_all[BETA_SIZE+1];
TH1D* HSiPM_Channel[BETA_SIZE+1];
TH1D* HSiPM;

TH1D* HSiPMHigh_False;
TH1D* HSiPMLow_False;
TH1D* HSiPMHigh_Time_False;
TH1D* HSiPMLow_Time_False;
TH1D* HSiPMHigh_SiPM_False;
TH1D* HSiPMLow_SiPM_False;

TH1D* HSiPMHigh_New;
TH1D* HSiPMLow_New;
TH1D* HSiPMHigh_Time_New;
TH1D* HSiPMLow_Time_New;
TH1D* HSiPMHigh_SiPM_New;
TH1D* HSiPMLow_SiPM_New;



TDirectory *dir_Chi2SiPMLowHigh;
TDirectory *dir_Chi2SiPMs;
TDirectory *dir_SiPMsReconstruction;


////////////////////////////////////

void ProgressBar(ULong64_t cEntry, ULong64_t TotalEntries, clock_t start, clock_t Current)
{
    if (cEntry % 100000 == 0 && cEntry > 2*100000)
        {
            Current = clock();
            const Char_t *Color;
            Double_t Frac = 1.0 * cEntry / TotalEntries;
            Double_t Timeclock = ((double)(Current - start) / CLOCKS_PER_SEC);
            Double_t TimeLeft = Timeclock * (1 / Frac - 1.);
            Color = "\e[1;31m";

            cout << Form("\r%s <SAM> Entry : ")
                 << TotalEntries
                 << " --- "
                 << Form("%4.2f", 100. * cEntry / TotalEntries) << " %"
                 << " --- "
                 << " Time Left : " << Form("%2d min ", (int)TimeLeft / 60)
                 << Form("%02d sec", (int)TimeLeft % 60)
                 << "           "<<flush;
        }
}



inline double Chi2TreeHist(const Double_t *par)
{
    chi2 = 0;
    Reader->Restart();
    TTreeReaderValue<double> Channel(*Reader, "Channel");
    TH1D *TreeHist = (TH1D *)Ref_Hist->Clone((to_string(par[0])).c_str());
    Ref_Hist->GetXaxis()->SetRangeUser(range_low, range_high);
    TreeHist->GetXaxis()->SetRangeUser(range_low, range_high);
    TreeHist->Reset();
    
    while (Reader->Next())
    {
        TreeHist->Fill((*Channel) * par[0]);
    }

    if (TreeHist->GetEntries() < 100 || Ref_Hist->GetEntries() < 100)
        return 0;

    chi2 = Ref_Hist->Chi2Test(TreeHist, "UU CHI2/NDF");
    graph->SetPoint(counter, par[0], chi2);
    counter++;
    return chi2;
}

inline double Chi2SiPM(const Double_t *par)
{
    chi2 = 0;
    Reader->Restart();
    TTreeReaderArray<Signal> *SiPM = new TTreeReaderArray<Signal>(*Reader, "Tree_SiPM");
    TH1D *TreeHist = (TH1D *)Ref_Hist->Clone((to_string(par[0])).c_str());
    Ref_Hist->GetXaxis()->SetRangeUser(range_low, range_high);
    TreeHist->GetXaxis()->SetRangeUser(range_low, range_high);
    TreeHist->Reset();
    
    while (Reader->Next())
    {
        for (size_t i = 0; i < SiPM->GetSize(); i++)
        {
            if (current_SiPM == GetDetectorChannel((*SiPM)[i].Label))
                TreeHist->Fill((*SiPM)[i].Channel * par[0]);
        }
    }

    if (TreeHist->GetEntries() == 0)
        return 0;

    chi2 = Ref_Hist->Chi2Test(TreeHist, "WW CHI2/NDF");
    graph->SetPoint(counter, par[0], chi2);
    counter++;
    return chi2;
}

inline int InitHistograms()
{
    
    for (size_t i = 0; i < detectorNum; ++i)
    {
        if (IsDetectorSiliStrip(i))
        {
            HStrip_Channel[i] = new TH1D(("HStri_Channel_" + detectorName[i]).c_str(), ("HStri_Channel_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
            HStrip_Channel[i]->GetXaxis()->SetTitle("Strip [Channel]");
            HStrip_Channel[i]->GetYaxis()->SetTitle("Counts");
            HStrip_Channel[i]->GetXaxis()->CenterTitle();
            HStrip_Channel[i]->GetYaxis()->CenterTitle();
        }
        if (IsDetectorSiliBack(i))
        {
            HRear_Channel[i] = new TH1D(("HRea_Channel_" + detectorName[i]).c_str(), ("HRea_Channel_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
            HRear_Channel[i]->GetXaxis()->SetTitle("Rear [Channel]");
            HRear_Channel[i]->GetYaxis()->SetTitle("Counts");
            HRear_Channel[i]->GetXaxis()->CenterTitle();
            HRear_Channel[i]->GetYaxis()->CenterTitle();
        }
    }

    for (int i = 1; i <= BETA_SIZE; i++)
    {
            HSiPMHigh_Channel[i] = new TH1D(("HSiPMHigh_" + to_string(i) + "_Channel").c_str(), ("HSiPMHigh_" + to_string(i) + "hannel").c_str(), eLowN/10, eLowMin, eLowMax);
            HSiPMHigh_Channel[i]->GetXaxis()->SetTitle("SiPM [Channel]");
            HSiPMHigh_Channel[i]->GetYaxis()->SetTitle("Counts");
            HSiPMHigh_Channel[i]->GetXaxis()->CenterTitle();
            HSiPMHigh_Channel[i]->GetYaxis()->CenterTitle();

            HSiPMHigh_Channel_all[i] = new TH1D(("HSiPMHigh_all_" + to_string(i) + "_Channel").c_str(), ("HSiPMHigh_all_" + to_string(i) + "hannel").c_str(), eLowN/10, eLowMin, eLowMax);
            HSiPMHigh_Channel_all[i]->GetXaxis()->SetTitle("SiPM [Channel]");
            HSiPMHigh_Channel_all[i]->GetYaxis()->SetTitle("Counts");
            HSiPMHigh_Channel_all[i]->GetXaxis()->CenterTitle();
            HSiPMHigh_Channel_all[i]->GetYaxis()->CenterTitle();

            HSiPMLow_Channel[i] = new TH1D(("HSiPMLow_" + to_string(i) + "_Channel").c_str(), ("HSiPMLow_" + to_string(i) + "hannel").c_str(), eLowN/10, eLowMin, eLowMax);
            HSiPMLow_Channel[i]->GetXaxis()->SetTitle("SiPM [Channel]");
            HSiPMLow_Channel[i]->GetYaxis()->SetTitle("Counts");
            HSiPMLow_Channel[i]->GetXaxis()->CenterTitle();
            HSiPMLow_Channel[i]->GetYaxis()->CenterTitle();

            HSiPMLow_Channel_all[i] = new TH1D(("HSiPMLow_all_" + to_string(i) + "_Channel").c_str(), ("HSiPMLow_all_" + to_string(i) + "hannel").c_str(), eLowN/10, eLowMin, eLowMax);
            HSiPMLow_Channel_all[i]->GetXaxis()->SetTitle("SiPM [Channel]");
            HSiPMLow_Channel_all[i]->GetYaxis()->SetTitle("Counts");
            HSiPMLow_Channel_all[i]->GetXaxis()->CenterTitle();
            HSiPMLow_Channel_all[i]->GetYaxis()->CenterTitle();

            HSiPM_Channel[i] = new TH1D(("HSiPM_" + to_string(i) + "_Channel").c_str(), ("HSiPM_" + to_string(i) + "hannel").c_str(), eLowN/10, eLowMin, eLowMax);
            HSiPM_Channel[i]->GetXaxis()->SetTitle("SiPM [Channel]");
            HSiPM_Channel[i]->GetYaxis()->SetTitle("Counts");
            HSiPM_Channel[i]->GetXaxis()->CenterTitle();
            HSiPM_Channel[i]->GetYaxis()->CenterTitle();


    }

    HSiPM = new TH1D("HSiPM", "HSiPM", eLowN/5, eLowMin, eLowMax*5);
    HSiPM->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPM->GetYaxis()->SetTitle("Counts");
    HSiPM->GetXaxis()->CenterTitle();
    HSiPM->GetYaxis()->CenterTitle();
    
    HSiPMHigh_False = new TH1D("HSiPMHigh_False", "HSiPMHigh_False", eLowN/10, eLowMin, eLowMax);
    HSiPMHigh_False->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMHigh_False->GetYaxis()->SetTitle("Counts");
    HSiPMHigh_False->GetXaxis()->CenterTitle();
    HSiPMHigh_False->GetYaxis()->CenterTitle();

    HSiPMLow_False = new TH1D("HSiPMLow_False", "HSiPMLow_False", eLowN/10, eLowMin, eLowMax);
    HSiPMLow_False->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMLow_False->GetYaxis()->SetTitle("Counts");
    HSiPMLow_False->GetXaxis()->CenterTitle();
    HSiPMLow_False->GetYaxis()->CenterTitle();

    HSiPMHigh_Time_False = new TH1D("HSiPMHigh_Time_False", "HSiPMHigh_Time_False", winHighN, winHighMin, winHighMax);
    HSiPMHigh_Time_False->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMHigh_Time_False->GetYaxis()->SetTitle("Counts");
    HSiPMHigh_Time_False->GetXaxis()->CenterTitle();
    HSiPMHigh_Time_False->GetYaxis()->CenterTitle();

    HSiPMLow_Time_False = new TH1D("HSiPMLow_Time_False", "HSiPMLow_Time_False", winLowN, winLowMin, winLowMax);
    HSiPMLow_Time_False->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMLow_Time_False->GetYaxis()->SetTitle("Counts");
    HSiPMLow_Time_False->GetXaxis()->CenterTitle();
    HSiPMLow_Time_False->GetYaxis()->CenterTitle();

    HSiPMHigh_SiPM_False = new TH1D("HSiPMHigh_SiPM_False", "HSiPMHigh_SiPM_False", BETA_SIZE, 1, BETA_SIZE+1);
    HSiPMHigh_SiPM_False->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMHigh_SiPM_False->GetYaxis()->SetTitle("Counts");
    HSiPMHigh_SiPM_False->GetXaxis()->CenterTitle();
    HSiPMHigh_SiPM_False->GetYaxis()->CenterTitle();

    HSiPMLow_SiPM_False = new TH1D("HSiPMLow_SiPM_False", "HSiPMLow_SiPM_False", BETA_SIZE, 1, BETA_SIZE+1);
    HSiPMLow_SiPM_False->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMLow_SiPM_False->GetYaxis()->SetTitle("Counts");
    HSiPMLow_SiPM_False->GetXaxis()->CenterTitle();
    HSiPMLow_SiPM_False->GetYaxis()->CenterTitle();

    HSiPMHigh_New = new TH1D("HSiPMHigh_New", "HSiPMHigh_New", 500, eLowMin, eLowMax);
    HSiPMHigh_New->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMHigh_New->GetYaxis()->SetTitle("Counts");
    HSiPMHigh_New->GetXaxis()->CenterTitle();
    HSiPMHigh_New->GetYaxis()->CenterTitle();
    
    HSiPMLow_New = new TH1D("HSiPMLow_New", "HSiPMLow_New", 500, eLowMin, eLowMax);
    HSiPMLow_New->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMLow_New->GetYaxis()->SetTitle("Counts");
    HSiPMLow_New->GetXaxis()->CenterTitle();
    HSiPMLow_New->GetYaxis()->CenterTitle();

    HSiPMHigh_Time_New = new TH1D("HSiPMHigh_Time_New", "HSiPMHigh_Time_New", winHighN, winHighMin, winHighMax);
    HSiPMHigh_Time_New->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMHigh_Time_New->GetYaxis()->SetTitle("Counts");
    HSiPMHigh_Time_New->GetXaxis()->CenterTitle();
    HSiPMHigh_Time_New->GetYaxis()->CenterTitle();

    HSiPMLow_Time_New = new TH1D("HSiPMLow_Time_New", "HSiPMLow_Time_New", winLowN, winLowMin, winLowMax);
    HSiPMLow_Time_New->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMLow_Time_New->GetYaxis()->SetTitle("Counts");
    HSiPMLow_Time_New->GetXaxis()->CenterTitle();
    HSiPMLow_Time_New->GetYaxis()->CenterTitle();

    HSiPMHigh_SiPM_New = new TH1D("HSiPMHigh_SiPM_New", "HSiPMHigh_SiPM_New", BETA_SIZE, 1, BETA_SIZE+1);
    HSiPMHigh_SiPM_New->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMHigh_SiPM_New->GetYaxis()->SetTitle("Counts");
    HSiPMHigh_SiPM_New->GetXaxis()->CenterTitle();
    HSiPMHigh_SiPM_New->GetYaxis()->CenterTitle();

    HSiPMLow_SiPM_New = new TH1D("HSiPMLow_SiPM_New", "HSiPMLow_SiPM_New", BETA_SIZE, 1, BETA_SIZE+1);
    HSiPMLow_SiPM_New->GetXaxis()->SetTitle("SiPM [Channel]");
    HSiPMLow_SiPM_New->GetYaxis()->SetTitle("Counts");
    HSiPMLow_SiPM_New->GetXaxis()->CenterTitle();
    HSiPMLow_SiPM_New->GetYaxis()->CenterTitle();


    return 0;
}

inline int WriteHistograms()
{
    for (size_t i = 0; i < detectorNum; ++i)
    {
        // if (IsDetectorSiliStrip(i))
        // {
        //     HStri_Channel[i]->Write();
        // }
        // if (IsDetectorSiliBack(i))
        // {
        //     HRea_Channel[i]->Write();
        // }
    }

    for (int i = 1; i <= BETA_SIZE; i++)
    {
        dir_Chi2SiPMLowHigh->cd();
        TCanvas *c = new TCanvas(("HSiPM_Merged" + to_string(i) + "_Channel").c_str(), ("HSiPM_" + to_string(i) + "_Channel").c_str(), 800, 600);
        c->cd();
        HSiPMHigh_Channel[i]->SetLineColor(kRed);
        HSiPMHigh_Channel[i]->Draw();
        HSiPMLow_Channel[i]->SetLineColor(kBlue);
        HSiPMLow_Channel[i]->Draw("SAME");
        HSiPM_Channel[i]->SetLineColor(kBlack);
        HSiPM_Channel[i]->Draw("SAME");
        c->Write();
        delete c;

        dir_Chi2SiPMs->cd();
        TCanvas *c1 = new TCanvas(("HSiPMHigh_1vs" + to_string(i)).c_str(), ("HSiPMHigh_1vs" + to_string(i)).c_str(), 800, 600);
        c1->cd();
        HSiPMHigh_Channel_all[i]->SetLineColor(kBlue);
        HSiPMHigh_Channel_all[i]->Draw("HIST");
        HSiPMHigh_Channel_all[1]->SetLineColor(kBlack);
        HSiPMHigh_Channel_all[1]->Draw("SAME");
        HSiPMHigh_Channel_all[i]->SetTitle(("HSiPMHigh_1vs" + to_string(i)).c_str());
        HSiPMHigh_Channel_all[1]->SetTitle(("HSiPMHigh_1vs" + to_string(i)).c_str());
        c1->Write();
        delete c1;

        TCanvas *c2 = new TCanvas(("HSiPMLow_1vs" + to_string(i)).c_str(), ("HSiPMLow_1vs" + to_string(i)).c_str(), 800, 600);
        c2->cd();
        HSiPMLow_Channel_all[i]->SetLineColor(kBlue);
        HSiPMLow_Channel_all[i]->Draw("HIST");
        HSiPMLow_Channel_all[1]->SetLineColor(kBlack);
        HSiPMLow_Channel_all[1]->Draw("SAME");
        HSiPMLow_Channel_all[i]->SetTitle(("HSiPMLow_1vs" + to_string(i)).c_str());
        HSiPMLow_Channel_all[1]->SetTitle(("HSiPMLow_1vs" + to_string(i)).c_str());
        c2->Write();
        delete c2;
    }
        dir_SiPMsReconstruction->cd();
        HSiPMHigh_False->Write();
        HSiPMLow_False->Write();
        HSiPMHigh_Time_False->Write();
        HSiPMLow_Time_False->Write();
        HSiPMHigh_SiPM_False->Write();
        HSiPMLow_SiPM_False->Write();
        
        HSiPMHigh_New->Write();
        HSiPMLow_New->Write();
        HSiPMHigh_Time_New->Write();
        HSiPMLow_Time_New->Write();
        HSiPMHigh_SiPM_New->Write();
        HSiPMLow_SiPM_New->Write();

    

    HSiPM->Write();
    return 0;
}

#endif