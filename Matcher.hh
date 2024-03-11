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
TH1I *Ref_Hist;

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


double detectorMatching[SIGNAL_MAX];

TTree* Tree_Read;
TTree* Tree_Write; 

TH1I* HStrip_Channel[SIGNAL_MAX];
TH1I* HRear_Channel[SIGNAL_MAX];
TH1I* HSiPMHigh_Channel[BETA_SIZE+1];
TH1I* HSiPMLow_Channel[BETA_SIZE+1];
TH1I* HSiPM_Channel[BETA_SIZE+1];

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
    TTreeReaderValue<int> Channel(*Reader, "Channel");
    TH1I *TreeHist = (TH1I *)Ref_Hist->Clone((to_string(par[0])).c_str());
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
    TH1I *TreeHist = (TH1I *)Ref_Hist->Clone((to_string(par[0])).c_str());
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

    chi2 = Ref_Hist->Chi2Test(TreeHist, "UU CHI2/NDF");
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
            HStrip_Channel[i] = new TH1I(("HStrip_Channel_" + detectorName[i]).c_str(), ("HStrip_Channel_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
            HStrip_Channel[i]->GetXaxis()->SetTitle("Strip [Channel]");
            HStrip_Channel[i]->GetYaxis()->SetTitle("Counts");
            HStrip_Channel[i]->GetXaxis()->CenterTitle();
            HStrip_Channel[i]->GetYaxis()->CenterTitle();
        }
        if (IsDetectorSiliBack(i))
        {
            HRear_Channel[i] = new TH1I(("HRear_Channel_" + detectorName[i]).c_str(), ("HRear_Channel_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
            HRear_Channel[i]->GetXaxis()->SetTitle("Rear [Channel]");
            HRear_Channel[i]->GetYaxis()->SetTitle("Counts");
            HRear_Channel[i]->GetXaxis()->CenterTitle();
            HRear_Channel[i]->GetYaxis()->CenterTitle();
        }
    }

    for (int i = 1; i <= BETA_SIZE; i++)
    {
            HSiPMHigh_Channel[i] = new TH1I(("HSiPMHigh_" + to_string(i) + "_Channel").c_str(), ("HSiPMHigh_" + to_string(i) + "_Channel").c_str(), eLowN/10, eLowMin, eLowMax);
            HSiPMHigh_Channel[i]->GetXaxis()->SetTitle("SiPM [Channel]");
            HSiPMHigh_Channel[i]->GetYaxis()->SetTitle("Counts");
            HSiPMHigh_Channel[i]->GetXaxis()->CenterTitle();
            HSiPMHigh_Channel[i]->GetYaxis()->CenterTitle();

            HSiPMLow_Channel[i] = new TH1I(("HSiPMLow_" + to_string(i) + "_Channel").c_str(), ("HSiPMLow_" + to_string(i) + "_Channel").c_str(), eLowN/10, eLowMin, eLowMax);
            HSiPMLow_Channel[i]->GetXaxis()->SetTitle("SiPM [Channel]");
            HSiPMLow_Channel[i]->GetYaxis()->SetTitle("Counts");
            HSiPMLow_Channel[i]->GetXaxis()->CenterTitle();
            HSiPMLow_Channel[i]->GetYaxis()->CenterTitle();

            HSiPM_Channel[i] = new TH1I(("HSiPM_" + to_string(i) + "_Channel").c_str(), ("HSiPM_" + to_string(i) + "_Channel").c_str(), eLowN/10, eLowMin, eLowMax);
            HSiPM_Channel[i]->GetXaxis()->SetTitle("SiPM [Channel]");
            HSiPM_Channel[i]->GetYaxis()->SetTitle("Counts");
            HSiPM_Channel[i]->GetXaxis()->CenterTitle();
            HSiPM_Channel[i]->GetYaxis()->CenterTitle();


    }
    return 0;
}

inline int WriteHistograms()
{
    for (size_t i = 0; i < detectorNum; ++i)
    {
        // if (IsDetectorSiliStrip(i))
        // {
        //     HStrip_Channel[i]->Write();
        // }
        // if (IsDetectorSiliBack(i))
        // {
        //     HRear_Channel[i]->Write();
        // }
    }

    for (int i = 1; i <= BETA_SIZE; i++)
    {
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

    }
    return 0;
}
