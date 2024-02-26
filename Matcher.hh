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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRandom.h"

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

    chi2 = Ref_Hist->Chi2Test(TreeHist, "WW CHI2/NDF");
    graph->SetPoint(counter, par[0], chi2);
    counter++;
    return chi2;
}