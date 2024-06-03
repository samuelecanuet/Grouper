#ifndef GROUPER_HH
#define GROUPER_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream> // Pour l'écriture dans un fichier

#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TTree.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRandom.h"

#include "Signal.h"
#include "Detectors.hh"

#include "/home/local1/Documents/lib/GTools1.0/include/GString.hh"

using namespace std;

int Verbose = 0;
int RunCal;

pair<double, double> detectorCleaning[SIGNAL_MAX];

double Event = 0;
int counter_graph[BETA_SIZE + 1];

string baseFileName;
string dirNameGrouped;
string dirNameCleaned;
string dirNameMatched;

TFile *File_Grouped;
TTree *Tree_Grouped;
int Tree_Event;
vector<Signal> Tree_Silicon;
vector<Signal> Tree_SiPMHigh;
vector<Signal> Tree_SiPMLow;

TFile *File_Cleaned;
TTree *Tree_Cleaned;
int Tree_Cleaned_Event;
vector<Signal> *Tree_Cleaned_Silicon;
vector<Signal> *Tree_Cleaned_SiPMHigh;
vector<Signal> *Tree_Cleaned_SiPMLow;
double Tree_Channel;

//////////////GROUPED////////////////
/// SiPM High
TH1D *HSiPMHigh_Channel_G[SIGNAL_MAX];
TH1D *HSiPMHighRear_Time_G[SIGNAL_MAX];
TH2D *HSiPMHighRear_TimeChannel_G[SIGNAL_MAX];

/// SiPM Low
TH1D *HSiPMLow_Channel_G[SIGNAL_MAX];
TH1D *HSiPMLowRear_Time_G[SIGNAL_MAX];
TH2D *HSiPMLowRear_TimeChannel_G[SIGNAL_MAX];

/// SIPMs
TH1D *HSiPMHigh_False_G;
TH1D *HSiPMLow_False_G;
TH1D *HSiPMHigh_Time_False_G;
TH1D *HSiPMLow_Time_False_G;
TH1D *HSiPMHigh_SiPM_False_G;
TH1D *HSiPMLow_SiPM_False_G;

double counter_true = 0;
double counter_high = 0;
double counter_low = 0;
TH1D *HSiPM_Counter_False_G;

TH1D *HSiPMHigh_Multiplicity_G;
TH1D *HSiPMLow_Multiplicity_G;
TH2D *HSiPM_Multiplicities_G;
/////////////////////////////////////
//////////////CLEANED////////////////
/// Silicon
TH1D *HStrip_Channel_C[SIGNAL_MAX];
TH1D *HRear_Channel_C[SIGNAL_MAX];
TH1D *HStripRear_Time_C[SIGNAL_MAX];
TH2D *HStripRear_Channel_C[SIGNAL_MAX];
TH1D *HFracRearStripChannel_C[SIGNAL_MAX];

TH2D *HStripsMultiplicity_C[SIGNAL_MAX];

/// SiPM High
TH1D *HSiPMHigh_Channel_C[SIGNAL_MAX];
TH1D *HSiPMHighRear_Time_C[SIGNAL_MAX];
TH2D *HSiPMHighRear_TimeChannel_C[SIGNAL_MAX];

/// SiPM Low
TH1D *HSiPMLow_Channel_C[SIGNAL_MAX];
TH1D *HSiPMLowRear_Time_C[SIGNAL_MAX];
TH2D *HSiPMLowRear_TimeChannel_C[SIGNAL_MAX];

/// SIPMs
TH1D *HSiPMHigh_Multiplicity_C;
TH1D *HSiPMLow_Multiplicity_C;
TH2D *HSiPM_Multiplicities_C;
TGraph *GSiPM_Channel_C[BETA_SIZE + 1];
TProfile *PSiPM_Channel_C[BETA_SIZE + 1];

TTree *Tree_Silicons[SIGNAL_MAX];

TH2D *HSiPM12_Channel_C;
TH2D *HSiPM13_Channel_C;
TH2D *HSiPM23_Channel_C;

TH2D *HSiPM12LOW_Channel_C;
TH2D *HSiPM13LOW_Channel_C;
TH2D *HSiPM23LOW_Channel_C;




//////////////////////////////////////

inline int InitHistograms_Grouped()
{
  for (size_t i = 0; i < SIGNAL_MAX; ++i)
  {
    // HSiPMRearHigh_TimeChannel_G[i] = NULL;
    HSiPMHigh_Channel_G[i] = NULL;
    HSiPMHigh_Multiplicity_G = NULL;
    HSiPMHighRear_Time_G[i] = NULL;
    HSiPMHighRear_TimeChannel_G[i] = NULL;
    HSiPMLow_Channel_G[i] = NULL;
    HSiPMLow_Multiplicity_G = NULL;
    HSiPMLowRear_Time_G[i] = NULL;
    HSiPMLowRear_TimeChannel_G[i] = NULL;
    HSiPM_Multiplicities_G = NULL;

    HSiPMHigh_False_G = NULL;
    HSiPMLow_False_G = NULL;
    HSiPMHigh_Time_False_G = NULL;
    HSiPMLow_Time_False_G = NULL;
    HSiPMHigh_SiPM_False_G = NULL;
    HSiPMLow_SiPM_False_G = NULL;

    HSiPM_Counter_False_G = NULL;
  }

  for (size_t i = 0; i < detectorNum; ++i)
  {
    if (IsDetectorBetaHigh(i))
    {
      HSiPMHigh_Channel_G[i] = new TH1D(("HSiPMHigh_Channel_G_" + detectorName[i]).c_str(), ("HSiPMHigh_Channel_G_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
      HSiPMHigh_Channel_G[i]->GetXaxis()->SetTitle("SiPM High [Channel]");
      HSiPMHigh_Channel_G[i]->GetYaxis()->SetTitle("Counts");
      HSiPMHigh_Channel_G[i]->GetXaxis()->CenterTitle();
      HSiPMHigh_Channel_G[i]->GetYaxis()->CenterTitle();

      HSiPMHighRear_Time_G[i] = new TH1D(("HSiPMHighRear_Time_G_" + detectorName[i]).c_str(), ("HSiPMHighRear_Time_G_" + detectorName[i]).c_str(), winHighN, winHighMin, winHighMax);
      HSiPMHighRear_Time_G[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMHighRear_Time_G[i]->GetYaxis()->SetTitle("Counts");
      HSiPMHighRear_Time_G[i]->GetXaxis()->CenterTitle();
      HSiPMHighRear_Time_G[i]->GetYaxis()->CenterTitle();

      HSiPMHighRear_TimeChannel_G[i] = new TH2D(("HSiPMHighRear_TimeChannel_G_" + detectorName[i]).c_str(), ("HSiPMHighRear_TimeChannel_G_" + detectorName[i]).c_str(), winHighN, winHighMin, winHighMax, eHighN, eHighMin, eHighMax);
      HSiPMHighRear_TimeChannel_G[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMHighRear_TimeChannel_G[i]->GetYaxis()->SetTitle("SiPM High [Channel]");
      HSiPMHighRear_TimeChannel_G[i]->GetXaxis()->CenterTitle();
      HSiPMHighRear_TimeChannel_G[i]->GetYaxis()->CenterTitle();
      HSiPMHighRear_TimeChannel_G[i]->SetDrawOption("COLZ");
    }

    if (IsDetectorBetaLow(i))
    {
      HSiPMLow_Channel_G[i] = new TH1D(("HSiPMLow_Channel_G_" + detectorName[i]).c_str(), ("HSiPMLow_Channel_G_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax);
      HSiPMLow_Channel_G[i]->GetXaxis()->SetTitle("SiPM Low [Channel]");
      HSiPMLow_Channel_G[i]->GetYaxis()->SetTitle("Counts");
      HSiPMLow_Channel_G[i]->GetXaxis()->CenterTitle();
      HSiPMLow_Channel_G[i]->GetYaxis()->CenterTitle();

      HSiPMLowRear_Time_G[i] = new TH1D(("HSiPMLowRear_Time_G_" + detectorName[i]).c_str(), ("HSiPMLowRear_Time_G_" + detectorName[i]).c_str(), winLowN, winLowMin, winLowMax);
      HSiPMLowRear_Time_G[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMLowRear_Time_G[i]->GetYaxis()->SetTitle("Counts");
      HSiPMLowRear_Time_G[i]->GetXaxis()->CenterTitle();
      HSiPMLowRear_Time_G[i]->GetYaxis()->CenterTitle();

      HSiPMLowRear_TimeChannel_G[i] = new TH2D(("HSiPMLowRear_TimeChannel_G_" + detectorName[i]).c_str(), ("HSiPMLowRear_TimeChannel_G_" + detectorName[i]).c_str(), winLowN, winLowMin, winLowMax, eLowN, eLowMin, eLowMax);
      HSiPMLowRear_TimeChannel_G[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMLowRear_TimeChannel_G[i]->GetYaxis()->SetTitle("SiPM Low [Channel]");
      HSiPMLowRear_TimeChannel_G[i]->GetXaxis()->CenterTitle();
      HSiPMLowRear_TimeChannel_G[i]->GetYaxis()->CenterTitle();
      HSiPMLowRear_TimeChannel_G[i]->SetDrawOption("COLZ");
    }
  }

  HSiPMHigh_Multiplicity_G = new TH1D(("HSiPMHigh_Multiplicity_G"), ("HSiPMHigh_Multiplicity_G"), 10, 0, 10);
  HSiPMHigh_Multiplicity_G->GetXaxis()->SetTitle("Multiplicity");
  HSiPMHigh_Multiplicity_G->GetYaxis()->SetTitle("Counts");
  HSiPMHigh_Multiplicity_G->GetXaxis()->CenterTitle();
  HSiPMHigh_Multiplicity_G->GetYaxis()->CenterTitle();

  HSiPMLow_Multiplicity_G = new TH1D(("HSiPMLow_Multiplicity_G"), ("HSiPMLow_Multiplicity_G"), 10, 0, 10);
  HSiPMLow_Multiplicity_G->GetXaxis()->SetTitle("Multiplicity");
  HSiPMLow_Multiplicity_G->GetYaxis()->SetTitle("Counts");
  HSiPMLow_Multiplicity_G->GetXaxis()->CenterTitle();
  HSiPMLow_Multiplicity_G->GetYaxis()->CenterTitle();

  HSiPM_Multiplicities_G = new TH2D(("HSiPM_Multiplicities_G"), ("HSiPM_Multiplicities_G"), 10, 0, 10, 10, 0, 10);
  HSiPM_Multiplicities_G->GetXaxis()->SetTitle("Multiplicity High");
  HSiPM_Multiplicities_G->GetYaxis()->SetTitle("Multiplicity Low");
  HSiPM_Multiplicities_G->GetXaxis()->CenterTitle();
  HSiPM_Multiplicities_G->GetYaxis()->CenterTitle();
  HSiPM_Multiplicities_G->SetDrawOption("COLZ");

  HSiPMLow_False_G = new TH1D(("HSiPMLow_False_G"), ("HSiPMLow_False_G"), eLowN, eLowMin, eLowMax);
  HSiPMLow_False_G->GetXaxis()->SetTitle("SiPM Low [Channel]");
  HSiPMLow_False_G->GetYaxis()->SetTitle("Counts");
  HSiPMLow_False_G->GetXaxis()->CenterTitle();
  HSiPMLow_False_G->GetYaxis()->CenterTitle();

  HSiPMHigh_False_G = new TH1D(("HSiPMHigh_False_G"), ("HSiPMHigh_False_G"), eHighN, eHighMin, eHighMax);
  HSiPMHigh_False_G->GetXaxis()->SetTitle("SiPM High [Channel]");
  HSiPMHigh_False_G->GetYaxis()->SetTitle("Counts");
  HSiPMHigh_False_G->GetXaxis()->CenterTitle();
  HSiPMHigh_False_G->GetYaxis()->CenterTitle();

  HSiPMLow_Time_False_G = new TH1D(("HSiPMLow_Time_False_G"), ("HSiPMLow_Time_False_G"), winLowN, winLowMin, winLowMax);
  HSiPMLow_Time_False_G->GetXaxis()->SetTitle("Time [ns]");
  HSiPMLow_Time_False_G->GetYaxis()->SetTitle("Counts");
  HSiPMLow_Time_False_G->GetXaxis()->CenterTitle();
  HSiPMLow_Time_False_G->GetYaxis()->CenterTitle();

  HSiPMHigh_Time_False_G = new TH1D(("HSiPMHigh_Time_False_G"), ("HSiPMHigh_Time_False_G"), winHighN, winHighMin, winHighMax);
  HSiPMHigh_Time_False_G->GetXaxis()->SetTitle("Time [ns]");
  HSiPMHigh_Time_False_G->GetYaxis()->SetTitle("Counts");
  HSiPMHigh_Time_False_G->GetXaxis()->CenterTitle();
  HSiPMHigh_Time_False_G->GetYaxis()->CenterTitle();

  HSiPMLow_SiPM_False_G = new TH1D(("HSiPMLow_SiPM_False_G"), ("HSiPMLow_SiPM_False_G"), 10, 0, 10);
  HSiPMLow_SiPM_False_G->GetXaxis()->SetTitle("SiPM");
  HSiPMLow_SiPM_False_G->GetYaxis()->SetTitle("Counts");
  HSiPMLow_SiPM_False_G->GetXaxis()->CenterTitle();
  HSiPMLow_SiPM_False_G->GetYaxis()->CenterTitle();

  HSiPMHigh_SiPM_False_G = new TH1D(("HSiPMHigh_SiPM_False_G"), ("HSiPMHigh_SiPM_False_G"), 10, 0, 10);
  HSiPMHigh_SiPM_False_G->GetXaxis()->SetTitle("SiPM");
  HSiPMHigh_SiPM_False_G->GetYaxis()->SetTitle("Counts");
  HSiPMHigh_SiPM_False_G->GetXaxis()->CenterTitle();
  HSiPMHigh_SiPM_False_G->GetYaxis()->CenterTitle();

  HSiPM_Counter_False_G = new TH1D(("HSiPM_Counter_False_G"), ("HSiPM_Counter_False_G"), 5, -2, 2);
  HSiPM_Counter_False_G->GetXaxis()->SetTitle("False Low \t True \t False High");
  HSiPM_Counter_False_G->GetYaxis()->SetTitle("Counts");
  HSiPM_Counter_False_G->GetXaxis()->CenterTitle();
  HSiPM_Counter_False_G->GetYaxis()->CenterTitle();

  for (int i = 1; i <= BETA_SIZE; ++i)
  {
    GSiPM_Channel_C[i] = new TGraph();

    PSiPM_Channel_C[i] = new TProfile(("PSiPM_" + to_string(i)).c_str(), ("PSiPM_" + to_string(i)).c_str(), 500, 0, eLowMax/2, 0, eHighMax);
    PSiPM_Channel_C[i]->GetXaxis()->SetTitle("SiPM Low [Channel]");
    PSiPM_Channel_C[i]->GetYaxis()->SetTitle("SiPM High [Channel]");
    PSiPM_Channel_C[i]->GetXaxis()->CenterTitle();
    PSiPM_Channel_C[i]->GetYaxis()->CenterTitle(); 
  }

  HSiPM12_Channel_C = new TH2D(("HSiPM12_Channel_C"), ("HSiPM12_Channel_C"), eLowN, eLowMin, eLowMax, eHighN, eHighMin, eHighMax);
  HSiPM12_Channel_C->GetXaxis()->SetTitle("SiPM 1 [Channel]");
  HSiPM12_Channel_C->GetYaxis()->SetTitle("SiPM 2 [Channel]");
  HSiPM12_Channel_C->GetXaxis()->CenterTitle();
  HSiPM12_Channel_C->GetYaxis()->CenterTitle();
  HSiPM12_Channel_C->SetDrawOption("COLZ");

  HSiPM13_Channel_C = new TH2D(("HSiPM13_Channel_C"), ("HSiPM13_Channel_C"), eLowN, eLowMin, eLowMax, eHighN, eHighMin, eHighMax);
  HSiPM13_Channel_C->GetXaxis()->SetTitle("SiPM 1 [Channel]");
  HSiPM13_Channel_C->GetYaxis()->SetTitle("SiPM 3 [Channel]");
  HSiPM13_Channel_C->GetXaxis()->CenterTitle();
  HSiPM13_Channel_C->GetYaxis()->CenterTitle();
  HSiPM13_Channel_C->SetDrawOption("COLZ");

  HSiPM23_Channel_C = new TH2D(("HSiPM23_Channel_C"), ("HSiPM23_Channel_C"), eLowN, eLowMin, eLowMax, eHighN, eHighMin, eHighMax);
  HSiPM23_Channel_C->GetXaxis()->SetTitle("SiPM 2 [Channel]");
  HSiPM23_Channel_C->GetYaxis()->SetTitle("SiPM 3 [Channel]");
  HSiPM23_Channel_C->GetXaxis()->CenterTitle();
  HSiPM23_Channel_C->GetYaxis()->CenterTitle();
  HSiPM23_Channel_C->SetDrawOption("COLZ");

  HSiPM12LOW_Channel_C = new TH2D(("HSiPM12LOW_Channel_C"), ("HSiPM12LOW_Channel_C"), eLowN, eLowMin, eLowMax, eLowN, eLowMin, eLowMax);
  HSiPM12LOW_Channel_C->GetXaxis()->SetTitle("SiPM 1 [Channel]");
  HSiPM12LOW_Channel_C->GetYaxis()->SetTitle("SiPM 2 [Channel]");
  HSiPM12LOW_Channel_C->GetXaxis()->CenterTitle();
  HSiPM12LOW_Channel_C->GetYaxis()->CenterTitle();
  HSiPM12LOW_Channel_C->SetDrawOption("COLZ");

  HSiPM13LOW_Channel_C = new TH2D(("HSiPM13LOW_Channel_C"), ("HSiPM13LOW_Channel_C"), eLowN, eLowMin, eLowMax, eLowN, eLowMin, eLowMax);
  HSiPM13LOW_Channel_C->GetXaxis()->SetTitle("SiPM 1 [Channel]");
  HSiPM13LOW_Channel_C->GetYaxis()->SetTitle("SiPM 3 [Channel]");
  HSiPM13LOW_Channel_C->GetXaxis()->CenterTitle();
  HSiPM13LOW_Channel_C->GetYaxis()->CenterTitle();
  HSiPM13LOW_Channel_C->SetDrawOption("COLZ");

  HSiPM23LOW_Channel_C = new TH2D(("HSiPM23LOW_Channel_C"), ("HSiPM23LOW_Channel_C"), eLowN, eLowMin, eLowMax, eLowN, eLowMin, eLowMax);
  HSiPM23LOW_Channel_C->GetXaxis()->SetTitle("SiPM 2 [Channel]");
  HSiPM23LOW_Channel_C->GetYaxis()->SetTitle("SiPM 3 [Channel]");
  HSiPM23LOW_Channel_C->GetXaxis()->CenterTitle();
  HSiPM23LOW_Channel_C->GetYaxis()->CenterTitle();
  HSiPM23LOW_Channel_C->SetDrawOption("COLZ");





  return 0;
}

inline int InitHistograms_Cleaned()
{
  for (size_t i = 0; i < SIGNAL_MAX; ++i)
  {
    HStrip_Channel_C[i] = NULL;
    HRear_Channel_C[i] = NULL;
    HStripRear_Time_C[i] = NULL;
    HStripRear_Channel_C[i] = NULL;
    HFracRearStripChannel_C[i] = NULL;
    HStripsMultiplicity_C[i] = NULL;
    HSiPMHigh_Channel_C[i] = NULL;
    HSiPMHigh_Multiplicity_C = NULL;
    HSiPMHighRear_Time_C[i] = NULL;
    HSiPMHighRear_TimeChannel_C[i] = NULL;
    HSiPMLow_Channel_C[i] = NULL;
    HSiPMLow_Multiplicity_C = NULL;
    HSiPMLowRear_Time_C[i] = NULL;
    HSiPMLowRear_TimeChannel_C[i] = NULL;
    HSiPM_Multiplicities_C = NULL;

    GSiPM_Channel_C[i] = NULL;
    PSiPM_Channel_C[i] = NULL;
  }

  for (size_t i = 0; i < detectorNum; ++i)
  {
    if (IsDetectorSiliStrip(i))
    {
      HStrip_Channel_C[i] = new TH1D(("HStrip_Channel_C_" + detectorName[i]).c_str(), ("HStrip_Channel_C_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      HStrip_Channel_C[i]->GetXaxis()->SetTitle("Strips [Channel]");
      HStrip_Channel_C[i]->GetYaxis()->SetTitle("Counts");
      HStrip_Channel_C[i]->GetXaxis()->CenterTitle();
      HStrip_Channel_C[i]->GetYaxis()->CenterTitle();

      HStripRear_Time_C[i] = new TH1D(("HStripRear_Time_C_" + detectorName[i]).c_str(), ("HStripRear_Time_C_" + detectorName[i]).c_str(), winSiliN, winSiliMin, winSiliMax);
      HStripRear_Time_C[i]->GetXaxis()->SetTitle("Time [ns]");
      HStripRear_Time_C[i]->GetYaxis()->SetTitle("Counts");
      HStripRear_Time_C[i]->GetXaxis()->CenterTitle();
      HStripRear_Time_C[i]->GetYaxis()->CenterTitle();

      HStripRear_Channel_C[i] = new TH2D(("HStripRear_Channel_C_" + detectorName[i]).c_str(), ("HStripRear_Channel_C_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      HStripRear_Channel_C[i]->GetXaxis()->SetTitle("Rear [Channel]");
      HStripRear_Channel_C[i]->GetYaxis()->SetTitle("Strip [Channel]");
      HStripRear_Channel_C[i]->GetXaxis()->CenterTitle();
      HStripRear_Channel_C[i]->GetYaxis()->CenterTitle();
      HStripRear_Channel_C[i]->SetDrawOption("COLZ");

      HFracRearStripChannel_C[i] = new TH1D(("HFracRearStripChannel_C_" + detectorName[i]).c_str(), ("HFracRearStripChannel_C_" + detectorName[i]).c_str(), 2000, 0, 2);
      HFracRearStripChannel_C[i]->GetXaxis()->SetTitle("Rear/Strip");
      HFracRearStripChannel_C[i]->GetYaxis()->SetTitle("Counts");
      HFracRearStripChannel_C[i]->GetXaxis()->CenterTitle();
      HFracRearStripChannel_C[i]->GetYaxis()->CenterTitle();
    }

    if (IsDetectorSiliBack(i))
    {
      HStripsMultiplicity_C[i] = new TH2D(("HStripsMultiplicity_C_" + detectorName[i]).c_str(), ("HStripsMultiplicity_C_" + detectorName[i]).c_str(), 6, 0, 6, 6, 0, 6);
      HStripsMultiplicity_C[i]->GetXaxis()->SetTitle("Strip A");
      HStripsMultiplicity_C[i]->GetYaxis()->SetTitle("Strip B");
      HStripsMultiplicity_C[i]->GetXaxis()->CenterTitle();
      HStripsMultiplicity_C[i]->GetYaxis()->CenterTitle();
      HStripsMultiplicity_C[i]->SetDrawOption("COLZ");

      HRear_Channel_C[i] = new TH1D(("HRear_Channel_C_" + detectorName[i]).c_str(), ("HRear_Channel_C_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      HRear_Channel_C[i]->GetXaxis()->SetTitle("Rear [Channel]");
      HRear_Channel_C[i]->GetYaxis()->SetTitle("Counts");
      HRear_Channel_C[i]->GetXaxis()->CenterTitle();
      HRear_Channel_C[i]->GetYaxis()->CenterTitle();
    }

    if (IsDetectorBetaHigh(i))
    {
      HSiPMHigh_Channel_C[i] = new TH1D(("HSiPMHigh_Channel_C_" + detectorName[i]).c_str(), ("HSiPMHigh_Channel_C_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
      HSiPMHigh_Channel_C[i]->GetXaxis()->SetTitle("SiPM High [Channel]");
      HSiPMHigh_Channel_C[i]->GetYaxis()->SetTitle("Counts");
      HSiPMHigh_Channel_C[i]->GetXaxis()->CenterTitle();
      HSiPMHigh_Channel_C[i]->GetYaxis()->CenterTitle();

      HSiPMHighRear_Time_C[i] = new TH1D(("HSiPMHighRear_Time_C_" + detectorName[i]).c_str(), ("HSiPMHighRear_Time_C_" + detectorName[i]).c_str(), winHighN, winHighMin, winHighMax);
      HSiPMHighRear_Time_C[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMHighRear_Time_C[i]->GetYaxis()->SetTitle("Counts");
      HSiPMHighRear_Time_C[i]->GetXaxis()->CenterTitle();
      HSiPMHighRear_Time_C[i]->GetYaxis()->CenterTitle();

      HSiPMHighRear_TimeChannel_C[i] = new TH2D(("HSiPMHighRear_TimeChannel_C_" + detectorName[i]).c_str(), ("HSiPMHighRear_TimeChannel_C_" + detectorName[i]).c_str(), winHighN, winHighMin, winHighMax, eHighN, eHighMin, eHighMax);
      HSiPMHighRear_TimeChannel_C[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMHighRear_TimeChannel_C[i]->GetYaxis()->SetTitle("SiPM High [Channel]");
      HSiPMHighRear_TimeChannel_C[i]->GetXaxis()->CenterTitle();
      HSiPMHighRear_TimeChannel_C[i]->GetYaxis()->CenterTitle();
      HSiPMHighRear_TimeChannel_C[i]->SetDrawOption("COLZ");
    }

    if (IsDetectorBetaLow(i))
    {
      HSiPMLow_Channel_C[i] = new TH1D(("HSiPMLow_Channel_C_" + detectorName[i]).c_str(), ("HSiPMLow_Channel_C_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax);
      HSiPMLow_Channel_C[i]->GetXaxis()->SetTitle("SiPM Low [Channel]");
      HSiPMLow_Channel_C[i]->GetYaxis()->SetTitle("Counts");
      HSiPMLow_Channel_C[i]->GetXaxis()->CenterTitle();
      HSiPMLow_Channel_C[i]->GetYaxis()->CenterTitle();

      HSiPMLowRear_Time_C[i] = new TH1D(("HSiPMLowRear_Time_C_" + detectorName[i]).c_str(), ("HSiPMLowRear_Time_C_" + detectorName[i]).c_str(), winLowN, winLowMin, winLowMax);
      HSiPMLowRear_Time_C[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMLowRear_Time_C[i]->GetYaxis()->SetTitle("Counts");
      HSiPMLowRear_Time_C[i]->GetXaxis()->CenterTitle();
      HSiPMLowRear_Time_C[i]->GetYaxis()->CenterTitle();

      HSiPMLowRear_TimeChannel_C[i] = new TH2D(("HSiPMLowRear_TimeChannel_C_" + detectorName[i]).c_str(), ("HSiPMLowRear_TimeChannel_C_" + detectorName[i]).c_str(), winLowN, winLowMin, winLowMax, eLowN, eLowMin, eLowMax);
      HSiPMLowRear_TimeChannel_C[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMLowRear_TimeChannel_C[i]->GetYaxis()->SetTitle("SiPM Low [Channel]");
      HSiPMLowRear_TimeChannel_C[i]->GetXaxis()->CenterTitle();
      HSiPMLowRear_TimeChannel_C[i]->GetYaxis()->CenterTitle();
      HSiPMLowRear_TimeChannel_C[i]->SetDrawOption("COLZ");
    }
  }

  for (int i = 1; i <= BETA_SIZE; ++i)
  {
    GSiPM_Channel_C[i] = new TGraph();

    PSiPM_Channel_C[i] = new TProfile(("PSiPM_" + to_string(i)).c_str(), ("PSiPM_" + to_string(i)).c_str(), 500, 0, eLowMax/2, 0, eHighMax);
    PSiPM_Channel_C[i]->GetXaxis()->SetTitle("SiPM Low [Channel]");
    PSiPM_Channel_C[i]->GetYaxis()->SetTitle("SiPM High [Channel]");
    PSiPM_Channel_C[i]->GetXaxis()->CenterTitle();
    PSiPM_Channel_C[i]->GetYaxis()->CenterTitle(); 
  }

  HSiPMHigh_Multiplicity_C = new TH1D(("HSiPMHigh_Multiplicity_C"), ("HSiPMHigh_Multiplicity_C"), 10, 0, 10);
  HSiPMHigh_Multiplicity_C->GetXaxis()->SetTitle("Multiplicity");
  HSiPMHigh_Multiplicity_C->GetYaxis()->SetTitle("Counts");
  HSiPMHigh_Multiplicity_C->GetXaxis()->CenterTitle();
  HSiPMHigh_Multiplicity_C->GetYaxis()->CenterTitle();

  HSiPMLow_Multiplicity_C = new TH1D(("HSiPMLow_Multiplicity_C"), ("HSiPMLow_Multiplicity_C"), 10, 0, 10);
  HSiPMLow_Multiplicity_C->GetXaxis()->SetTitle("Multiplicity");
  HSiPMLow_Multiplicity_C->GetYaxis()->SetTitle("Counts");
  HSiPMLow_Multiplicity_C->GetXaxis()->CenterTitle();
  HSiPMLow_Multiplicity_C->GetYaxis()->CenterTitle();

  HSiPM_Multiplicities_C = new TH2D(("HSiPM_Multiplicities_C"), ("HSiPM_Multiplicities_C"), 10, 0, 10, 10, 0, 10);
  HSiPM_Multiplicities_C->GetXaxis()->SetTitle("Multiplicity High");
  HSiPM_Multiplicities_C->GetYaxis()->SetTitle("Multiplicity Low");
  HSiPM_Multiplicities_C->GetXaxis()->CenterTitle();
  HSiPM_Multiplicities_C->GetYaxis()->CenterTitle();
  HSiPM_Multiplicities_C->SetDrawOption("COLZ");
  return 0;
}

inline int InitTree_Grouped()
{
  Tree_Grouped = new TTree("Tree", "Tree");
  Tree_Grouped->Branch("Tree_Event", &Tree_Event, "Tree_Event/I");
  Tree_Grouped->Branch("Tree_Silicon", &Tree_Silicon);
  Tree_Grouped->Branch("Tree_SiPMHigh", &Tree_SiPMHigh);
  Tree_Grouped->Branch("Tree_SiPMLow", &Tree_SiPMLow);
  return 0;
}

inline int InitTree_Cleaned()
{
  Tree_Cleaned = new TTree("Tree", "Tree");
  Tree_Cleaned->Branch("Tree_Event", &Tree_Cleaned_Event, "Tree_Cleaned_Event/I");
  Tree_Cleaned->Branch("Tree_Silicon", &Tree_Cleaned_Silicon);
  Tree_Cleaned->Branch("Tree_SiPMHigh", &Tree_Cleaned_SiPMHigh);
  Tree_Cleaned->Branch("Tree_SiPMLow", &Tree_Cleaned_SiPMLow);

  for (size_t i = 0; i < SIGNAL_MAX; ++i)
  {
    if (IsDetectorSili(i))
    {
      Tree_Silicons[i] = new TTree(("Tree_Silicon_" + detectorName[i]).c_str(), ("Tree_Silicon_" + detectorName[i]).c_str());
      Tree_Silicons[i]->Branch("Channel", &Tree_Channel, "Channel/D");
    }
  }

  return 0;
}

inline int WriteHistograms_Grouped()
{
  File_Grouped->cd();
  TDirectory *SiPMHigh_Channel_G = File_Grouped->mkdir("SiPMHigh_Channel");
  TDirectory *SiPMHigh_Time_G = File_Grouped->mkdir("SiPMHigh_Time");
  TDirectory *SiPMHighRear_TimeChannel_G = File_Grouped->mkdir("SiPMHighRear_TimeChannel");
  TDirectory *SiPMHigh_Multiplicity_G = File_Grouped->mkdir("SiPMHigh_Multiplicity");

  TDirectory *SiPMLow_Channel_G = File_Grouped->mkdir("SiPMLow_Channel");
  TDirectory *SiPMLow_Time_G = File_Grouped->mkdir("SiPMLow_Time");
  TDirectory *SiPMLowRear_TimeChannel_G = File_Grouped->mkdir("SiPMLowRear_TimeChannel");
  TDirectory *SiPMLow_Multiplicity_G = File_Grouped->mkdir("SiPMLow_Multiplicity");

  TDirectory *SiPM_Multiplicities_G = File_Grouped->mkdir("SiPM_Multiplicities");

  for (size_t i = 0; i < detectorNum; ++i)
  {
    if (IsDetectorBetaHigh(i))
    {
      SiPMHigh_Channel_G->cd();
      HSiPMHigh_Channel_G[i]->Write();
      delete HSiPMHigh_Channel_G[i];
      SiPMHigh_Time_G->cd();
      HSiPMHighRear_Time_G[i]->Write();
      delete HSiPMHighRear_Time_G[i];
      SiPMHighRear_TimeChannel_G->cd();
      HSiPMHighRear_TimeChannel_G[i]->Write();
      delete HSiPMHighRear_TimeChannel_G[i];
    }

    if (IsDetectorBetaLow(i))
    {
      SiPMLow_Channel_G->cd();
      HSiPMLow_Channel_G[i]->Write();
      delete HSiPMLow_Channel_G[i];
      SiPMLow_Time_G->cd();
      HSiPMLowRear_Time_G[i]->Write();
      delete HSiPMLowRear_Time_G[i];
      SiPMLowRear_TimeChannel_G->cd();
      HSiPMLowRear_TimeChannel_G[i]->Write();
      delete HSiPMLowRear_TimeChannel_G[i];
    }
  }

  SiPMHigh_Multiplicity_G->cd();
  HSiPMHigh_Multiplicity_G->Write();
  delete HSiPMHigh_Multiplicity_G;
  SiPMLow_Multiplicity_G->cd();
  HSiPMLow_Multiplicity_G->Write();
  delete HSiPMLow_Multiplicity_G;
  SiPM_Multiplicities_G->cd();
  HSiPM_Multiplicities_G->Write();
  delete HSiPM_Multiplicities_G;

  SiPM_Multiplicities_G->cd();
  HSiPMHigh_False_G->Write();
  delete HSiPMHigh_False_G;
  HSiPMLow_False_G->Write();
  delete HSiPMLow_False_G;
  HSiPMHigh_Time_False_G->Write();
  delete HSiPMHigh_Time_False_G;
  HSiPMLow_Time_False_G->Write();
  delete HSiPMLow_Time_False_G;
  HSiPMHigh_SiPM_False_G->Write();
  delete HSiPMHigh_SiPM_False_G;
  HSiPMLow_SiPM_False_G->Write();
  delete HSiPMLow_SiPM_False_G;

  HSiPM_Counter_False_G->Fill(-1., counter_low);
  HSiPM_Counter_False_G->Fill(0., counter_true);
  HSiPM_Counter_False_G->Fill(1., counter_high);
  HSiPM_Counter_False_G->Write();
  delete HSiPM_Counter_False_G;

  SiPM_Multiplicities_G->cd();
  delete HSiPM_Multiplicities_C;

  for (int i = 1; i <= BETA_SIZE; ++i)
  {
    TCanvas *canvas = new TCanvas(("SiPM_" + to_string(i)).c_str(), ("SiPM_" + to_string(i)).c_str(), 200, 10, 700, 500);
    canvas->cd();
    GSiPM_Channel_C[i]->SetTitle(("SiPM_" + to_string(i)).c_str());
    GSiPM_Channel_C[i]->SetName(("SiPM_" + to_string(i)).c_str());
    GSiPM_Channel_C[i]->Draw("AP");
    canvas->Write();
    delete canvas;

    PSiPM_Channel_C[i]->Write();
    delete PSiPM_Channel_C[i];
  }

  HSiPM12_Channel_C->Write();
  HSiPM13_Channel_C->Write();
  HSiPM23_Channel_C->Write();

  HSiPM12LOW_Channel_C->Write();
  HSiPM13LOW_Channel_C->Write();
  HSiPM23LOW_Channel_C->Write();

  return 0;
}

inline int WriteHistograms_Cleaned()
{
  File_Cleaned->cd();
  TDirectory *Strip_Channel_C = File_Cleaned->mkdir("Strip_Channel");
  TDirectory *Rear_Channel_C = File_Cleaned->mkdir("Rear_Channel");
  TDirectory *Strip_Time_C = File_Cleaned->mkdir("Strip_Time");
  TDirectory *Rear_Strip_Channel_C = File_Cleaned->mkdir("Rear_Strip_Channel");
  TDirectory *fracRear_Strip_Channel_C = File_Cleaned->mkdir("Frac_Rear_Strip_Channel");
  TDirectory *Strip_Multiplicities_C = File_Cleaned->mkdir("Strip_Multiplicities");

  TDirectory *SiPMHigh_Channel_C = File_Cleaned->mkdir("SiPMHigh_Channel");
  TDirectory *SiPMHigh_Time_C = File_Cleaned->mkdir("SiPMHigh_Time");
  TDirectory *SiPMHighRear_TimeChannel_C = File_Cleaned->mkdir("SiPMHighRear_TimeChannel");
  TDirectory *SiPMHigh_Multiplicity_C = File_Cleaned->mkdir("SiPMHigh_Multiplicity");

  TDirectory *SiPMLow_Channel_C = File_Cleaned->mkdir("SiPMLow_Channel");
  TDirectory *SiPMLow_Time_C = File_Cleaned->mkdir("SiPMLow_Time");
  TDirectory *SiPMLowRear_TimeChannel_C = File_Cleaned->mkdir("SiPMLowRear_TimeChannel");
  TDirectory *SiPMLow_Multiplicity_C = File_Cleaned->mkdir("SiPMLow_Multiplicity");

  TDirectory *SiPM_Multiplicities_C = File_Cleaned->mkdir("SiPM_Multiplicities");

  for (size_t i = 0; i < detectorNum; ++i)
  {
    if (IsDetectorSiliStrip(i))
    {
      Strip_Channel_C->cd();
      HStrip_Channel_C[i]->Write();
      delete HStrip_Channel_C[i];
      Strip_Time_C->cd();
      HStripRear_Time_C[i]->Write();
      delete HStripRear_Time_C[i];
      Rear_Strip_Channel_C->cd();
      HStripRear_Channel_C[i]->Write();
      delete HStripRear_Channel_C[i];
      fracRear_Strip_Channel_C->cd();
      HFracRearStripChannel_C[i]->Write();
      delete HFracRearStripChannel_C[i];
    }

    if (IsDetectorSiliBack(i))
    {
      Strip_Multiplicities_C->cd();
      HStripsMultiplicity_C[i]->Write();
      delete HStripsMultiplicity_C[i];
      Rear_Channel_C->cd();
      HRear_Channel_C[i]->Write();
      delete HRear_Channel_C[i];
    }

    if (IsDetectorBetaHigh(i))
    {
      SiPMHigh_Channel_C->cd();
      HSiPMHigh_Channel_C[i]->Write();
      delete HSiPMHigh_Channel_C[i];
      SiPMHigh_Time_C->cd();
      HSiPMHighRear_Time_C[i]->Write();
      delete HSiPMHighRear_Time_C[i];
      SiPMHighRear_TimeChannel_C->cd();
      HSiPMHighRear_TimeChannel_C[i]->Write();
      delete HSiPMHighRear_TimeChannel_C[i];
    }

    if (IsDetectorBetaLow(i))
    {
      SiPMLow_Channel_C->cd();
      HSiPMLow_Channel_C[i]->Write();
      delete HSiPMLow_Channel_C[i];
      SiPMLow_Time_C->cd();
      HSiPMLowRear_Time_C[i]->Write();
      delete HSiPMLowRear_Time_C[i];
      SiPMLowRear_TimeChannel_C->cd();
      HSiPMLowRear_TimeChannel_C[i]->Write();
      delete HSiPMLowRear_TimeChannel_C[i];
    }
  }

  SiPMHigh_Multiplicity_C->cd();
  HSiPMHigh_Multiplicity_C->Write();
  delete HSiPMHigh_Multiplicity_C;

  SiPMLow_Multiplicity_C->cd();
  HSiPMLow_Multiplicity_C->Write();
  delete HSiPMLow_Multiplicity_C;

  ////TCANVAS

  SiPM_Multiplicities_C->cd();
  HSiPM_Multiplicities_C->Write();
  delete HSiPM_Multiplicities_C;

  for (int i = 1; i <= BETA_SIZE; ++i)
  {
    TCanvas *canvas = new TCanvas(("SiPM_" + to_string(i)).c_str(), ("SiPM_" + to_string(i)).c_str(), 200, 10, 700, 500);
    canvas->cd();
    GSiPM_Channel_C[i]->SetTitle(("SiPM_" + to_string(i)).c_str());
    GSiPM_Channel_C[i]->SetName(("SiPM_" + to_string(i)).c_str());
    GSiPM_Channel_C[i]->Draw("AP");
    canvas->Write();
    delete canvas;

    PSiPM_Channel_C[i]->Write();
    delete PSiPM_Channel_C[i];
  }

  return 0;
}

inline int WriteTree_Grouped()
{
  File_Grouped->cd();
  Tree_Grouped->Write();
  delete Tree_Grouped;
  return 0;
}

inline int WriteTree_Cleaned()
{
  File_Cleaned->cd();
  Tree_Cleaned->Write();

  for (size_t i = 0; i < SIGNAL_MAX; ++i)
  {
    if (IsDetectorSili(i))
    {
      Tree_Silicons[i]->Write();
      delete Tree_Silicons[i];
    }
  }

  delete Tree_Cleaned;
  return 0;
}

void ProgressBar(ULong64_t cEntry, ULong64_t TotalEntries, clock_t start, clock_t Current)
{
  Current = clock();
  const Char_t *Color;
  Double_t Frac = 1.0 * cEntry / TotalEntries;
  Double_t Timeclock = ((double)(Current - start) / CLOCKS_PER_SEC);
  Double_t TimeLeft = Timeclock * (1 / Frac - 1.);
  Color = "\e[1;31m";
  cout << Form("\r%sEntry : %9llu", Color, cEntry)
       << "/" << TotalEntries
       << " --- "
       << Form("%4.2f", 100. * cEntry / TotalEntries) << " %"
       << " --- "
       << Form("%7.00f RunEvt/sec", cEntry / Timeclock)
       << " --- "
       << " Time Left : " << Form("%2d min ", (int)TimeLeft / 60)
       << Form("%02d sec", (int)TimeLeft % 60)
       << flush;
}

void SavingData_Grouped(Signal TrigSignal, vector<Signal> SiPM_High, vector<Signal> SiPM_Low)
{
  /// Saving in Hist
  if (TrigSignal.isValid)
  {
    for (auto it = SiPM_High.begin(); it != SiPM_High.end(); ++it)
    {
      HSiPMHigh_Channel_G[it->Label]->Fill(it->Channel);
      HSiPMHighRear_Time_G[it->Label]->Fill(it->Time - TrigSignal.Time);
      HSiPMHighRear_TimeChannel_G[it->Label]->Fill(it->Time - TrigSignal.Time, it->Channel);
    }

    for (auto it = SiPM_Low.begin(); it != SiPM_Low.end(); ++it)
    {
      HSiPMLow_Channel_G[it->Label]->Fill(it->Channel);
      HSiPMLowRear_Time_G[it->Label]->Fill(it->Time - TrigSignal.Time);
      HSiPMLowRear_TimeChannel_G[it->Label]->Fill(it->Time - TrigSignal.Time, it->Channel);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    SiPM_High.erase(std::remove_if(SiPM_High.begin(), SiPM_High.end(), [](const Signal &s)
                                   { return GetDetectorChannel(s.Label) == 4; }),
                    SiPM_High.end());
    SiPM_Low.erase(std::remove_if(SiPM_Low.begin(), SiPM_Low.end(), [](const Signal &s)
                                  { return GetDetectorChannel(s.Label) == 9; }),
                   SiPM_Low.end());


    for (auto ith = SiPM_High.begin(); ith != SiPM_High.end(); ++ith)
    {
      for (auto itl = SiPM_Low.begin(); itl != SiPM_Low.end(); ++itl)
      {
        if (GetDetectorChannel(ith->Label) == GetDetectorChannel(itl->Label))
        {
          GSiPM_Channel_C[GetDetectorChannel(ith->Label)]->SetPoint(counter_graph[GetDetectorChannel(ith->Label)], itl->Channel, ith->Channel);
          PSiPM_Channel_C[GetDetectorChannel(ith->Label)]->Fill(itl->Channel, ith->Channel);
          counter_graph[GetDetectorChannel(ith->Label)]++;
        }
      }
    }

    Signal one;
    Signal two;
    Signal three;
    for (auto ith = SiPM_High.begin(); ith != SiPM_High.end(); ++ith)
    {
      if (GetDetectorChannel(ith->Label) == 1)
      {
        one = *ith;
      }

      if (GetDetectorChannel(ith->Label) == 2)
      {
        two = *ith;
      }

      if (GetDetectorChannel(ith->Label) == 3)
      {
        three = *ith;
      }
    }

    if (one.isValid && two.isValid)
    {
      HSiPM12_Channel_C->Fill(one.Channel, two.Channel);
    }
    if (one.isValid && three.isValid)
    {
      HSiPM13_Channel_C->Fill(one.Channel, three.Channel);
    }
    if (two.isValid && three.isValid)
    {
      HSiPM23_Channel_C->Fill(two.Channel, three.Channel);
    }

    Signal onel;
    Signal twol;
    Signal threel;
    for (auto ith = SiPM_Low.begin(); ith != SiPM_Low.end(); ++ith)
    {
      if (GetDetectorChannel(ith->Label) == 1)
      {
        onel = *ith;
      }

      if (GetDetectorChannel(ith->Label) == 2)
      {
        twol = *ith;
      }

      if (GetDetectorChannel(ith->Label) == 3)
      {
        threel = *ith;
      }
    }

    if (onel.isValid && twol.isValid)
    {
      HSiPM12LOW_Channel_C->Fill(onel.Channel, twol.Channel);
    }
    if (onel.isValid && threel.isValid)
    {
      HSiPM13LOW_Channel_C->Fill(onel.Channel, threel.Channel);
    }
    if (twol.isValid && threel.isValid)
    {
      HSiPM23LOW_Channel_C->Fill(twol.Channel, threel.Channel);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////

    HSiPMHigh_Multiplicity_G->Fill(SiPM_High.size());
    HSiPMLow_Multiplicity_G->Fill(SiPM_Low.size());
    HSiPM_Multiplicities_G->Fill(SiPM_High.size(), SiPM_Low.size());
  }

  else
  {
    //////TO DO
    // Analyse no-coincidence
  }

  /// Saving in Tree
  

    Event++;
    Tree_Event = Event;
    Tree_SiPMHigh = SiPM_High;
    Tree_SiPMLow = SiPM_Low;
    Tree_Grouped->Fill();
  
}

void SavingData_Cleaned(vector<Signal> SiliconRear, vector<Signal> SiPM_High, vector<Signal> SiPM_Low)
{
  Signal TrigSignal = SiliconRear[1];
  Signal Silicon = SiliconRear[0];
  /// Saving in Hist
  double Total = 0;

  if (Silicon.isValid)
  {
    HStrip_Channel_C[Silicon.Label]->Fill(Silicon.Channel);
    HRear_Channel_C[TrigSignal.Label]->Fill(TrigSignal.Channel);
    HStripRear_Time_C[Silicon.Label]->Fill(Silicon.Time - TrigSignal.Time);
    HStripRear_Channel_C[Silicon.Label]->Fill(TrigSignal.Channel, Silicon.Channel);
    HFracRearStripChannel_C[Silicon.Label]->Fill(static_cast<double>(Silicon.Channel) / TrigSignal.Channel);

    for (auto it = SiPM_High.begin(); it != SiPM_High.end(); ++it)
    {
      HSiPMHigh_Channel_C[it->Label]->Fill(it->Channel);
      HSiPMHighRear_Time_C[it->Label]->Fill(it->Time - TrigSignal.Time);
      HSiPMHighRear_TimeChannel_C[it->Label]->Fill(it->Time - TrigSignal.Time, it->Channel);
    }

    for (auto it = SiPM_Low.begin(); it != SiPM_Low.end(); ++it)
    {
      HSiPMLow_Channel_C[it->Label]->Fill(it->Channel);
      HSiPMLowRear_Time_C[it->Label]->Fill(it->Time - TrigSignal.Time);
      HSiPMLowRear_TimeChannel_C[it->Label]->Fill(it->Time - TrigSignal.Time, it->Channel);
    }

    HSiPMHigh_Multiplicity_C->Fill(SiPM_High.size());
    HSiPMLow_Multiplicity_C->Fill(SiPM_Low.size());
    HSiPM_Multiplicities_C->Fill(SiPM_High.size(), SiPM_Low.size());

    for (auto ith = SiPM_High.begin(); ith != SiPM_High.end(); ++ith)
    {
      for (auto itl = SiPM_Low.begin(); itl != SiPM_Low.end(); ++itl)
      {
        if (GetDetectorChannel(ith->Label) == GetDetectorChannel(itl->Label))
        {
          GSiPM_Channel_C[GetDetectorChannel(ith->Label)]->SetPoint(counter_graph[GetDetectorChannel(ith->Label)], itl->Channel, ith->Channel);
          PSiPM_Channel_C[GetDetectorChannel(ith->Label)]->Fill(itl->Channel, ith->Channel);
          counter_graph[GetDetectorChannel(ith->Label)]++;
        }
      }
    }
  }

  

  else
  {
    //////TO DO
    // Analyse no-coincidence
  }
}

// Signal ProcessSilicon(vector<Signal> Signal_Vector, Signal Trigger)
// {
//   vector<Signal> Output_Vector = vector<Signal>();
//   /// Search in Signal Vector between 1 and n
//   for (auto current = Signal_Vector.begin(); current != Signal_Vector.end(); ++current)
//   {
//     if (current->Time - Trigger.Time < winSiliMax && Trigger.Time - current->Time > winSiliMin && IsSameSiliDetector(Trigger.Label, current->Label) && Trigger != *current) //&& current->Channel < 42000 && current->Channel > 40000)
//     {
//       Output_Vector.push_back(*current);
//     }
//   }

//   if (Output_Vector.size() != 1)
//   {
//     if (Output_Vector.size() == 2)
//     {
//       if (Output_Vector[0].Label < Output_Vector[1].Label)
//         HStripsMultiplicity_G[Trigger.Label]->Fill(GetDetectorChannel(Output_Vector[0].Label), GetDetectorChannel(Output_Vector[1].Label));
//       else
//         HStripsMultiplicity_G[Trigger.Label]->Fill(GetDetectorChannel(Output_Vector[1].Label), GetDetectorChannel(Output_Vector[0].Label));
//     }

//     return Signal();
//   }

//   return Output_Vector[0];
// }

vector<Signal> ProcessSiPM(vector<Signal> Signal_Vector, Signal Trigger)
{
  vector<Signal> Output_Vector = vector<Signal>();

  /// Search in Signal Vector between 1 and n
  for (auto current = Signal_Vector.begin(); current != Signal_Vector.end(); ++current)
  {
    if (current->Time - Trigger.Time < winHighMax && current->Time - Trigger.Time > winHighMin)
    {
      Output_Vector.push_back(*current);
    }
  }

  return Output_Vector;
}

void SearchForCoincidence(TTreeReader *Reader, TTreeReaderValue<double> &Times, TTreeReaderValue<int> &Labels, TTreeReaderValue<double> &Channels)
{
  if (Verbose > 1)
    cout << "### Starting a Coincidence group on " << detectorName[*Labels] << " DATA nb : " << Reader->GetCurrentEntry() << "  TIME : " << *Times << endl;

  // Vector save
  vector<Signal> SiPMSignalsHigh;
  vector<Signal> SiPMSignalsLow;

  // Triger variables
  Signal TrigSignal = Signal(*Labels, *Times, *Channels);
  int TrigEntry = Reader->GetCurrentEntry();

  // Searching the first event of the group
  int counter = 0;
  while (*Times >= TrigSignal.Time + winTotalMin && TrigEntry - counter > 0)
  {
    counter++;
    Reader->SetEntry(TrigEntry - counter);
  }

  // Starting at the first event of the group
  int GroupEntry = TrigEntry - counter + 1;
  Reader->SetEntry(GroupEntry);

  while (*Times <= TrigSignal.Time + winTotalMax && Reader->GetCurrentEntry() < Reader->GetEntries())
  {
    Signal current = Signal(*Labels, *Times, *Channels); /// Current Signal
    if (IsDetectorBetaLow(*Labels))
    {
      SiPMSignalsLow.push_back(current);
    }
    else if (IsDetectorBetaHigh(*Labels))
    {
      SiPMSignalsHigh.push_back(current);
    }

    Reader->Next();
  }

  int EndofGroup = Reader->GetCurrentEntry();
  Reader->SetEntry(EndofGroup - 1);

  ////Raw Save

  /// Selecting Signals in time window
  vector<Signal> SiPM_High = ProcessSiPM(SiPMSignalsHigh, TrigSignal); /// Time Window
  vector<Signal> SiPM_Low = ProcessSiPM(SiPMSignalsLow, TrigSignal);   /// Time Window

  if (IsDetectorBetaHigh(TrigSignal.Label))
  {
    SiPM_High.push_back(TrigSignal);
  }
  else if (IsDetectorBetaLow(TrigSignal.Label))
  {
    SiPM_Low.push_back(TrigSignal);
  }

  SavingData_Grouped(TrigSignal, SiPM_High, SiPM_Low);
}

void InitCleaning()
{
  for (size_t i = 0; i < detectorNum; ++i)
  {
    if (IsDetectorSiliStrip(i))
    {
      TH1D *H = (TH1D *)File_Grouped->Get(("Frac_Rear_Strip_Channel/HFracRearStripChannel_G_" + detectorName[i]).c_str());
      H->GetXaxis()->SetRangeUser(1.0, 1.5);

      detectorCleaning[i] = make_pair(H->GetBinCenter(H->GetMaximumBin()), H->GetRMS());
    }
  }
}

#endif