#ifndef MERGER_GROUPER_HH
#define MERGER_GROUPER_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream> // Pour l'écriture dans un fichier

#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TTree.h"
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
int Tree_Channel;

//////////////GROUPED////////////////
/// Silicon
TH1I *HStrip_Channel_G[SIGNAL_MAX];
TH1I *HStripRear_Time_G[SIGNAL_MAX];
TH2I *HStripRear_Channel_G[SIGNAL_MAX];
TH1D *HFracRearStripChannel_G[SIGNAL_MAX];

TH2I *HStripsMultiplicity_G[SIGNAL_MAX];

/// SiPM High
TH1I *HSiPMHigh_Channel_G[SIGNAL_MAX];
TH1I *HSiPMHighRear_Time_G[SIGNAL_MAX];
TH2I *HSiPMHighRear_TimeChannel_G[SIGNAL_MAX];

/// SiPM Low
TH1I *HSiPMLow_Channel_G[SIGNAL_MAX];
TH1I *HSiPMLowRear_Time_G[SIGNAL_MAX];
TH2I *HSiPMLowRear_TimeChannel_G[SIGNAL_MAX];

/// SIPMs
TH1I *HSiPMHigh_Multiplicity_G;
TH1I *HSiPMLow_Multiplicity_G;
TH2I *HSiPM_Multiplicities_G;
/////////////////////////////////////
//////////////CLEANED////////////////
/// Silicon
TH1I *HStrip_Channel_C[SIGNAL_MAX];
TH1I *HRear_Channel_C[SIGNAL_MAX];
TH1I *HStripRear_Time_C[SIGNAL_MAX];
TH2I *HStripRear_Channel_C[SIGNAL_MAX];
TH1D *HFracRearStripChannel_C[SIGNAL_MAX];

TH2I *HStripsMultiplicity_C[SIGNAL_MAX];

/// SiPM High
TH1I *HSiPMHigh_Channel_C[SIGNAL_MAX];
TH1I *HSiPMHighRear_Time_C[SIGNAL_MAX];
TH2I *HSiPMHighRear_TimeChannel_C[SIGNAL_MAX];

/// SiPM Low
TH1I *HSiPMLow_Channel_C[SIGNAL_MAX];
TH1I *HSiPMLowRear_Time_C[SIGNAL_MAX];
TH2I *HSiPMLowRear_TimeChannel_C[SIGNAL_MAX];

/// SIPMs
TH1I *HSiPMHigh_Multiplicity_C;
TH2I *HSiPMHigh_Correlation_C[BETA_SIZE][BETA_SIZE];
TH2I *HSiPMHigh_SUMCorrelation_C[BETA_SIZE];
TCanvas *CanvasHigh;
TH1I *HSiPMLow_Multiplicity_C;
TH2I *HSiPMLow_Correlation_C[BETA_SIZE][BETA_SIZE];
TH2I *HSiPMLow_SUMCorrelation_C[BETA_SIZE];
TCanvas *CanvasLow;
TH2I *HSiPM_Multiplicities_C;

TTree *Tree_Silicons[SIGNAL_MAX];
TTree *Tree_SiPMs[SIGNAL_MAX];
/////////////////////////////////////

inline int InitHistograms_Grouped()
{
  for (size_t i = 0; i < SIGNAL_MAX; ++i)
  {
    HStrip_Channel_G[i] = NULL;
    HStripRear_Time_G[i] = NULL;
    HStripRear_Channel_G[i] = NULL;
    HFracRearStripChannel_G[i] = NULL;
    HStripsMultiplicity_G[i] = NULL;
    HSiPMHigh_Channel_G[i] = NULL;
    HSiPMHigh_Multiplicity_G = NULL;
    HSiPMHighRear_Time_G[i] = NULL;
    HSiPMHighRear_TimeChannel_G[i] = NULL;
    HSiPMLow_Channel_G[i] = NULL;
    HSiPMLow_Multiplicity_G = NULL;
    HSiPMLowRear_Time_G[i] = NULL;
    HSiPMLowRear_TimeChannel_G[i] = NULL;
    HSiPM_Multiplicities_G = NULL;
  }

  for (size_t i = 0; i < detectorNum; ++i)
  {
    if (IsDetectorSiliStrip(i))
    {
      HStrip_Channel_G[i] = new TH1I(("HStrip_Channel_G_" + detectorName[i]).c_str(), ("HStrip_Channel_G_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      HStrip_Channel_G[i]->GetXaxis()->SetTitle("Strips [Channel]");
      HStrip_Channel_G[i]->GetYaxis()->SetTitle("Counts");
      HStrip_Channel_G[i]->GetXaxis()->CenterTitle();
      HStrip_Channel_G[i]->GetYaxis()->CenterTitle();

      HStripRear_Time_G[i] = new TH1I(("HStripRear_Time_G_" + detectorName[i]).c_str(), ("HStripRear_Time_G_" + detectorName[i]).c_str(), winSiliN, winSiliMin, winSiliMax);
      HStripRear_Time_G[i]->GetXaxis()->SetTitle("Time [ns]");
      HStripRear_Time_G[i]->GetYaxis()->SetTitle("Counts");
      HStripRear_Time_G[i]->GetXaxis()->CenterTitle();
      HStripRear_Time_G[i]->GetYaxis()->CenterTitle();

      HStripRear_Channel_G[i] = new TH2I(("HStripRear_Channel_G_" + detectorName[i]).c_str(), ("HStripRear_Channel_G_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
      HStripRear_Channel_G[i]->GetXaxis()->SetTitle("Rear [Channel]");
      HStripRear_Channel_G[i]->GetYaxis()->SetTitle("Strip [Channel]");
      HStripRear_Channel_G[i]->GetXaxis()->CenterTitle();
      HStripRear_Channel_G[i]->GetYaxis()->CenterTitle();
      HStripRear_Channel_G[i]->SetDrawOption("COLZ");

      HFracRearStripChannel_G[i] = new TH1D(("HFracRearStripChannel_G_" + detectorName[i]).c_str(), ("HFracRearStripChannel_G_" + detectorName[i]).c_str(), 2000, 0, 2);
      HFracRearStripChannel_G[i]->GetXaxis()->SetTitle("Rear/Strip");
      HFracRearStripChannel_G[i]->GetYaxis()->SetTitle("Counts");
      HFracRearStripChannel_G[i]->GetXaxis()->CenterTitle();
      HFracRearStripChannel_G[i]->GetYaxis()->CenterTitle();
    }

    if (IsDetectorSiliBack(i))
    {
      HStripsMultiplicity_G[i] = new TH2I(("HStripsMultiplicity_G_" + detectorName[i]).c_str(), ("HStripsMultiplicity_G_" + detectorName[i]).c_str(), 6, 0, 6, 6, 0, 6);
      HStripsMultiplicity_G[i]->GetXaxis()->SetTitle("Strip A");
      HStripsMultiplicity_G[i]->GetYaxis()->SetTitle("Strip B");
      HStripsMultiplicity_G[i]->GetXaxis()->CenterTitle();
      HStripsMultiplicity_G[i]->GetYaxis()->CenterTitle();
      HStripsMultiplicity_G[i]->SetDrawOption("COLZ");
    }

    if (IsDetectorBetaHigh(i))
    {
      HSiPMHigh_Channel_G[i] = new TH1I(("HSiPMHigh_Channel_G_" + detectorName[i]).c_str(), ("HSiPMHigh_Channel_G_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
      HSiPMHigh_Channel_G[i]->GetXaxis()->SetTitle("SiPM High [Channel]");
      HSiPMHigh_Channel_G[i]->GetYaxis()->SetTitle("Counts");
      HSiPMHigh_Channel_G[i]->GetXaxis()->CenterTitle();
      HSiPMHigh_Channel_G[i]->GetYaxis()->CenterTitle();

      HSiPMHighRear_Time_G[i] = new TH1I(("HSiPMHighRear_Time_G_" + detectorName[i]).c_str(), ("HSiPMHighRear_Time_G_" + detectorName[i]).c_str(), winHighN, winHighMin, winHighMax);
      HSiPMHighRear_Time_G[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMHighRear_Time_G[i]->GetYaxis()->SetTitle("Counts");
      HSiPMHighRear_Time_G[i]->GetXaxis()->CenterTitle();
      HSiPMHighRear_Time_G[i]->GetYaxis()->CenterTitle();

      HSiPMHighRear_TimeChannel_G[i] = new TH2I(("HSiPMHighRear_TimeChannel_G_" + detectorName[i]).c_str(), ("HSiPMHighRear_TimeChannel_G_" + detectorName[i]).c_str(), winHighN, winHighMin, winHighMax, eHighN, eHighMin, eHighMax);
      HSiPMHighRear_TimeChannel_G[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMHighRear_TimeChannel_G[i]->GetYaxis()->SetTitle("SiPM High [Channel]");
      HSiPMHighRear_TimeChannel_G[i]->GetXaxis()->CenterTitle();
      HSiPMHighRear_TimeChannel_G[i]->GetYaxis()->CenterTitle();
      HSiPMHighRear_TimeChannel_G[i]->SetDrawOption("COLZ");
    }

    if (IsDetectorBetaLow(i))
    {
      HSiPMLow_Channel_G[i] = new TH1I(("HSiPMLow_Channel_G_" + detectorName[i]).c_str(), ("HSiPMLow_Channel_G_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax);
      HSiPMLow_Channel_G[i]->GetXaxis()->SetTitle("SiPM Low [Channel]");
      HSiPMLow_Channel_G[i]->GetYaxis()->SetTitle("Counts");
      HSiPMLow_Channel_G[i]->GetXaxis()->CenterTitle();
      HSiPMLow_Channel_G[i]->GetYaxis()->CenterTitle();

      HSiPMLowRear_Time_G[i] = new TH1I(("HSiPMLowRear_Time_G_" + detectorName[i]).c_str(), ("HSiPMLowRear_Time_G_" + detectorName[i]).c_str(), winLowN, winLowMin, winLowMax);
      HSiPMLowRear_Time_G[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMLowRear_Time_G[i]->GetYaxis()->SetTitle("Counts");
      HSiPMLowRear_Time_G[i]->GetXaxis()->CenterTitle();
      HSiPMLowRear_Time_G[i]->GetYaxis()->CenterTitle();

      HSiPMLowRear_TimeChannel_G[i] = new TH2I(("HSiPMLowRear_TimeChannel_G_" + detectorName[i]).c_str(), ("HSiPMLowRear_TimeChannel_G_" + detectorName[i]).c_str(), winLowN, winLowMin, winLowMax, eLowN, eLowMin, eLowMax);
      HSiPMLowRear_TimeChannel_G[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMLowRear_TimeChannel_G[i]->GetYaxis()->SetTitle("SiPM Low [Channel]");
      HSiPMLowRear_TimeChannel_G[i]->GetXaxis()->CenterTitle();
      HSiPMLowRear_TimeChannel_G[i]->GetYaxis()->CenterTitle();
      HSiPMLowRear_TimeChannel_G[i]->SetDrawOption("COLZ");
    }
  }
  HSiPMHigh_Multiplicity_G = new TH1I(("HSiPMHigh_Multiplicity_G"), ("HSiPMHigh_Multiplicity_G"), 10, 0, 10);
  HSiPMHigh_Multiplicity_G->GetXaxis()->SetTitle("Multiplicity");
  HSiPMHigh_Multiplicity_G->GetYaxis()->SetTitle("Counts");
  HSiPMHigh_Multiplicity_G->GetXaxis()->CenterTitle();
  HSiPMHigh_Multiplicity_G->GetYaxis()->CenterTitle();

  HSiPMLow_Multiplicity_G = new TH1I(("HSiPMLow_Multiplicity_G"), ("HSiPMLow_Multiplicity_G"), 10, 0, 10);
  HSiPMLow_Multiplicity_G->GetXaxis()->SetTitle("Multiplicity");
  HSiPMLow_Multiplicity_G->GetYaxis()->SetTitle("Counts");
  HSiPMLow_Multiplicity_G->GetXaxis()->CenterTitle();
  HSiPMLow_Multiplicity_G->GetYaxis()->CenterTitle();

  HSiPM_Multiplicities_G = new TH2I(("HSiPM_Multiplicities_G"), ("HSiPM_Multiplicities_G"), 10, 0, 10, 10, 0, 10);
  HSiPM_Multiplicities_G->GetXaxis()->SetTitle("Multiplicity High");
  HSiPM_Multiplicities_G->GetYaxis()->SetTitle("Multiplicity Low");
  HSiPM_Multiplicities_G->GetXaxis()->CenterTitle();
  HSiPM_Multiplicities_G->GetYaxis()->CenterTitle();
  HSiPM_Multiplicities_G->SetDrawOption("COLZ");
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
  }

  for (size_t i = 0; i < BETA_SIZE; ++i)
  {
    HSiPMHigh_SUMCorrelation_C[i] = NULL;
    HSiPMLow_SUMCorrelation_C[i] = NULL;
    for (size_t j = 0; j < BETA_SIZE; ++j)
    {
      HSiPMHigh_Correlation_C[i][j] = NULL;
      HSiPMLow_Correlation_C[i][j] = NULL;
    }
  }

  for (size_t i = 0; i < detectorNum; ++i)
  {
    if (IsDetectorSiliStrip(i))
    {
      HStrip_Channel_C[i] = new TH1I(("HStrip_Channel_C_" + detectorName[i]).c_str(), ("HStrip_Channel_C_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      HStrip_Channel_C[i]->GetXaxis()->SetTitle("Strips [Channel]");
      HStrip_Channel_C[i]->GetYaxis()->SetTitle("Counts");
      HStrip_Channel_C[i]->GetXaxis()->CenterTitle();
      HStrip_Channel_C[i]->GetYaxis()->CenterTitle();

      HStripRear_Time_C[i] = new TH1I(("HStripRear_Time_C_" + detectorName[i]).c_str(), ("HStripRear_Time_C_" + detectorName[i]).c_str(), winSiliN, winSiliMin, winSiliMax);
      HStripRear_Time_C[i]->GetXaxis()->SetTitle("Time [ns]");
      HStripRear_Time_C[i]->GetYaxis()->SetTitle("Counts");
      HStripRear_Time_C[i]->GetXaxis()->CenterTitle();
      HStripRear_Time_C[i]->GetYaxis()->CenterTitle();

      HStripRear_Channel_C[i] = new TH2I(("HStripRear_Channel_C_" + detectorName[i]).c_str(), ("HStripRear_Channel_C_" + detectorName[i]).c_str(), eSiliN / 10, eSiliMin, eSiliMax, eSiliN / 10, eSiliMin, eSiliMax);
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
      HStripsMultiplicity_C[i] = new TH2I(("HStripsMultiplicity_C_" + detectorName[i]).c_str(), ("HStripsMultiplicity_C_" + detectorName[i]).c_str(), 6, 0, 6, 6, 0, 6);
      HStripsMultiplicity_C[i]->GetXaxis()->SetTitle("Strip A");
      HStripsMultiplicity_C[i]->GetYaxis()->SetTitle("Strip B");
      HStripsMultiplicity_C[i]->GetXaxis()->CenterTitle();
      HStripsMultiplicity_C[i]->GetYaxis()->CenterTitle();
      HStripsMultiplicity_C[i]->SetDrawOption("COLZ");

      HRear_Channel_C[i] = new TH1I(("HRear_Channel_C_" + detectorName[i]).c_str(), ("HRear_Channel_C_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      HRear_Channel_C[i]->GetXaxis()->SetTitle("Rear [Channel]");
      HRear_Channel_C[i]->GetYaxis()->SetTitle("Counts");
      HRear_Channel_C[i]->GetXaxis()->CenterTitle();
      HRear_Channel_C[i]->GetYaxis()->CenterTitle();
    }

    if (IsDetectorBetaHigh(i))
    {
      HSiPMHigh_Channel_C[i] = new TH1I(("HSiPMHigh_Channel_C_" + detectorName[i]).c_str(), ("HSiPMHigh_Channel_C_" + detectorName[i]).c_str(), eHighN, eHighMin, eHighMax);
      HSiPMHigh_Channel_C[i]->GetXaxis()->SetTitle("SiPM High [Channel]");
      HSiPMHigh_Channel_C[i]->GetYaxis()->SetTitle("Counts");
      HSiPMHigh_Channel_C[i]->GetXaxis()->CenterTitle();
      HSiPMHigh_Channel_C[i]->GetYaxis()->CenterTitle();

      HSiPMHighRear_Time_C[i] = new TH1I(("HSiPMHighRear_Time_C_" + detectorName[i]).c_str(), ("HSiPMHighRear_Time_C_" + detectorName[i]).c_str(), winHighN, winHighMin, winHighMax);
      HSiPMHighRear_Time_C[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMHighRear_Time_C[i]->GetYaxis()->SetTitle("Counts");
      HSiPMHighRear_Time_C[i]->GetXaxis()->CenterTitle();
      HSiPMHighRear_Time_C[i]->GetYaxis()->CenterTitle();

      HSiPMHighRear_TimeChannel_C[i] = new TH2I(("HSiPMHighRear_TimeChannel_C_" + detectorName[i]).c_str(), ("HSiPMHighRear_TimeChannel_C_" + detectorName[i]).c_str(), winHighN, winHighMin, winHighMax, eHighN, eHighMin, eHighMax);
      HSiPMHighRear_TimeChannel_C[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMHighRear_TimeChannel_C[i]->GetYaxis()->SetTitle("SiPM High [Channel]");
      HSiPMHighRear_TimeChannel_C[i]->GetXaxis()->CenterTitle();
      HSiPMHighRear_TimeChannel_C[i]->GetYaxis()->CenterTitle();
      HSiPMHighRear_TimeChannel_C[i]->SetDrawOption("COLZ");

      HSiPMHigh_SUMCorrelation_C[i] = new TH2I(("HSiPMHigh_SUMCorrelation_C_" + detectorName[i]).c_str(), ("HSiPMHigh_SUMCorrelation_C_" + detectorName[i]).c_str(), eHighN / 10, eHighMin, eHighMax, eHighN/10, eHighMin, eHighMax);
      HSiPMHigh_SUMCorrelation_C[i]->GetXaxis()->SetTitle((detectorName[i] + " [Channel]").c_str());
      HSiPMHigh_SUMCorrelation_C[i]->GetYaxis()->SetTitle("All High SiPMs [Channel]");
      HSiPMHigh_SUMCorrelation_C[i]->GetXaxis()->CenterTitle();
      HSiPMHigh_SUMCorrelation_C[i]->GetYaxis()->CenterTitle();
      HSiPMHigh_SUMCorrelation_C[i]->SetDrawOption("COLZ");
    }

    if (IsDetectorBetaLow(i))
    {
      HSiPMLow_Channel_C[i] = new TH1I(("HSiPMLow_Channel_C_" + detectorName[i]).c_str(), ("HSiPMLow_Channel_C_" + detectorName[i]).c_str(), eLowN, eLowMin, eLowMax);
      HSiPMLow_Channel_C[i]->GetXaxis()->SetTitle("SiPM Low [Channel]");
      HSiPMLow_Channel_C[i]->GetYaxis()->SetTitle("Counts");
      HSiPMLow_Channel_C[i]->GetXaxis()->CenterTitle();
      HSiPMLow_Channel_C[i]->GetYaxis()->CenterTitle();

      HSiPMLowRear_Time_C[i] = new TH1I(("HSiPMLowRear_Time_C_" + detectorName[i]).c_str(), ("HSiPMLowRear_Time_C_" + detectorName[i]).c_str(), winLowN, winLowMin, winLowMax);
      HSiPMLowRear_Time_C[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMLowRear_Time_C[i]->GetYaxis()->SetTitle("Counts");
      HSiPMLowRear_Time_C[i]->GetXaxis()->CenterTitle();
      HSiPMLowRear_Time_C[i]->GetYaxis()->CenterTitle();

      HSiPMLowRear_TimeChannel_C[i] = new TH2I(("HSiPMLowRear_TimeChannel_C_" + detectorName[i]).c_str(), ("HSiPMLowRear_TimeChannel_C_" + detectorName[i]).c_str(), winLowN, winLowMin, winLowMax, eLowN, eLowMin, eLowMax);
      HSiPMLowRear_TimeChannel_C[i]->GetXaxis()->SetTitle("Time [ns]");
      HSiPMLowRear_TimeChannel_C[i]->GetYaxis()->SetTitle("SiPM Low [Channel]");
      HSiPMLowRear_TimeChannel_C[i]->GetXaxis()->CenterTitle();
      HSiPMLowRear_TimeChannel_C[i]->GetYaxis()->CenterTitle();
      HSiPMLowRear_TimeChannel_C[i]->SetDrawOption("COLZ");

      HSiPMLow_SUMCorrelation_C[i] = new TH2I(("HSiPMLow_SUMCorrelation_C_" + detectorName[i]).c_str(), ("HSiPMLow_SUMCorrelation_C_" + detectorName[i]).c_str(), eLowN / 10, eLowMin, eLowMax, eLowN /10 , eLowMin, eLowMax);
      HSiPMLow_SUMCorrelation_C[i]->GetXaxis()->SetTitle((detectorName[i] + " [Channel]").c_str());
      HSiPMLow_SUMCorrelation_C[i]->GetYaxis()->SetTitle("All Low SiPMs [Channel]");
      HSiPMLow_SUMCorrelation_C[i]->GetXaxis()->CenterTitle();
      HSiPMLow_SUMCorrelation_C[i]->GetYaxis()->CenterTitle();
      HSiPMLow_SUMCorrelation_C[i]->SetDrawOption("COLZ");
    }
  }
  HSiPMHigh_Multiplicity_C = new TH1I(("HSiPMHigh_Multiplicity_C"), ("HSiPMHigh_Multiplicity_C"), 10, 0, 10);
  HSiPMHigh_Multiplicity_C->GetXaxis()->SetTitle("Multiplicity");
  HSiPMHigh_Multiplicity_C->GetYaxis()->SetTitle("Counts");
  HSiPMHigh_Multiplicity_C->GetXaxis()->CenterTitle();
  HSiPMHigh_Multiplicity_C->GetYaxis()->CenterTitle();

  for (int j = 1; j <= 9; ++j)
  {
    for (int k = j + 1; k <= 9; ++k)
    {
      std::string histName = "HBetaHiCorrelation" + std::to_string(j) + "_" + std::to_string(k);
      std::string histTitle = "BetaHigh Correlation " + std::to_string(j) + "/" + std::to_string(k);
      HSiPMHigh_Correlation_C[j - 1][k - 1] = new TH2I(histName.c_str(), histTitle.c_str(), eHighN / 10, eHighMin, eHighMax, eHighN / 10, eHighMin, eHighMax);

      HSiPMHigh_Correlation_C[j - 1][k - 1]->GetXaxis()->SetTitle(("Beta " + std::to_string(j) + "[Channel]").c_str());
      HSiPMHigh_Correlation_C[j - 1][k - 1]->GetYaxis()->SetTitle(("Beta " + std::to_string(k) + "[Channel]").c_str());
      HSiPMHigh_Correlation_C[j - 1][k - 1]->GetXaxis()->CenterTitle();
      HSiPMHigh_Correlation_C[j - 1][k - 1]->GetYaxis()->CenterTitle();
    }
  }

  HSiPMLow_Multiplicity_C = new TH1I(("HSiPMLow_Multiplicity_C"), ("HSiPMLow_Multiplicity_C"), 10, 0, 10);
  HSiPMLow_Multiplicity_C->GetXaxis()->SetTitle("Multiplicity");
  HSiPMLow_Multiplicity_C->GetYaxis()->SetTitle("Counts");
  HSiPMLow_Multiplicity_C->GetXaxis()->CenterTitle();
  HSiPMLow_Multiplicity_C->GetYaxis()->CenterTitle();

  for (int j = 1; j <= 9; ++j)
  {
    for (int k = j + 1; k <= 9; ++k)
    {
      std::string histName = "HBetaLowCorrelation" + std::to_string(j) + "_" + std::to_string(k);
      std::string histTitle = "BetaLow Correlation " + std::to_string(j) + "/" + std::to_string(k);
      HSiPMLow_Correlation_C[j - 1][k - 1] = new TH2I(histName.c_str(), histTitle.c_str(), eLowN / 10, eLowMin, eLowMax, eLowN / 10, eLowMin, eLowMax);

      HSiPMLow_Correlation_C[j - 1][k - 1]->GetXaxis()->SetTitle(("Beta " + std::to_string(j) + "[Channel]").c_str());
      HSiPMLow_Correlation_C[j - 1][k - 1]->GetYaxis()->SetTitle(("Beta " + std::to_string(k) + "[Channel]").c_str());
      HSiPMLow_Correlation_C[j - 1][k - 1]->GetXaxis()->CenterTitle();
      HSiPMLow_Correlation_C[j - 1][k - 1]->GetYaxis()->CenterTitle();
    }
  }

  HSiPM_Multiplicities_C = new TH2I(("HSiPM_Multiplicities_C"), ("HSiPM_Multiplicities_C"), 10, 0, 10, 10, 0, 10);
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
      Tree_Silicons[i]->Branch("Channel", &Tree_Channel, "Channel/I");
    }

    if (IsDetectorBeta(i))
    {
      Tree_SiPMs[i] = new TTree(("Tree_SiPM_" + detectorName[i]).c_str(), ("Tree_SiPM_" + detectorName[i]).c_str());
      Tree_SiPMs[i]->Branch("Channel", &Tree_Channel, "Channel/I");
    }
  }

  return 0;
}

inline int WriteHistograms_Grouped()
{
  File_Grouped->cd();
  TDirectory *Strip_Channel_G = File_Grouped->mkdir("Strip_Channel");
  TDirectory *Strip_Time_G = File_Grouped->mkdir("Strip_Time");
  TDirectory *Rear_Strip_Channel_G = File_Grouped->mkdir("Rear_Strip_Channel");
  TDirectory *fracRear_Strip_Channel_G = File_Grouped->mkdir("Frac_Rear_Strip_Channel");
  TDirectory *Strip_Multiplicities_G = File_Grouped->mkdir("Strip_Multiplicities");

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
    if (IsDetectorSiliStrip(i))
    {
      Strip_Channel_G->cd();
      HStrip_Channel_G[i]->Write();
      delete HStrip_Channel_G[i];
      Strip_Time_G->cd();
      HStripRear_Time_G[i]->Write();
      delete HStripRear_Time_G[i];
      Rear_Strip_Channel_G->cd();
      HStripRear_Channel_G[i]->Write();
      delete HStripRear_Channel_G[i];
      fracRear_Strip_Channel_G->cd();
      HFracRearStripChannel_G[i]->Write();
      delete HFracRearStripChannel_G[i];
    }

    if (IsDetectorSiliBack(i))
    {
      Strip_Multiplicities_G->cd();
      HStripsMultiplicity_G[i]->Write();
      delete HStripsMultiplicity_G[i];
    }

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
  TDirectory *SiPMHigh_Correlation_C = File_Cleaned->mkdir("SiPMHigh_Correlation");
  TDirectory *SiPMHigh_SUMCorrelation_C = File_Cleaned->mkdir("SiPMHigh_SUMCorrelation");

  TDirectory *SiPMLow_Channel_C = File_Cleaned->mkdir("SiPMLow_Channel");
  TDirectory *SiPMLow_Time_C = File_Cleaned->mkdir("SiPMLow_Time");
  TDirectory *SiPMLowRear_TimeChannel_C = File_Cleaned->mkdir("SiPMLowRear_TimeChannel");
  TDirectory *SiPMLow_Multiplicity_C = File_Cleaned->mkdir("SiPMLow_Multiplicity");
  TDirectory *SiPMLow_Correlation_C = File_Cleaned->mkdir("SiPMLow_Correlation");
  TDirectory *SiPMLow_SUMCorrelation_C = File_Cleaned->mkdir("SiPMLow_SUMCorrelation");

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
      SiPMHigh_SUMCorrelation_C->cd();
      HSiPMHigh_SUMCorrelation_C[i]->Write(); // delete HSiPMHigh_SUMCorrelation_C[i];
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
      SiPMLow_SUMCorrelation_C->cd();
      HSiPMLow_SUMCorrelation_C[i]->Write(); // delete HSiPMLow_SUMCorrelation_C[i];
    }
  }

  SiPMHigh_Multiplicity_C->cd();
  HSiPMHigh_Multiplicity_C->Write();
  delete HSiPMHigh_Multiplicity_C;

  SiPMHigh_Correlation_C->cd();
  ///// TCANVAS
  for (int j = 1; j <= 9; ++j)
  {
    CanvasHigh = new TCanvas(("High" + to_string(j)).c_str(), ("High" + to_string(j)).c_str(), 1600, 800);
    CanvasHigh->Divide(3, 3);

    for (int k = 1; k <= 9; ++k)
    {

      if (HSiPMHigh_Correlation_C[j - 1][k - 1] != NULL)
      {
        CanvasHigh->cd(k);
        HSiPMHigh_Correlation_C[j - 1][k - 1]->Draw("COLZ");
      }

      else if (HSiPMHigh_Correlation_C[k - 1][j - 1] != NULL)
      {
        CanvasHigh->cd(k);
        HSiPMHigh_Correlation_C[k - 1][j - 1]->Draw("COLZ");
      }
    }

    CanvasHigh->Write();
  }
  for (int j = 1; j <= 9; ++j)
  {
    for (int k = j + 1; k <= 9; ++k)
    {
      delete HSiPMHigh_Correlation_C[j - 1][k - 1];
    }
  }

  SiPMLow_Multiplicity_C->cd();
  HSiPMLow_Multiplicity_C->Write();
  delete HSiPMLow_Multiplicity_C;

  SiPMLow_Correlation_C->cd();
  ///// TCANVAS
  for (int j = 1; j <= 9; ++j)
  {
    CanvasLow = new TCanvas(("Low" + to_string(j)).c_str(), ("Low" + to_string(j)).c_str(), 1600, 800);
    CanvasLow->Divide(3, 3);

    for (int k = 1; k <= 9; ++k)
    {

      if (HSiPMLow_Correlation_C[j - 1][k - 1] != NULL)
      {
        CanvasLow->cd(k);
        HSiPMLow_Correlation_C[j - 1][k - 1]->Draw("COLZ");
      }

      else if (HSiPMLow_Correlation_C[k - 1][j - 1] != NULL)
      {
        CanvasLow->cd(j);
        HSiPMLow_Correlation_C[k - 1][j - 1]->Draw("COLZ");
      }
    }

    CanvasLow->Write();
  }
  for (int j = 1; j <= 9; ++j)
  {
    for (int k = j + 1; k <= 9; ++k)
    {
      delete HSiPMLow_Correlation_C[j - 1][k - 1];
    }
  }

  ////TCANVAS

  SiPM_Multiplicities_C->cd();
  HSiPM_Multiplicities_C->Write();
  delete HSiPM_Multiplicities_C;
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

    if (IsDetectorBeta(i))
    {
      Tree_SiPMs[i]->Write();
      delete Tree_SiPMs[i];
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

void SavingData_Grouped(vector<Signal> SiliconRear, vector<Signal> SiPM_High, vector<Signal> SiPM_Low)
{
  Signal TrigSignal = SiliconRear[1];
  Signal Silicon = SiliconRear[0];
  /// Saving in Hist

  if (Silicon.isValid)
  {
    HStrip_Channel_G[Silicon.Label]->Fill(Silicon.Channel);
    HStripRear_Time_G[Silicon.Label]->Fill(Silicon.Time - TrigSignal.Time);
    HStripRear_Channel_G[Silicon.Label]->Fill(TrigSignal.Channel, Silicon.Channel);
    HFracRearStripChannel_G[Silicon.Label]->Fill(static_cast<double>(Silicon.Channel) / TrigSignal.Channel);

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
  if (Silicon.isValid)
  {

    Event++;
    Tree_Event = Event;
    Tree_Silicon = SiliconRear;
    Tree_SiPMHigh = SiPM_High;
    Tree_SiPMLow = SiPM_Low;
    Tree_Grouped->Fill();
  }
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

      Total = 0;
      for (auto it2 = SiPM_High.begin(); it2 != SiPM_High.end(); ++it2)
      {
        if (it->Label != it2->Label)
        {
          if (it->Label < it2->Label)
          {
            HSiPMHigh_Correlation_C[GetDetectorChannel(it->Label) - 1][GetDetectorChannel(it2->Label) - 1]->Fill(it->Channel, it2->Channel);
          }
          else
            HSiPMHigh_Correlation_C[GetDetectorChannel(it2->Label) - 1][GetDetectorChannel(it->Label) - 1]->Fill(it2->Channel, it->Channel);
        }

        Total += it2->Channel;
      }
      HSiPMHigh_SUMCorrelation_C[it->Label]->Fill(it->Channel, Total/SiPM_High.size());
    }
    for (auto it = SiPM_Low.begin(); it != SiPM_Low.end(); ++it)
    {
      HSiPMLow_Channel_C[it->Label]->Fill(it->Channel);
      HSiPMLowRear_Time_C[it->Label]->Fill(it->Time - TrigSignal.Time);
      HSiPMLowRear_TimeChannel_C[it->Label]->Fill(it->Time - TrigSignal.Time, it->Channel);

      Total = 0;
      for (auto it2 = SiPM_Low.begin(); it2 != SiPM_Low.end(); ++it2)
      {
        if (it->Label != it2->Label)
        {
          if (it->Label < it2->Label)
          {
            HSiPMLow_Correlation_C[GetDetectorChannel(it->Label) - 1][GetDetectorChannel(it2->Label) - 1]->Fill(it->Channel, it2->Channel);
          }
          else
            HSiPMLow_Correlation_C[GetDetectorChannel(it2->Label) - 1][GetDetectorChannel(it->Label) - 1]->Fill(it2->Channel, it->Channel);
        }
        Total += it2->Channel;
      }
      HSiPMLow_SUMCorrelation_C[GetDetectorChannel(it->Label) - 1]->Fill(it->Channel, Total/SiPM_Low.size());
    }

    HSiPMHigh_Multiplicity_C->Fill(SiPM_High.size());
    HSiPMLow_Multiplicity_C->Fill(SiPM_Low.size());
    HSiPM_Multiplicities_C->Fill(SiPM_High.size(), SiPM_Low.size());
  }

  else
  {
    //////TO DO
    // Analyse no-coincidence
  }
}

Signal ProcessSilicon(vector<Signal> Signal_Vector, Signal Trigger)
{
  vector<Signal> Output_Vector = vector<Signal>();
  /// Search in Signal Vector between 1 and n
  for (auto current = Signal_Vector.begin(); current != Signal_Vector.end(); ++current)
  {
    if (current->Time - Trigger.Time < winSiliMax && Trigger.Time - current->Time > winSiliMin && IsSameSiliDetector(Trigger.Label, current->Label) && Trigger != *current)
    {
      Output_Vector.push_back(*current);
    }
  }

  if (Output_Vector.size() != 1)
  {
    if (Output_Vector.size() == 2)
    {
      if (Output_Vector[0].Label < Output_Vector[1].Label)
        HStripsMultiplicity_G[Trigger.Label]->Fill(Output_Vector[0].Label, Output_Vector[1].Label);
      else
        HStripsMultiplicity_G[Trigger.Label]->Fill(Output_Vector[1].Label, Output_Vector[0].Label);
    }

    return Signal();
  }

  return Output_Vector[0];
}

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

void SearchForCoincidence(TTreeReader *Reader, TTreeReaderValue<double> &Times, TTreeReaderValue<int> &Labels, TTreeReaderValue<int> &Channels)
{
  if (Verbose > 1)
    cout << "### Starting a Coincidence group on " << detectorName[*Labels] << endl;

  // Vector save
  vector<Signal> SiPMSignalsHigh;
  vector<Signal> SiPMSignalsLow;
  vector<Signal> SiSignals;

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

  while (*Times <= TrigSignal.Time + winTotalMax)
  {
    Signal current = Signal(*Labels, *Times, *Channels); /// Current Signal

    if (IsDetectorSili(*Labels))
    {
      SiSignals.push_back(current);
    }
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

  ////Raw Save

  /// Selecting Signals in time window
  Signal Silicon = ProcessSilicon(SiSignals, TrigSignal);              /// Time Window + single silicon signal on same detector
  vector<Signal> SiPM_High = ProcessSiPM(SiPMSignalsHigh, TrigSignal); /// Time Window
  vector<Signal> SiPM_Low = ProcessSiPM(SiPMSignalsLow, TrigSignal);   /// Time Window

  vector<Signal> SiliconRear = {Silicon, TrigSignal};

  SavingData_Grouped(SiliconRear, SiPM_High, SiPM_Low);
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