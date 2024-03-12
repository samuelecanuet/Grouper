#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <boost/filesystem.hpp>

#include "TFile.h"
#include "TCanvas.h"
#include <TKey.h>
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
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
namespace fs = boost::filesystem;

map<string, vector<string>> fileMap;
string Catcher;
vector<int> RunNumbers; 
vector<int> runs;
string baseFileName;
double detectorCalib[SIGNAL_MAX];
double detectorCalib_Fermi[SIGNAL_MAX];
double detectorCalib_GT[SIGNAL_MAX];
double detectorMatching[SIGNAL_MAX];
vector<pair<string, string>> file_string_Time;
vector<pair<double, double>> file_Time;

TFile *Merged_File;
TTree *Tree_Merged;
vector<Signal> Tree_Merged_Silicon;
vector<Signal> Tree_Merged_SiPM;

TFile *Matched_File;
TFile *Cleaned_File;
TTree *Tree_Read;
TH1D *HSilicon_Channel_Unmatched[SIGNAL_MAX];
TH1D *HSilicon_Channel_Matched[SIGNAL_MAX];
TH1D *HSilicon[SIGNAL_MAX];
TH1D *HSilicon_coinc[SIGNAL_MAX];
TH1D *HSilicon_no_coinc[SIGNAL_MAX];

TH1D *HSiPM[BETA_SIZE+1];
TH1D *HSiPM_F[BETA_SIZE+1];
TH1D *HSiPM_GT[BETA_SIZE+1];

TH2D* HSiPMRun[BETA_SIZE+1];

TH2D *HSiliconRun_Channel_Unmatched[SIGNAL_MAX];
TH2D *HSiliconRun_Channel_Matched[SIGNAL_MAX];
TH2D *HSiliconRun[SIGNAL_MAX];

TH1D *HSilicon_Channel_Unmatched_temp[SIGNAL_MAX];
TH1D *HSilicon_Channel_Matched_temp[SIGNAL_MAX];
TH1D *HSilicon_temp[SIGNAL_MAX];
TH1D *to_save_unmatched;
TH1D *to_save_matched;

TF1 *linear = new TF1("linear", "[0]*x");


int Verbose = 0;

bool Calibration = false;

int counter;

double bin_time;
int n_bin_time;
double start_time;
double stop_time;


////////////////////////////////////

void ProgressBar(ULong64_t cEntry, ULong64_t TotalEntries, clock_t start, clock_t Current)
{
  if (cEntry % 100000 == 0 && cEntry > 2 * 100000)
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
         << "           " << flush;
  }
}


double Convert_DatetoTime(string datestring, string refstring)
{
    std::tm date = {}, referenceDate = {};
    istringstream dateStream(datestring);
    dateStream >> get_time(&date, "%d-%m-%Y %H:%M:%S");

    istringstream referenceDateStream(refstring);
    referenceDateStream >> get_time(&referenceDate, "%d-%m-%Y %H:%M:%S");

    time_t timeInSeconds = mktime(&date);
    time_t referenceTimeInSeconds = mktime(&referenceDate);

    double diffInSeconds = difftime(timeInSeconds, referenceTimeInSeconds);

    double diffInHours = diffInSeconds / 3600.0; 

    return diffInHours;
}

double GetOnlyHour(string refstring)
{
    std::tm referenceDate = {};
    istringstream referenceDateStream(refstring);
    referenceDateStream >> get_time(&referenceDate, "%d-%m-%Y %H:%M:%S");

    time_t referenceTimeInSeconds = mktime(&referenceDate);

    double diffInHours = referenceDate.tm_hour + referenceDate.tm_min/60.0 + referenceDate.tm_sec/3600.0; 

    return diffInHours;
}

void InitRunHistograms()
{
  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      HSilicon_Channel_Unmatched_temp[i] = new TH1D(("HSilicon_Channel_Unmatched_temp_" + detectorName[i]).c_str(), ("HSilicon_Channel_Unmatched_temp_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      HSilicon_Channel_Matched_temp[i] = new TH1D(("HSilicon_Channel_Matched_temp_" + detectorName[i]).c_str(), ("HSilicon_Channel_Matched_temp_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      HSilicon_temp[i] = new TH1D(("HStrip_temp_" + detectorName[i]).c_str(), ("HStrip_temp_" + detectorName[i]).c_str(), 10000, 0, 10000);
    }
  }
}

void InitHistograms()
{
  //// BINING TIME /////////////////////////
  bin_time = std::numeric_limits<double>::max();
  for (const auto& p : file_Time) 
  {
      double diff = std::abs(p.first - p.second);
      bin_time = std::min(bin_time, diff);
  }

  double t_min = file_Time[0].first;
  int t_max = (int)file_Time[file_Time.size()-1].second + 1;
  n_bin_time = (int)(t_max - t_min) / bin_time;
  //////////////////////////////////////////

  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      HSiliconRun_Channel_Unmatched[i] = new TH2D(("HSiliconRun_Channel_Unmatched_" + detectorName[i]).c_str(), ("HSiliconRun_Channel_Unmatched_" + detectorName[i]).c_str(), t_max, t_min, t_max,  eSiliN, eSiliMin, eSiliMax);
      HSiliconRun_Channel_Unmatched[i]->GetXaxis()->SetTitle("Time (h)");
      HSiliconRun_Channel_Unmatched[i]->GetYaxis()->SetTitle("Channel");
      HSiliconRun_Channel_Unmatched[i]->GetXaxis()->CenterTitle();
      HSiliconRun_Channel_Unmatched[i]->GetYaxis()->CenterTitle();

      HSiliconRun_Channel_Matched[i] = new TH2D(("HSiliconRun_Channel_Matched_" + detectorName[i]).c_str(), ("HSiliconRun_Channel_Matched_" + detectorName[i]).c_str(), t_max, t_min, t_max,  eSiliN, eSiliMin, eSiliMax);
      HSiliconRun_Channel_Matched[i]->GetXaxis()->SetTitle("Time (h)");
      HSiliconRun_Channel_Matched[i]->GetYaxis()->SetTitle("Channel");
      HSiliconRun_Channel_Matched[i]->GetXaxis()->CenterTitle();
      HSiliconRun_Channel_Matched[i]->GetYaxis()->CenterTitle();

      HSiliconRun[i] = new TH2D(("HSiliconRun"+ detectorName[i]).c_str(), ("HSiliconRun"+ detectorName[i]).c_str(), t_max, t_min, t_max, 5000, 0, 10000);
      HSiliconRun[i]->GetXaxis()->SetTitle("Time (h)");
      HSiliconRun[i]->GetYaxis()->SetTitle("Energy [keV]");
      HSiliconRun[i]->GetXaxis()->CenterTitle();
      HSiliconRun[i]->GetYaxis()->CenterTitle();

      HSilicon_Channel_Unmatched[i] = new TH1D(("HSilicon_Channel_Unmatched_" + detectorName[i]).c_str(), ("HSilicon_Channel_Unmatched_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      HSilicon_Channel_Unmatched[i]->GetXaxis()->SetTitle("Channel");
      HSilicon_Channel_Unmatched[i]->GetYaxis()->SetTitle("Counts");
      HSilicon_Channel_Unmatched[i]->GetXaxis()->CenterTitle();
      HSilicon_Channel_Unmatched[i]->GetYaxis()->CenterTitle();

      HSilicon_Channel_Matched[i] = new TH1D(("HSilicon_Channel_Matched_" + detectorName[i]).c_str(), ("HSilicon_Channel_Matched_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      HSilicon_Channel_Matched[i]->GetXaxis()->SetTitle("Channel");
      HSilicon_Channel_Matched[i]->GetYaxis()->SetTitle("Counts");
      HSilicon_Channel_Matched[i]->GetXaxis()->CenterTitle();
      HSilicon_Channel_Matched[i]->GetYaxis()->CenterTitle();

      HSilicon[i] = new TH1D(("HStrip_" + detectorName[i]).c_str(), ("HStrip_" + detectorName[i]).c_str(), 30000, 0, 10000);
      HSilicon[i]->GetXaxis()->SetTitle("Energy [keV]");
      HSilicon[i]->GetYaxis()->SetTitle("Counts");
      HSilicon[i]->GetXaxis()->CenterTitle();
      HSilicon[i]->GetYaxis()->CenterTitle();

      HSilicon_coinc[i] = new TH1D(("HStrip_coinc_" + detectorName[i]).c_str(), ("HStrip_coinc_" + detectorName[i]).c_str(), 30000, 0, 10000);
      HSilicon_coinc[i]->GetXaxis()->SetTitle("Energy [keV]");
      HSilicon_coinc[i]->GetYaxis()->SetTitle("Counts");
      HSilicon_coinc[i]->GetXaxis()->CenterTitle();
      HSilicon_coinc[i]->GetYaxis()->CenterTitle();

      HSilicon_no_coinc[i] = new TH1D(("HStrip_no_coinc_" + detectorName[i]).c_str(), ("HStrip_no_coinc_" + detectorName[i]).c_str(), 30000, 0, 10000);
      HSilicon_no_coinc[i]->GetXaxis()->SetTitle("Energy [keV]");
      HSilicon_no_coinc[i]->GetYaxis()->SetTitle("Counts");
      HSilicon_no_coinc[i]->GetXaxis()->CenterTitle();
      HSilicon_no_coinc[i]->GetYaxis()->CenterTitle();
    }
  }
    for (int i = 1; i <= BETA_SIZE; i++)
    {
      HSiPMRun[GetDetectorChannel(i)] = new TH2D(("HSiPMRun_M"+to_string(i)).c_str(), ("HSiPMRun_M"+to_string(i)).c_str(), t_max, t_min, t_max, eLowN/10, eLowMin, eLowMax*8);
      HSiPMRun[GetDetectorChannel(i)]->GetXaxis()->SetTitle("Time (h)");
      HSiPMRun[GetDetectorChannel(i)]->GetYaxis()->SetTitle("Channel");
      HSiPMRun[GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
      HSiPMRun[GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

      HSiPM[GetDetectorChannel(i)] = new TH1D(("HSiPM_M"+to_string(i)).c_str(), ("HSiPM_M"+to_string(i)).c_str(), eLowN/10, eLowMin, eLowMax*5);
      HSiPM[GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
      HSiPM[GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
      HSiPM[GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
      HSiPM[GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

      HSiPM_F[GetDetectorChannel(i)] = new TH1D(("HSiPM_F_M"+to_string(i)).c_str(), ("HSiPM_F_M"+to_string(i)).c_str(), eLowN/10, eLowMin, eLowMax*5);
      HSiPM_F[GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
      HSiPM_F[GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
      HSiPM_F[GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
      HSiPM_F[GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

      HSiPM_GT[GetDetectorChannel(i)] = new TH1D(("HSiPM_GT_M"+to_string(i)).c_str(), ("HSiPM_GT_M"+to_string(i)).c_str(), eLowN/10, eLowMin, eLowMax*5);
      HSiPM_GT[GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
      HSiPM_GT[GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
      HSiPM_GT[GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
      HSiPM_GT[GetDetectorChannel(i)]->GetYaxis()->CenterTitle();
    }
  
  
}

void WriteHistograms()
{
  Merged_File->cd();
  TDirectory *dir_HSilicon_Channel_Unmatched = Merged_File->mkdir("Silicon_Unmatched");
  TDirectory *dir_HSilicon_Channel_Matched = Merged_File->mkdir("Silicon_Matched");
  TDirectory *dir_HSilicon = Merged_File->mkdir("Silicon_Calibrated");

  TDirectory *dir_HSiPM = Merged_File->mkdir("SiPM");

  double first_time = GetOnlyHour(file_string_Time[0].first);
  cout<<GetOnlyHour(file_string_Time[0].first)<<endl;
  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      dir_HSilicon_Channel_Unmatched->cd();
      HSilicon_Channel_Unmatched[i]->Write();
      HSiliconRun_Channel_Unmatched[i]->Write();

      dir_HSilicon_Channel_Matched->cd();
      HSilicon_Channel_Matched[i]->Write();

      TCanvas *c1 = new TCanvas("c", "c", 800, 600);
      c1->cd();
      TAxis *xaxis = HSiliconRun_Channel_Matched[i]->GetXaxis();     
      for (int bin = 1; bin <= xaxis->GetNbins(); ++bin) {
        int hourmodulo = (first_time + bin) / 24 - static_cast<int>(first_time+bin/24);;
        xaxis->SetBinLabel(bin, Form("%d", hourmodulo));
      }
      HSiliconRun_Channel_Matched[i]->GetXaxis()->SetTitle("Time (h modulo 24)");
      HSiliconRun_Channel_Matched[i]->Draw("HIST");
      c1->Write();
      delete c1;
      
      dir_HSilicon->cd();
      HSilicon[i]->Write();
      HSiliconRun[i]->Write();
      TCanvas *c = new TCanvas((detectorName[i]).c_str(), (detectorName[i]).c_str(), 800, 600);
      c->cd();
      HSilicon[i]->SetLineColor(kBlack);
      HSilicon[i]->Draw("HIST");
      HSilicon_coinc[i]->SetLineColor(kRed);
      HSilicon_coinc[i]->Draw("SAME");
      HSilicon_no_coinc[i]->SetLineColor(kBlue);
      HSilicon_no_coinc[i]->Draw("SAME");
      c->Write();
      delete c;
    }
  }

  /////Multiplicity
  dir_HSiPM->cd();
  for (int i = 0; i <= BETA_SIZE; i++)
  {
    HSiPM[GetDetectorChannel(i)]->Write();
    HSiPM_F[GetDetectorChannel(i)]->Write();
    HSiPM_GT[GetDetectorChannel(i)]->Write();
  }
}

void InitFiles()
{

  if (Catcher == "Al")
  {
    RunNumbers = {11, 12, 13, 14, 15, 16, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 36, 37, 38, 39, 40};
  }
  else if (Catcher == "Mylar")
  {
    RunNumbers = {49, 20, 51, 52, 53, 55, 56, 57, 58, 59, 60, 61, 62};
  }

  for (int run : RunNumbers)
  {
    runs.push_back(run);
  }

}

void InitMatching()
{
  TGraphErrors *MatchingGraph = nullptr;
  TIter next(Matched_File->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *)next()))
  {
    if (strcmp(key->GetClassName(), "TGraphErrors") == 0)
    {
      MatchingGraph = (TGraphErrors *)Matched_File->Get(key->GetName());
      break;
    }
  }
  
  double label;
  double value;
  for (int counter = 1; counter < MatchingGraph->GetN(); counter++)
  {
    MatchingGraph->GetPoint(counter, label, value);
    detectorMatching[static_cast<int>(label)] = value;
  }
}

int MakeCalibration() ///////////////////////////////++++++FIT
{
  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      // HSilicon_Channel_Matched[i]->GetXaxis()->SetRangeUser(detectorRangeMin_Fermi[i], detectorRangeMax_Fermi[i]);
      // double fermi[2] = {HSilicon_Channel_Matched[i]->GetMean(), detectorCalib_Fermi[i]};
      // HSilicon_Channel_Matched[i]->GetXaxis()->SetRangeUser(detectorRangeMin_GT[i], detectorRangeMax_GT[i]);
      // double gt[2] = {HSilicon_Channel_Matched[i]->GetMean(), detectorCalib_GT[i]};

      TGraphErrors* graph = new TGraphErrors();

      HSilicon_Channel_Matched[i]->GetXaxis()->SetRangeUser(40000, 44000);
      graph->SetPoint(0, detectorCalib_Fermi[i], HSilicon_Channel_Matched[i]->GetMean());
      graph->SetPointError(0, 0, HSilicon_Channel_Matched[i]->GetMeanError());

      HSilicon_Channel_Matched[i]->GetXaxis()->SetRangeUser(25000, 27000);
      graph->SetPoint(0, detectorCalib_GT[i], HSilicon_Channel_Matched[i]->GetMean());
      graph->SetPointError(0, 0, HSilicon_Channel_Matched[i]->GetMeanError());

      linear->SetParameters(0, 12);
      graph->Fit("linear");

      detectorCalib[i] = linear->GetParameter(0);
    }
  }
}

int InitCalib()
{

  ///////////////////FOR FERMI
  int error = 0;
  string CalibDir = "../Calib/";
  string CalibFileName = "Sim_32Ar_MylarAl_F.txt";
  ifstream file(CalibDir + CalibFileName);

  for (auto &value : detectorCalib_Fermi)
  {
    value = 0.;
  }

  if (file.is_open())
  {
    GLogMessage("<SAM> Calibration file found");

    Calibration = true;

    string line;
    int i = 0;
    while (getline(file, line))
    {
      istringstream iss(line);
      double E;
      string DetName;
      iss >> DetName >> E;

      int index = 0;
      int strip;
      string stripstr(1, DetName.back());
      if (stripstr == "R")
      {
        strip = 6;
      }
      else
      {
        strip = stoi(stripstr);
      }
      if (DetName.find("Up") != string::npos)
      {
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
          if (GetDetector(i) <= 4 and GetDetectorChannel(i) == strip)
          {
            detectorCalib_Fermi[i] = E;
          }
        }
      }
      else if (DetName.find("Down") != string::npos)
      {
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
          if (GetDetector(i) > 4 and GetDetectorChannel(i) == strip)
          {
            detectorCalib_Fermi[i] = E;
          }
        }
      }
    }
  }
  else
  {
    GLogMessage("<SAM> No Calibration file found");
    error = 1;
  }

  ////////////////////////////////////////
  CalibFileName = "Sim_32Ar_MylarAl_GT1.txt";
  ifstream fileGT(CalibDir + CalibFileName);

  for (auto &value : detectorCalib_Fermi)
  {
    value = 0.;
  }

  if (fileGT.is_open())
  {
    GLogMessage("<SAM> Calibration file found");

    Calibration = true;

    string line;
    int i = 0;
    while (getline(fileGT, line))
    {
      istringstream iss(line);
      double E;
      string DetName;
      iss >> DetName >> E;

      int index = 0;
      int strip;
      string stripstr(1, DetName.back());
      if (stripstr == "R")
      {
        strip = 6;
      }
      else
      {
        strip = stoi(stripstr);
      }
      if (DetName.find("Up") != string::npos)
      {
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
          if (GetDetector(i) <= 4 and GetDetectorChannel(i) == strip)
          {
            detectorCalib_GT[i] = E;
          }
        }
      }
      else if (DetName.find("Down") != string::npos)
      {
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
          if (GetDetector(i) > 4 and GetDetectorChannel(i) == strip)
          {
            detectorCalib_GT[i] = E;
          }
        }
      }
    }
  }
  else
  {
    GLogMessage("<SAM> No Calibration file found");
    error = 1;
  }

  // for (int i = 0; i < SIGNAL_MAX; i++)
  // {
  //   if (IsDetectorBetaHigh(i))
  //   {
  //     detectorCalib[i] = 1390./24000000;
  //   }
  //   if (IsDetectorBetaLow(i))
  //   {
  //     detectorCalib[i] = 6000./9000000;
  //   }
  // }

  // FOR CHECK
  // int index = 0;
  // for (auto &name : detectorName)
  // {
  //   cout << name << "       " << detectorCalib_Fermi[index] << endl;
  //   index++;
  // }

  file.close();
  return error;
}

bool Is_F(double E, int label)
{
  if (E > 3100 && E < 3300)
  {
    return true;
  }
  return false;
}

bool Is_GT(double E, int label)
{
  if (E > 1900 && E < 2200)
  {
    return true;
  }
  return false;
}

void SetTime(TFile* file)
{
  file_string_Time.push_back(GetTime(file));

  double start_time = Convert_DatetoTime(file_string_Time[file_string_Time.size()-1].first, file_string_Time[0].first);
  double stop_time = Convert_DatetoTime(file_string_Time[file_string_Time.size()-1].second, file_string_Time[0].first);
  file_Time.push_back(make_pair(start_time, stop_time));
}


