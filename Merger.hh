#ifndef MERGER_HH
#define MERGER_HH

#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <random>
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
#include "TLine.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRandom.h"

#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"


#include "Signal.h"
#include "Detectors.hh"
#include "/home/local1/Documents/lib/GTools1.0/include/GString.hh"

using namespace std;
namespace fs = boost::filesystem;
using namespace ROOT::Math;
random_device rd;
mt19937 gen(rd());

map<string, vector<string>> fileMap;
string Catcher;
vector<int> RunNumbers;
vector<int> runs;
string baseFileName;
pair<double, double> SiliconCalib[100][SIGNAL_MAX] = {{make_pair(0., 0.)}};
double SiliconCalib_Fermi[SIGNAL_MAX];
double SiliconCalib_GT[SIGNAL_MAX];
pair<double, double> SiPMCalib[100] = {{make_pair(0., 0.)}};
double detectorMatching[SIGNAL_MAX];
pair<int, int> peaks_window_F[SIGNAL_MAX];
pair<int, int> peaks_window_GT[SIGNAL_MAX];
bool RunDetectorSelection[100][SIGNAL_MAX] = {{false}};
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
TH1D *HSilicon_Channel_Matched_CURRENTRUN[SIGNAL_MAX];
TH1D *HSilicon[SIGNAL_MAX];
TH1D *HSilicon_coinc[SIGNAL_MAX];
TH1D *HSilicon_no_coinc[SIGNAL_MAX];

TH1D *HSiPM[BETA_SIZE + 1];
TH1D *HSiPM_F[BETA_SIZE + 1];
TH1D *HSiPM_GT[BETA_SIZE + 1];

TH2D *HSiPMRun[BETA_SIZE + 1];

TH2D *HSiliconRun_Channel_Unmatched[SIGNAL_MAX];
TH2D *HSiliconRun_Channel_Matched[SIGNAL_MAX];
TH2D *HSiliconRun[SIGNAL_MAX];
TH2D *HSiPMRun_Channel_Unmatched;

TH1D *Ref_Hist;
TTree *SiPM_Runs_Tree;
TGraph* graph;
TTreeReader *Reader_calib_sipm;
double chi2;
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

  double diffInHours = referenceDate.tm_hour + referenceDate.tm_min / 60.0 + referenceDate.tm_sec / 3600.0;

  return diffInHours;
}

bool Is_F(double E, int label)
{
  if (E < 1/SiliconCalib[99][label].first * peaks_window_F[label].first && E < 1/SiliconCalib[99][label].first * peaks_window_F[label].second)
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

void InitHistograms()
{
  //// BINING TIME /////////////////////////
  bin_time = std::numeric_limits<double>::max();
  for (const auto &p : file_Time)
  {
    double diff = std::abs(p.first - p.second);
    bin_time = std::min(bin_time, diff);
  }

  double t_min = file_Time[0].first;
  int t_max = (int)file_Time[file_Time.size() - 1].second + 1;
  n_bin_time = (int)(t_max - t_min) / bin_time;
  //////////////////////////////////////////

  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      HSiliconRun_Channel_Unmatched[i] = new TH2D(("HSiliconRun_Channel_Unmatched_" + detectorName[i]).c_str(), ("HSiliconRun_Channel_Unmatched_" + detectorName[i]).c_str(), t_max, t_min, t_max, eSiliN, eSiliMin, eSiliMax);
      HSiliconRun_Channel_Unmatched[i]->GetXaxis()->SetTitle("Time (h)");
      HSiliconRun_Channel_Unmatched[i]->GetYaxis()->SetTitle("Channel");
      HSiliconRun_Channel_Unmatched[i]->GetXaxis()->CenterTitle();
      HSiliconRun_Channel_Unmatched[i]->GetYaxis()->CenterTitle();

      HSiliconRun_Channel_Matched[i] = new TH2D(("HSiliconRun_Channel_Matched_" + detectorName[i]).c_str(), ("HSiliconRun_Channel_Matched_" + detectorName[i]).c_str(), t_max, t_min, t_max, eSiliN, eSiliMin, eSiliMax);
      HSiliconRun_Channel_Matched[i]->GetXaxis()->SetTitle("Time (h)");
      HSiliconRun_Channel_Matched[i]->GetYaxis()->SetTitle("Channel");
      HSiliconRun_Channel_Matched[i]->GetXaxis()->CenterTitle();
      HSiliconRun_Channel_Matched[i]->GetYaxis()->CenterTitle();

      HSiliconRun[i] = new TH2D(("HSiliconRun" + detectorName[i]).c_str(), ("HSiliconRun" + detectorName[i]).c_str(), t_max, t_min, t_max, 5000, 0, 10000);
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

      HSilicon_Channel_Matched_CURRENTRUN[i] = new TH1D(("HSilicon_Channel_Matched_CURRENTRUN_" + detectorName[i]).c_str(), ("HSilicon_Channel_Matched_CURRENTRUN_" + detectorName[i]).c_str(), eSiliN, eSiliMin, eSiliMax);
      HSilicon_Channel_Matched_CURRENTRUN[i]->GetXaxis()->SetTitle("Channel");
      HSilicon_Channel_Matched_CURRENTRUN[i]->GetYaxis()->SetTitle("Counts");
      HSilicon_Channel_Matched_CURRENTRUN[i]->GetXaxis()->CenterTitle();
      HSilicon_Channel_Matched_CURRENTRUN[i]->GetYaxis()->CenterTitle();

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
    HSiPMRun[GetDetectorChannel(i)] = new TH2D(("HSiPMRun_M" + to_string(i)).c_str(), ("HSiPMRun_M" + to_string(i)).c_str(), t_max, t_min, t_max, eLowN / 10, eLowMin, eLowMax * 8);
    HSiPMRun[GetDetectorChannel(i)]->GetXaxis()->SetTitle("Time (h)");
    HSiPMRun[GetDetectorChannel(i)]->GetYaxis()->SetTitle("Channel");
    HSiPMRun[GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
    HSiPMRun[GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

    HSiPM[GetDetectorChannel(i)] = new TH1D(("HSiPM_M" + to_string(i)).c_str(), ("HSiPM_M" + to_string(i)).c_str(), eLowN / 10, eLowMin, eLowMax * 5);
    HSiPM[GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
    HSiPM[GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
    HSiPM[GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
    HSiPM[GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

    HSiPM_F[GetDetectorChannel(i)] = new TH1D(("HSiPM_F_M" + to_string(i)).c_str(), ("HSiPM_F_M" + to_string(i)).c_str(), eLowN / 10, eLowMin, eLowMax * 5);
    HSiPM_F[GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
    HSiPM_F[GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
    HSiPM_F[GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
    HSiPM_F[GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

    HSiPM_GT[GetDetectorChannel(i)] = new TH1D(("HSiPM_GT_M" + to_string(i)).c_str(), ("HSiPM_GT_M" + to_string(i)).c_str(), eLowN / 10, eLowMin, eLowMax * 5);
    HSiPM_GT[GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
    HSiPM_GT[GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
    HSiPM_GT[GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
    HSiPM_GT[GetDetectorChannel(i)]->GetYaxis()->CenterTitle();
  }

  HSiPMRun_Channel_Unmatched = new TH2D("HSiPMRun_Channel_Unmatched", "HSiPMRun_Channel_Unmatched", t_max, t_min, t_max, eLowN / 10, eLowMin, eLowMax * 8);
  HSiPMRun_Channel_Unmatched->GetXaxis()->SetTitle("Time (h)");
  HSiPMRun_Channel_Unmatched->GetYaxis()->SetTitle("Channel");
  HSiPMRun_Channel_Unmatched->GetXaxis()->CenterTitle();
  HSiPMRun_Channel_Unmatched->GetYaxis()->CenterTitle();
}

void WriteHistograms()
{
  Merged_File->cd();
  TDirectory *dir_HSilicon_Channel_Unmatched = Merged_File->mkdir("Silicon_Unmatched");
  TDirectory *dir_HSilicon_Channel_Unmatched__Run = dir_HSilicon_Channel_Unmatched->mkdir("Run_Time");
  TDirectory *dir_HSilicon_Channel_Matched = Merged_File->mkdir("Silicon_Matched");
  TDirectory *dir_HSilicon_Channel_Matched__Run = dir_HSilicon_Channel_Matched->mkdir("Run_Time");
  TDirectory *dir_HSilicon = Merged_File->mkdir("Silicon_Calibrated");
  TDirectory *dir_HSilicon__Fitting = dir_HSilicon->mkdir("Fitting");
  TDirectory *dir_HSilicon__Run = dir_HSilicon->mkdir("Run_Time");
  TDirectory *dir_HSilicon__Hist = dir_HSilicon->mkdir("Hist");
  TDirectory *dir_HSiPM = Merged_File->mkdir("SiPM");

  double first_time = GetOnlyHour(file_string_Time[0].first);

  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      /////////////////UNMATCHED
      dir_HSilicon_Channel_Unmatched->cd();
      HSilicon_Channel_Unmatched[i]->Write();

      dir_HSilicon_Channel_Unmatched__Run->cd();
      TCanvas *c2 = new TCanvas(("RunsTime_Unmatched" + detectorName[i]).c_str(), ("RunsTime_Unmatched" + detectorName[i]).c_str(), 800, 600);
      c2->cd();
      TAxis *xaxis_unmatched = HSiliconRun_Channel_Unmatched[i]->GetXaxis();
      for (int bin = 1; bin <= xaxis_unmatched->GetNbins(); ++bin)
      {
        int hourmodulo = xaxis_unmatched->GetBinCenter(bin) + first_time;
        while (hourmodulo > 24)
        {
          hourmodulo -= 24;
        }
        if (hourmodulo % 2 == 0)
          xaxis_unmatched->SetBinLabel(bin, to_string(hourmodulo).c_str());
      }
      HSiliconRun_Channel_Unmatched[i]->GetXaxis()->SetTitle("Time (h)");
      HSiliconRun_Channel_Unmatched[i]->Draw("COLZ");
      c2->Write();
      delete c2;

      ////////////////MATCHED
      dir_HSilicon_Channel_Matched->cd();
      HSilicon_Channel_Matched[i]->Write();

      dir_HSilicon_Channel_Matched__Run->cd();
      TCanvas *c1 = new TCanvas(("RunsTime_Matched" + detectorName[i]).c_str(), ("RunsTime_Matched" + detectorName[i]).c_str(), 800, 600);
      c1->cd();
      TAxis *xaxis_matched = HSiliconRun_Channel_Matched[i]->GetXaxis();
      for (int bin = 1; bin <= xaxis_matched->GetNbins(); ++bin)
      {
        int hourmodulo = xaxis_matched->GetBinCenter(bin) + first_time;
        while (hourmodulo > 24)
        {
          hourmodulo -= 24;
        }
        if (hourmodulo % 2 == 0)
          xaxis_matched->SetBinLabel(bin, to_string(hourmodulo).c_str());
      }
      HSiliconRun_Channel_Matched[i]->GetXaxis()->LabelsOption("h");
      HSiliconRun_Channel_Matched[i]->GetXaxis()->SetTitle("Time (h)");
      HSiliconRun_Channel_Matched[i]->Draw("COLZ");
      c1->Write();
      delete c1;

      ////////////////CALIBRATED
      dir_HSilicon->cd();
      // HSilicon[i]->Write();

      dir_HSilicon__Run->cd();
      TCanvas *c0 = new TCanvas(("RunsTime_Calibrated" + detectorName[i]).c_str(), ("RunsTime_Calibrated" + detectorName[i]).c_str(), 800, 600);
      c0->cd();
      TAxis *xaxis_calibrated = HSiliconRun[i]->GetXaxis();
      for (int bin = 1; bin <= xaxis_calibrated->GetNbins(); ++bin)
      {
        int hourmodulo = xaxis_calibrated->GetBinCenter(bin) + first_time;
        while (hourmodulo > 24)
        {
          hourmodulo -= 24;
        }
        if (hourmodulo % 2 == 0)
          xaxis_calibrated->SetBinLabel(bin, to_string(hourmodulo).c_str());
      }
      HSiliconRun[i]->GetXaxis()->LabelsOption("h");
      HSiliconRun[i]->GetXaxis()->SetTitle("Time (h)");
      HSiliconRun[i]->Draw("COLZ");
      c0->Write();
      delete c0;

      dir_HSilicon__Hist->cd();
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

      dir_HSilicon__Fitting->cd();
      TCanvas *cc = new TCanvas(("Fitting_" + detectorName[i]).c_str(), ("Fitting_" + detectorName[i]).c_str(), 800, 600);
      cc->cd();
      TGraphErrors *graph1 = new TGraphErrors();
      graph1->SetTitle("Calibration Coefficients");
      graph1->GetXaxis()->SetTitle("Run");
      graph1->GetYaxis()->SetTitle("Calibration Coefficient [Channel/keV]");
      graph1->GetXaxis()->CenterTitle();
      graph1->GetYaxis()->CenterTitle();

      int counter = 0;
      for (int run = 0; run < 99; run++)
      {
        if (SiliconCalib[run][i].first != 0.)
        {
          graph1->SetPoint(counter, run, SiliconCalib[run][i].first);
          graph1->SetPointError(counter, 0, SiliconCalib[run][i].second);
          counter++;
        }
      }
      graph1->Draw("*AP");

      TLine *line = new TLine(graph1->GetPointX(0), SiliconCalib[99][i].first, graph1->GetPointX(graph1->GetN() - 1), SiliconCalib[99][i].first);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      line->Draw();

      TGraphErrors *graph2 = new TGraphErrors();
      graph2->SetPoint(0, graph1->GetPointX(0), SiliconCalib[99][i].first - SiliconCalib[99][i].second);
      graph2->SetPoint(1, graph1->GetPointX(graph1->GetN() - 1), SiliconCalib[99][i].first - SiliconCalib[99][i].second);
      graph2->SetPoint(2, graph1->GetPointX(graph1->GetN() - 1), SiliconCalib[99][i].first + SiliconCalib[99][i].second);
      graph2->SetPoint(3, graph1->GetPointX(0), SiliconCalib[99][i].first + SiliconCalib[99][i].second);
      graph2->SetFillColor(kRed);
      graph2->SetFillStyle(3005);
      graph2->Draw("F");

      cc->Write();
      delete cc;
    }
  }

  /////SiPM
  dir_HSiPM->cd();
  for (int i = 0; i <= BETA_SIZE; i++)
  {
    HSiPM[GetDetectorChannel(i)]->Write();
    HSiPM_F[GetDetectorChannel(i)]->Write();
    HSiPM_GT[GetDetectorChannel(i)]->Write();
  }

  TCanvas *c3 = new TCanvas("HSiPMRun", "HSiPMRun", 800, 600);
  c3->cd();
  TAxis *xaxis_sipm = HSiPMRun_Channel_Unmatched->GetXaxis();
  for (int bin = 1; bin <= xaxis_sipm->GetNbins(); ++bin)
  {
    int hourmodulo = xaxis_sipm->GetBinCenter(bin) + first_time;
    while (hourmodulo > 24)
    {
      hourmodulo -= 24;
    }
    if (hourmodulo % 2 == 0)
      xaxis_sipm->SetBinLabel(bin, to_string(hourmodulo).c_str());

    // TH1D *proj = HSiPMRun_Channel_Unmatched->ProjectionY("proj", bin, bin);
    // double integral = proj->Integral();
    // if (integral != 0)
    // {
    //   for (int j = 1; j <= HSiPMRun_Channel_Unmatched->GetNbinsY(); j++)
    //   {
    //     HSiPMRun_Channel_Unmatched->SetBinContent(bin, j, HSiPMRun_Channel_Unmatched->GetBinContent(bin, j) / integral);
    //   }
    // }
  }

  HSiPMRun_Channel_Unmatched->GetXaxis()->LabelsOption("h");
  HSiPMRun_Channel_Unmatched->GetXaxis()->SetTitle("Time (h)");
  HSiPMRun_Channel_Unmatched->Draw("COLZ");
  c3->Write();
  delete c3;
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
  TGraphErrors *MatchingGraph = (TGraphErrors *)Matched_File->Get("Silicon_RunMatching/Silicon_RunMatching");
  double label;
  double value;
  for (int counter = 1; counter < MatchingGraph->GetN(); counter++)
  {
    MatchingGraph->GetPoint(counter, label, value);
    detectorMatching[static_cast<int>(label)] = value;
  }
}

void MakeSiliconCalibration()
{
  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      TGraphErrors *graph = new TGraphErrors();
      HSilicon_Channel_Matched[i]->GetXaxis()->SetRangeUser(peaks_window_F[i].first, peaks_window_F[i].second);
      graph->SetPoint(0, SiliconCalib_Fermi[i], HSilicon_Channel_Matched[i]->GetMean());
      graph->SetPointError(0, 0, HSilicon_Channel_Matched[i]->GetMeanError());

      HSilicon_Channel_Matched[i]->GetXaxis()->SetRangeUser(peaks_window_GT[i].first, peaks_window_GT[i].second);
      graph->SetPoint(0, SiliconCalib_GT[i], HSilicon_Channel_Matched[i]->GetMean());
      graph->SetPointError(0, 0, HSilicon_Channel_Matched[i]->GetMeanError());

      linear->SetParameters(0, 12);
      graph->Fit("linear", "Q");

      SiliconCalib[99][i] = make_pair(linear->GetParameter(0), linear->GetParError(0));
    }
  }
}
void MakeSiliconCalibration(int run)
{
  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      TGraphErrors *graph = new TGraphErrors();
      HSilicon_Channel_Matched_CURRENTRUN[i]->GetXaxis()->SetRangeUser(peaks_window_F[i].first, peaks_window_F[i].second);
      graph->SetPoint(0, SiliconCalib_Fermi[i], HSilicon_Channel_Matched_CURRENTRUN[i]->GetMean());
      graph->SetPointError(0, 0, HSilicon_Channel_Matched_CURRENTRUN[i]->GetMeanError());

      HSilicon_Channel_Matched_CURRENTRUN[i]->GetXaxis()->SetRangeUser(peaks_window_GT[i].first, peaks_window_GT[i].second);
      graph->SetPoint(0, SiliconCalib_GT[i], HSilicon_Channel_Matched_CURRENTRUN[i]->GetMean());
      graph->SetPointError(0, 0, HSilicon_Channel_Matched_CURRENTRUN[i]->GetMeanError());

      linear->SetParameters(0, 12);
      graph->Fit("linear", "Q");

      SiliconCalib[run][i] = make_pair(linear->GetParameter(0), linear->GetParError(0));

      HSilicon_Channel_Matched_CURRENTRUN[i]->Reset();
    }
  }
}

inline double Chi2TreeHist(const Double_t *par)
{
    chi2 = 0;
    Reader_calib_sipm->Restart();
    TTreeReaderArray<double> Tree_SiPM(*Reader_calib_sipm, "PlasticScintillator_Deposit_Energy");
    TH1D *TreeHist = (TH1D *)Ref_Hist->Clone(("SiPM_" + to_string(par[0])).c_str());
    TreeHist->Reset();
    

    while (Reader_calib_sipm->Next())
    {
      normal_distribution<> resolution(0, par[1]*sqrt(par[0]*Tree_SiPM[0]));
      TreeHist->Fill(par[0] * Tree_SiPM[0] + resolution(gen));

    }
    cout<<"0 : " << par[0]<< "     1 : "<<par[1] <<endl;

    chi2 = Ref_Hist->Chi2Test(TreeHist, "CHI2/NDF");
    graph->SetPoint(counter, par[0], chi2);
    counter++;
    return chi2;
}

void MakeSiPMCalibration(int run)
{
  graph = new TGraph();
  graph->SetTitle("Calibration;Proportionality; #chi^{2}");
  counter = 0;
  int result_counter = 0;

  //////////////////// MINIMIZER ////////////////////
  ////// First
    // Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    // ROOT::Math::Functor functor(&Chi2TreeHist, 1);
    // minimizer->SetFunction(functor);
    // minimizer->SetLimitedVariable(0, "Proportionality", 0.25, 0.01, 0, 0.5);
    // minimizer->SetLimitedVariable(1, "res", 3, 0.5, 0, 30);
    // minimizer->SetPrecision(0.05);
    // minimizer->Minimize();
    // const double *bestPar = minimizer->X();

    // ////// Second
    // minimizer->SetLimitedVariable(0, "Proportionality", bestPar[0], 0.001, 0, 0.5);
    // minimizer->SetLimitedVariable(1, "res", bestPar[0], 0.5, 0, 30);
    // minimizer->SetPrecision(1e-5);
    // minimizer->Minimize();
    // bestPar = minimizer->X();

    // // ////// Third
    // // minimizer->SetLimitedVariable(0, "Proportionality", bestPar[0], .005, 0.8, 1.2);
    // // minimizer->SetPrecision(1e-10);
    // // minimizer->Minimize();
    // // bestPar = minimizer->X();
    // double error = minimizer->Errors()[0];
    // cout << "SiPM" << "\t CHI2 : " << chi2 << "\t\t Best Proportionality : " << bestPar[0] << " +/- " << error << std::endl;

  //////////////////// SAVE ////////////////////
  TCanvas *c = new TCanvas(("SiPM_"+to_string(run)).c_str(), ("SiPM_"+to_string(run)).c_str(), 800, 600);
  c->Divide(1, 2);
  c->cd(1);
  Ref_Hist->Draw("HIST");

  Reader_calib_sipm->Restart();
  TTreeReaderArray<double> Tree_SiPM(*Reader_calib_sipm, "PlasticScintillator_Deposit_Energy");

  TH1D *TreeHist = (TH1D *)Ref_Hist->Clone(("SiPM_" + to_string(0)).c_str());
  TreeHist->Reset();

  while (Reader_calib_sipm->Next())
  {
    normal_distribution<> resolution(0, 1*sqrt(0.25*Tree_SiPM[0]));
    TreeHist->Fill(0.25*Tree_SiPM[0] + resolution(gen));
  }

  TreeHist->SetLineColor(kRed);
  TreeHist->Draw("SAME");

  c->cd(2);
  graph->Draw("*AP");
  Merged_File->cd();
  c->Write();

  result_counter++;
  // Result_SiPM_Runs->SetPoint(result_counter, , bestPar[0]); /////SAVE FOR MERGED RUN
  SiPMCalib[run] = make_pair(0, 0);
  // Result_Si->SetPointError(result_counter, 0, error);
}

void InitSelection()
{
  string CalibFileName = "./Config_Files/Selection_32Ar_" + Catcher + ".txt"; //////////////USE MATERIAL
  ifstream file(CalibFileName);

  if (file.is_open())
  {
    GLogMessage("<SAM> Selection file found");

    string line;
    while (getline(file, line))
    {
      istringstream iss(line);
      string DetName;
      int run;
      vector<int> runs;
      iss >> DetName;

      for (size_t i = 0; i < detectorNum; ++i)
      {
        if (DetName == detectorName[i])
        {
          while (iss >> run)
          {
            RunDetectorSelection[run][i] = true;
          }
        }
      }
    }
  }
  else
  {
    GLogMessage("<SAM> No Run Selection file found");
    exit(0);
  }
}

void InitPeakWindow()
{
  string CalibFileName = "./Config_Files/Win_32Ar_" + Catcher + "_F.txt"; //////////////USE MATERIAL
  ifstream fileF(CalibFileName);

  for (auto &value : peaks_window_F)
  {
    value = make_pair(0, 0);
  }

  if (fileF.is_open())
  {
    GLogMessage("<SAM> Window file 1 found");

    string line;
    while (getline(fileF, line))
    {
      istringstream iss(line);
      int min;
      int max;
      string DetName;
      iss >> DetName >> min >> max;

      for (size_t i = 0; i < detectorNum; ++i)
      {
        if (DetName == detectorName[i])
        {
          peaks_window_F[i] = make_pair(min, max);
        }
      }
    }
  }
  else
  {
    GLogMessage("<SAM> No Window file found");
    exit(0);
  }

  CalibFileName = "./Config_Files/Win_32Ar_" + Catcher + "_GT.txt"; //////////////USE MATERIAL
  ifstream fileGT(CalibFileName);

  for (auto &value : peaks_window_GT)
  {
    value = make_pair(0, 0);
  }

  if (fileGT.is_open())
  {
    GLogMessage("<SAM> Window file 2 found");

    string line;
    while (getline(fileGT, line))
    {
      istringstream iss(line);
      int min;
      int max;
      string DetName;
      iss >> DetName >> min >> max;

      for (size_t i = 0; i < detectorNum; ++i)
      {
        if (DetName == detectorName[i])
        {
          peaks_window_GT[i] = make_pair(min, max);
        }
      }
    }
  }
  else
  {
    GLogMessage("<SAM> No Window file found");
    exit(0);
  }
}

int InitCalib()
{

  ///////////////////FOR FERMI
  int error = 0;
  string CalibDir = "./Calibration_Data/";
  string CalibFileName = "Sim_32Ar_AlMylar_F.txt"; //////////////USE MATERIAL
  ifstream file(CalibDir + CalibFileName);

  for (auto &value : SiliconCalib_Fermi)
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
            SiliconCalib_Fermi[i] = E;
          }
        }
      }
      else if (DetName.find("Down") != string::npos)
      {
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
          if (GetDetector(i) > 4 and GetDetectorChannel(i) == strip)
          {
            SiliconCalib_Fermi[i] = E;
          }
        }
      }
    }
  }
  else
  {
    GLogMessage("<SAM> No Calibration file found");
    error = 1;
    exit(0);
  }

  ////////////////////////////////////////
  CalibFileName = "Sim_32Ar_AlMylar_GT1.txt";
  ifstream fileGT(CalibDir + CalibFileName);

  for (auto &value : SiliconCalib_Fermi)
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
            SiliconCalib_GT[i] = E;
          }
        }
      }
      else if (DetName.find("Down") != string::npos)
      {
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
          if (GetDetector(i) > 4 and GetDetectorChannel(i) == strip)
          {
            SiliconCalib_GT[i] = E;
          }
        }
      }
    }
  }
  else
  {
    GLogMessage("<SAM> No Calibration file found");
    error = 1;
    exit(0);
  }

  // FOR CHECK
  //  int index = 0;
  //  for (auto &name : detectorName)
  //  {
  //    cout << name << "       " << SiliconCalib_Fermi[index] << endl;
  //    index++;
  //  }

  ///////FOR SIPM
  TFile *Calibration_File = new TFile("../../../../../../mnt/hgfs/shared-2/CALIB/32Ar__a1_b0.root");
  TTree* tree = (TTree *)Calibration_File->Get("Tree");
  Reader_calib_sipm = new TTreeReader(tree);

  file.close();
  return error;
}

void SetTime(TFile *file)
{
  file_string_Time.push_back(GetTime(file));

  double start_time = Convert_DatetoTime(file_string_Time[file_string_Time.size() - 1].first, file_string_Time[0].first);
  double stop_time = Convert_DatetoTime(file_string_Time[file_string_Time.size() - 1].second, file_string_Time[0].first);
  file_Time.push_back(make_pair(start_time, stop_time));
}

#endif