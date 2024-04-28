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
#include <array>

#include "TFile.h"
#include <TStyle.h>
#include "TCanvas.h"
#include <TKey.h>
#include <TMath.h>
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
#include "TPaveText.h"
#include "TLegend.h"

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

TH1D *Ref_Hist_F[BETA_SIZE + 1][BETA_SIZE + 1];
TH1D *Ref_Hist_Multi_F[BETA_SIZE + 1];
TH1D *TreeHist_F[BETA_SIZE + 1][BETA_SIZE + 1];
TH1D *TreeHist_Multi_F[BETA_SIZE + 1];
TH2D *TreeHist_Conv_F[BETA_SIZE + 1];
TH2D *TreeHist_Conv_F_av;
TH1D *Ref_Hist_F_SUM;
string current_option;
Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
pair<string, double> coefficents[7] = {make_pair("Calibration_OffSet", 0), make_pair("Calibration", 0), make_pair("Resolution_OffSet", 0), make_pair("Resolution_SQRT", 0), make_pair("Resolution_2", 0), make_pair("Threshold", 0), make_pair("Threshold_STD", 0)};
TTree *SiPM_Runs_Tree;
TGraph *graph_sqrt;
double bestCHI2 = 10000;
TTreeReader *Reader_calib_sipm_F;
TTreeReader *Reader_calib_sipm_GT;
double chi2;
vector<double> chi2_vec;
double chi2_multiplicity[BETA_SIZE + 1];
TTree *Ref_Tree_F[BETA_SIZE + 1][BETA_SIZE + 1];
TTree *Ref_Tree_Multi_F[BETA_SIZE + 1];
TF1 *Threshold;
TH1D* h_probability[BETA_SIZE + 1];
TH1D* h_ratio;
vector<double> Par;
vector<array<double, 7>> bestPar_sipm;
int current_sipm;
TH1D *checker_si[SIGNAL_MAX];
TF1 *linear = new TF1("linear", "[0]*x");
TH1D *histF = new TH1D("histF", "histF", 800, 0, 8000);
TH1D *histGT = new TH1D("histGT", "histGT", 1000, 0, 10000);
TGraph *graph_res[BETA_SIZE + 1];
TGraph *graph_cal[BETA_SIZE + 1];
TGraph *graph_calib;
TGraph *graph_resol;
TH1D *TreeHist_Conv_F_1D[8];
TF1 *dark;
TFile *dark_file;

int Verbose = 0;

bool Calibration = false;

int counter;

double bin_time;
int n_bin_time;
double start_time;
double stop_time;

void generate_combinations(vector<pair<int, double>> &elements, int combination_size, int start, vector<pair<int, double>> &current_combination, double Channel, double Label, TTree *Tree)
{
  if (current_combination.size() == combination_size)
  {
    for (auto i : current_combination)
    {
      Label = i.first;
      Channel = i.second;
      Tree->Fill();
    }
  }
  else
  {
    for (int i = start; i < elements.size(); ++i)
    {
      current_combination.push_back(elements[i]);
      generate_combinations(elements, combination_size, i + 1, current_combination, Channel, Label, Tree);
      current_combination.pop_back();
    }
  }
}

void generate_combinations_sim(vector<double> &elements, int combination_size, int start, vector<double> &current_combination, TH1D *Hist)
{
  if (current_combination.size() == combination_size)
  {
    double sum = 0;
    for (auto i : current_combination)
    {
      sum += i;
    }
    Hist->Fill(sum / combination_size);
  }
  else
  {
    for (int i = start; i < elements.size(); ++i)
    {
      current_combination.push_back(elements[i]);
      generate_combinations_sim(elements, combination_size, i + 1, current_combination, Hist);
      current_combination.pop_back();
    }
  }
}
///
////////////////////////////////////

double Convert_DatetoTime(string datestring, string refstring)
{
  tm date = {}, referenceDate = {};
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
  tm referenceDate = {};
  istringstream referenceDateStream(refstring);
  referenceDateStream >> get_time(&referenceDate, "%d-%m-%Y %H:%M:%S");

  time_t referenceTimeInSeconds = mktime(&referenceDate);

  double diffInHours = referenceDate.tm_hour + referenceDate.tm_min / 60.0 + referenceDate.tm_sec / 3600.0;

  return diffInHours;
}

bool Is_F(double E, int label)
{
  if (E > 1 / SiliconCalib[99][label].first * peaks_window_F[label].first && E < 1 / SiliconCalib[99][label].first * peaks_window_F[label].second)
  {
    return true;
  }
  return false;
}

bool Is_GT(double E, int label)
{
  if (E > 1 / SiliconCalib[99][label].first * peaks_window_GT[label].first && E < 1 / SiliconCalib[99][label].first * peaks_window_GT[label].second)
  {
    return true;
  }
  return false;
}

void InitHistograms()
{
  //// BINING TIME /////////////////////////
  bin_time = numeric_limits<double>::max();
  for (const auto &p : file_Time)
  {
    double diff = abs(p.first - p.second);
    bin_time = min(bin_time, diff);
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

      checker_si[i] = new TH1D(("checker_si_" + detectorName[i]).c_str(), ("checker_si_" + detectorName[i]).c_str(), 10000, 0, 10000);
      checker_si[i]->GetXaxis()->SetTitle("Energy [keV]");
      checker_si[i]->GetYaxis()->SetTitle("Counts");
      checker_si[i]->GetXaxis()->CenterTitle();
      checker_si[i]->GetYaxis()->CenterTitle();
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
      checker_si[i]->Write();
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
      HSilicon_coinc[i]->Draw("SAME HIST");
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
    RunNumbers = {25, 26, 31, 32, 33, 34, 36, 37, 38};
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
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// void FillHistogramsChi2_old(const Double_t *bestPar)
// {
//   for (int index = 0; index < minimizer->NDim(); index++)
//   {
//     for (auto &p : coefficents)
//     {
//       if (p.first == minimizer->VariableName(index))
//       {
//         p.second = bestPar[index];
//         // cout << p.first << " : " << p.second << endl;
//         break;
//       }
//     }
//   }

//   Reader_calib_sipm_F->Restart();
//   TTreeReaderValue<double> Tree_SiPM(*Reader_calib_sipm_F, "PlasticScintillator_Deposit_Energy");

//   ////////////// conversion exp en kev///////////

//   Ref_Hist_F[1][current_sipm]->Reset();
//   TTreeReader *Reader = new TTreeReader(Ref_Tree_F[1][current_sipm]);
//   TTreeReaderValue<double> value(*Reader, "Channel");
//   while (Reader->Next())
//   {
//     Ref_Hist_F[1][current_sipm]->Fill(coefficents[0].second + coefficents[1].second * (*value));
//   }

//   /////////////////////////////////////////////////

//   vector<double> vec_e;
//   double e_sum;
//   double sigma_resolution = 0;
//   double sigma_resolution_sum = 0;
//   bool negative_res = false;

//   for (int i = 0; i < 10; i++)
//   {
//     Reader_calib_sipm_F->Restart();
//     while (Reader_calib_sipm_F->Next())
//     {
//       double energy = (*Tree_SiPM);
//       int real_multiplicity = 0;
//       vec_e.clear();
//       e_sum = 0;
//       sigma_resolution = sqrt(pow(coefficents[2].second, 2) + pow(coefficents[3].second * sqrt(energy), 2) + pow(coefficents[4].second * pow(energy, 2), 2));

//       normal_distribution<> resolution(0, sigma_resolution);
//       energy += resolution(gen);
//       normal_distribution<> threshold(coefficents[5].second, coefficents[6].second);

//       if (energy > threshold(gen))
//       {
//         TreeHist_F[1][current_sipm]->Fill(energy);
//       }
//     }
//   }
// }

void FillHistogramsChi2(const Double_t *bestPar)
{
  for (int index = 0; index < minimizer->NDim(); index++)
  {
    for (auto &p : coefficents)
    {
      if (p.first == minimizer->VariableName(index))
      {
        p.second = bestPar[index];
        // cout << p.first << " : " << p.second << endl;
        break;
      }
    }
  }

  Reader_calib_sipm_F->Restart();
  TTreeReaderValue<double> Tree_SiPM(*Reader_calib_sipm_F, "PlasticScintillator_Deposit_Energy");

  ////////////// conversion exp en kev///////////

  Ref_Hist_F[1][current_sipm]->Reset();
  TTreeReader *Reader = new TTreeReader(Ref_Tree_F[1][current_sipm]);
  TTreeReaderValue<double> value(*Reader, "Channel");
  while (Reader->Next())
  {
    Ref_Hist_F[1][current_sipm]->Fill(coefficents[0].second + coefficents[1].second * (*value));
  }

  /////////////////////////////////////////////////

  TF1 *gauss;
  Merged_File->cd();
  TreeHist_F[1][current_sipm]->Reset();

  for (int i = 1; i <= histF->GetNbinsX(); i++)
  {
    double energy = histF->GetBinCenter(i);
    double sigma_resolution = sqrt(pow(coefficents[2].second, 2) + pow(coefficents[3].second * sqrt(energy), 2) + pow(coefficents[4].second * pow(energy, 2), 2));

    gauss = new TF1("gauss", "gaus", 0, 8000);
    gauss->SetNpx(800);
    gauss->SetParameters(histF->GetBinContent(i) / (sigma_resolution * sqrt(2 * M_PI)), energy, sigma_resolution); /// weighted gaussian

    TreeHist_F[1][current_sipm]->Add((TH1D *)gauss->GetHistogram());
    delete gauss;
  }

  Threshold = new TF1("Threshold", "gaus", 0, 8000);
  Threshold->SetParameters(1, coefficents[5].second, coefficents[6].second);
  for (int i = 1; i <= histF->GetNbinsX(); i++)
  {
    TreeHist_F[1][current_sipm]->SetBinContent(i, TreeHist_F[1][current_sipm]->GetBinContent(i) * Threshold->Integral(0, histF->GetBinCenter(i)));
  }
}

inline double Chi2TreeHist(const Double_t *bestPar)
{
  double chi2 = 0;
  double pvalue = 0;
  chi2_vec.clear();

  TreeHist_F[1][current_sipm] = (TH1D *)Ref_Hist_F[1][current_sipm]->Clone(("SiPM_FF" + to_string(current_sipm)).c_str());
  TreeHist_F[1][current_sipm]->Reset();

  FillHistogramsChi2(bestPar);

  int min = 250;
  int max = 2000;

  Ref_Hist_F[1][current_sipm]->GetXaxis()->SetRangeUser(0, 8000);
  TreeHist_F[1][current_sipm]->GetXaxis()->SetRangeUser(0, 8000);

  pvalue = Ref_Hist_F[1][current_sipm]->Chi2Test(TreeHist_F[1][current_sipm], "");
  chi2 = Ref_Hist_F[1][current_sipm]->Chi2Test(TreeHist_F[1][current_sipm], "CHI2/NDF");

  cout << setprecision(5) << "p-value = " << pvalue << "  Chi2 = " << setprecision(5) << chi2 << "    {";
  for (size_t i = 0; i < 7; ++i)
  {
    cout << bestPar[i] << ", ";
  }
  cout << "});";
  if (bestCHI2 > abs(chi2 - 1))
  {
    bestCHI2 = abs(chi2 - 1);
    cout << "    BEST" << endl;
  }
  else
  {
    cout << endl;
  }
  return abs(chi2 - 1);
}

void MakeSiPMCalibration(int run)
{
  counter = 0;
  int result_counter = 0;

  // double bestPar[6] = {0.21575, 1.96249, 1.62856e-06, 0., 13.9911, 2.74464};          /// BEST    Positive offset
  // double bestPar[7] = {0.5, 0.203913, 0, 1.85754, 6.62335e-05, 13.3651, 2.12369};       //BEST    Negative offset propop, res, res2, offset, th, th_sigma

  // double bestPar[7] = {-3, 0.200901, 0., 1.9503, 5e-5, 14., 3.26826, }; //// BEST FOR MULTIPLICITY*
  //  //// BEST FOR MULTIPLICITY

  TCanvas *cF = new TCanvas(("SiPM_F" + to_string(run)).c_str(), ("SiPM_F" + to_string(run)).c_str(), 1920, 1080);
  TLegend *legend = new TLegend(0.1, 0.5, 0.25, 0.9);
  cF->Divide(3, 3);

  vector<double> chi2_sipm(9, 0);
  int countercd = 0;
  double sigma = 0;
  double value = 0;

  // bestPar_sipm.push_back({ 0, 0,0,0,0,0,0});
  // bestPar_sipm.push_back({7.91006, 5.01133, 0, 4.6819, 1.99548e-05, 83.8734, 12.8974, });
  // bestPar_sipm.push_back({6.06931, 4.96007, 0, 3.82433, 2.24847e-05, 81.10    UWminimizer->SetPrecision(0.000000001);
  // bestPar_sipm.push_back({2.58732, 5.0408, 0, 2.90152, 2.9204e-05, 87.252, 16.7987, });
  // bestPar_sipm.push_back({ 0, 0,0,0,0,0,0});
  // bestPar_sipm.push_back({2.7675, 4.99999, 0, 2.9588, 2.78769e-05, 83.4966, 17.721, });
  // bestPar_sipm.push_back({14.9584, 4.83285, 0, 2.11295, 3.7311e-05, 74.5177, 12.5693, });
  // bestPar_sipm.push_back({8.32008, 5.00363, 0, 4.99757, 3.31788e-05, 83.5593, 15.6704, });
  // bestPar_sipm.push_back({9.31388, 4.9, 0, 2.0, 2.21952e-05, 81.7, 18, });

  bestPar_sipm.push_back({
      0,
      0,
      0,
      0,
      0,
      0,
      0,
  });
  bestPar_sipm.push_back({
      38.5796,
      4.98336,
      0,
      3.80067,
      2.12837e-05,
      106.496,
      11.0801,
  });
  bestPar_sipm.push_back({
      38.9349,
      4.85626,
      0,
      3.81807,
      1.93505e-05,
      105.7508,
      13.657,
  });
  bestPar_sipm.push_back({
      38.9742,
      4.97679,
      0,
      3.8,
      2.46495e-05,
      112.127,
      13.2677,
  });
  bestPar_sipm.push_back({
      0,
      0,
      0,
      0,
      0,
      0,
      0,
  });
  bestPar_sipm.push_back({
      37.9151,
      4.96607,
      0,
      3.8,
      2.82445e-05,
      111.761,
      15.992,
  });
  bestPar_sipm.push_back({
      40.066,
      4.85737,
      0,
      3.50675,
      3.34238e-05,
      95.9868,
      9.83632,
  });
  bestPar_sipm.push_back({
      37.6385,
      4.93281,
      0,
      3.7357,
      3.39992e-05,
      101.591,
      13.2062,
  });
  bestPar_sipm.push_back({
      40.3775,
      4.87366,
      0,
      3.50013,
      2.36671e-05,
      104.521,
      15.734,
  });

  //   bestPar_sipm.push_back({0, 0, 0, 0, 0, 0, 0, });
  // bestPar_sipm.push_back({39, 4.9, 0, 3.7, 2.3e-05, 106.496, 12.0801, });
  // bestPar_sipm.push_back({39, 4.9, 0, 3.7, 2.3e-05, 106.496, 12.0801, });
  // bestPar_sipm.push_back({39, 4.9, 0, 3.7, 2.3e-05, 106.496, 12.0801, });
  // bestPar_sipm.push_back({39, 4.9, 0, 3.7, 2.3e-05, 106.496, 12.0801, });
  // bestPar_sipm.push_back({39, 4.9, 0, 3.7, 2.3e-05, 106.496, 12.0801, });
  // bestPar_sipm.push_back({39, 4.9, 0, 3.7, 2.3e-05, 106.496, 12.0801, });
  // bestPar_sipm.push_back({39, 4.9, 0, 3.7, 2.3e-05, 106.496, 12.0801, });
  // bestPar_sipm.push_back({39, 4.9, 0, 3.7, 2.3e-05, 106.496, 12.0801, });

  for (int sipm = 0; sipm <= 9; sipm++)
  {
    if (sipm == 4 || sipm == 0 || sipm == 9)
    {
      bestPar_sipm[sipm] = {0, 0, 0, 0, 0, 0, 0};
      continue;
    }
    countercd++;
    cout << "SiPM : " << sipm << endl;
    current_sipm = sipm;

    ROOT::Math::Functor functor(&Chi2TreeHist, 7);
    minimizer->SetFunction(functor);
    minimizer->SetLimitedVariable(0, "Calibration_OffSet", bestPar_sipm[sipm][0], 5, 30, 45);
    minimizer->SetLimitedVariable(1, "Calibration", bestPar_sipm[sipm][1], 0.1, 4.85, 5.2);
    minimizer->SetFixedVariable(2, "Resolution_OffSet", bestPar_sipm[sipm][2]);
    minimizer->SetLimitedVariable(3, "Resolution_SQRT", bestPar_sipm[sipm][3], 1, 3.5, 3.9);
    minimizer->SetLimitedVariable(4, "Resolution_2", bestPar_sipm[sipm][4], 1e-6, 1e-5, 4e-5);
    minimizer->SetLimitedVariable(5, "Threshold", bestPar_sipm[sipm][5], 5, 80, 150);
    minimizer->SetLimitedVariable(6, "Threshold_STD", bestPar_sipm[sipm][6], 10, 5, 16);
    minimizer->SetPrecision(0.001);
    // minimizer->SetTolerance(0.000000001);
    minimizer->SetMaxFunctionCalls(10000000);
    minimizer->SetMaxIterations(10000000);
    bestCHI2 = 1e6;

    // minimizer->Minimize();
    // const double *bestPar = minimizer->X();
    // minimizer->PrintResults();

    // bestPar_sipm[sipm] = { bestPar[0], bestPar[1], bestPar[2], bestPar[3], bestPar[4], bestPar[5], bestPar[6] };

    const double *bestPar = bestPar_sipm[sipm].data();

    chi2_sipm[sipm] = 1 + Chi2TreeHist(bestPar);

    cF->cd(countercd);
    TreeHist_F[1][sipm]->GetXaxis()->SetRangeUser(0, 8000);
    Ref_Hist_F[1][sipm]->GetXaxis()->SetRangeUser(0, 8000);
    TreeHist_F[1][sipm]->Scale(Ref_Hist_F[1][sipm]->Integral() / TreeHist_F[1][sipm]->Integral());

    Ref_Hist_F[1][sipm]->GetXaxis()->SetRangeUser(0, 2000);
    TreeHist_F[1][sipm]->GetXaxis()->SetRangeUser(0, 2000);
    Ref_Hist_F[1][sipm]->Draw("HIST");
    TreeHist_F[1][sipm]->SetLineColor(kRed);
    TreeHist_F[1][sipm]->SetLineWidth(1);
    TreeHist_F[1][sipm]->Draw("SAME HIST");

    TPaveText *pt = new TPaveText(0.7, 0.65, 1, 0.7, "brNDC");
    pt->SetFillColor(0);
    pt->AddText(("#chi^{2} = " + to_string(chi2_sipm[sipm])).c_str());
    pt->Draw("SAME");

    cF->cd(8);
    graph_cal[sipm] = new TGraph();
    graph_cal[sipm]->SetTitle(("Calibration;Energy[keV]; Channel"));
    int counter_cal = 0;

    for (double e = 1; e <= 6000; e += 10)
    {
      counter_cal++;
      value = (e - bestPar_sipm[sipm][0]) / bestPar_sipm[sipm][1];
      graph_cal[sipm]->SetPoint(counter_cal, e, value);
    }
    graph_cal[sipm]->SetLineColor(sipm);
    if (sipm == 1)
    {
      graph_cal[sipm]->Draw("AL");
    }
    else
    {
      graph_cal[sipm]->Draw("SAME");
    }

    cF->cd(9);
    graph_res[sipm] = new TGraph();
    graph_res[sipm]->SetTitle(("Resolution;Energy[keV]; #sigma [keV]"));
    legend->AddEntry(graph_res[sipm], ("SiPM " + to_string(sipm)).c_str(), "l");
    int counter = 0;

    for (double e = 1; e <= 6000; e += 10)
    {
      counter++;
      sigma = sqrt(pow(bestPar_sipm[sipm][2], 2) + pow(bestPar_sipm[sipm][3] * sqrt(e), 2) + pow(bestPar_sipm[sipm][4] * pow(e, 2), 2));
      graph_res[sipm]->SetPoint(counter, e, sigma);
    }
    graph_res[sipm]->SetLineColor(sipm);
    if (sipm == 1)
    {
      graph_res[sipm]->Draw("AL");
    }
    else
    {
      graph_res[sipm]->Draw("SAME");
    }
  }
  cF->cd(9);
  legend->Draw("SAME");
  cF->Write();
  Merged_File->cd();

  for (int i = 0; i < 9; i++)
  {
    cout << "bestPar_sipm.push_back({";
    for (int j = 0; j < 7; j++)
    {
      cout << bestPar_sipm[i][j] << ", ";
    }
    cout << "});" << endl;
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// void MakeSiPM_MultiplicityPlots()
// {
//   double chi2 = 0;
//   double chi2_std = 0;
//   chi2_vec.clear();

//   for (int multiplicity = 0; multiplicity <= 7; multiplicity++)
//   {
//     TreeHist_Multi_F[multiplicity] = (TH1D *)Ref_Hist_Multi_F[multiplicity]->Clone(("SiPM_MF" + to_string(multiplicity)).c_str());
//     TreeHist_Multi_F[multiplicity]->Reset();
//   }

//   Reader_calib_sipm_F->Restart();
//   TTreeReaderValue<double> Tree_SiPM(*Reader_calib_sipm_F, "PlasticScintillator_Deposit_Energy");

//   ////////////// conversion exp en kev///////////
//   for (int i = 0; i < BETA_SIZE + 1; i++)
//   {
//     Ref_Hist_Multi_F[i]->Reset();
//     TTreeReader *Reader = new TTreeReader(Ref_Tree_Multi_F[i]);
//     TTreeReaderValue<double> value(*Reader, "Channel");
//     TTreeReaderValue<int> label(*Reader, "Label");
//     int counter_sum = 0;
//     double value_sum = 0;
//     while (Reader->Next())
//     {
//       // if (i == 0) /////mean of sipm
//       // {
//       //   counter_sum++;
//       //   value_sum += (bestPar_sipm[*label][0] + bestPar_sipm[*label][1] * (*value));

//       //   if (counter_sum == 7)
//       //   {
//       //     Ref_Hist_Multi_F[i]->Fill(value_sum / 7);
//       //     value_sum = 0;
//       //     counter_sum = 0;
//       //   }
//       // }

//         Ref_Hist_Multi_F[i]->Fill(bestPar_sipm[*label][0] + bestPar_sipm[*label][1] * (*value));

//     }
//   }
//   /////////////////////////////////////////////////

//   vector<double> vec_e;
//   double energy_sum = 0;
//   double sigma_resolution = 0;
//   double sigma_resolution_SUM = 0;
//   bool negative_res = false;

//   for (int i = 0; i < 10; i++)
//   {
//     Reader_calib_sipm_F->Restart();
//     while (Reader_calib_sipm_F->Next())
//     {
//       int real_multiplicity = 0;
//       vec_e.clear();
//       double energy_sum_final = 0;

//       for (int sipm = 1; sipm <= 8; sipm++)
//       {
//         double energy = (*Tree_SiPM);
//         energy_sum = (*Tree_SiPM);
//         if (sipm == 4)
//         {
//           continue;
//         }
//         sigma_resolution = sqrt(pow(bestPar_sipm[sipm][2], 2) + pow(bestPar_sipm[sipm][3] * sqrt(energy), 2) + pow(bestPar_sipm[sipm][4] * pow(energy, 2), 2));
//         normal_distribution<> resolution(0, sigma_resolution);

//         sigma_resolution_SUM = sqrt(pow(bestPar_sipm[sipm][2], 2) + pow(bestPar_sipm[sipm][3] * sqrt(energy), 2));
//         normal_distribution<> resolution_SUM(0, sigma_resolution_SUM);

//         energy += resolution(gen);
//         energy_sum += resolution_SUM(gen);

//         normal_distribution<> threshold(bestPar_sipm[sipm][5], bestPar_sipm[sipm][6]);
//         if (energy > threshold(gen))
//         {
//           real_multiplicity++;
//           vec_e.push_back(energy);
//           energy_sum_final += energy_sum;
//         }
//       }

//       for (int multi = 1; multi <= real_multiplicity; multi++)
//       {
//         for (int sipm = 1; sipm <= real_multiplicity; sipm++)
//         {
//           TreeHist_Multi_F[multi]->Fill(vec_e[sipm - 1]);
//         }
//       }

//       if (real_multiplicity == 7)
//       {
//         TreeHist_Multi_F[0]->Fill(energy_sum_final / 7);
//       }
//     }
//   }

//   int min = 250;
//   int max = 2000;

//   for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
//   {
//     Ref_Hist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 8000);
//     TreeHist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 8000);
//     chi2_vec.push_back(Ref_Hist_Multi_F[multiplicity]->Chi2Test(TreeHist_Multi_F[multiplicity], "CHI2/NDF"));
//     cout << Ref_Hist_Multi_F[multiplicity]->Chi2Test(TreeHist_Multi_F[multiplicity], "CHI2/NDF") << endl;

//   }

//   double sum = accumulate(chi2_vec.begin(), chi2_vec.end(), 0.0);
//   chi2 = sum / chi2_vec.size();

//   double sq_sum = inner_product(chi2_vec.begin(), chi2_vec.end(), chi2_vec.begin(), 0.0);
//   chi2_std = sqrt(sq_sum / chi2_vec.size() - chi2 * chi2);

//   cout << chi2 << " +/- " << chi2_std << endl;;

//   /////////////////PLOTTING
//   TCanvas *cF = new TCanvas("SiPM_F", "SiPM_F", 1920, 1080);
//   cF->Divide(3, 3);

//   for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
//   {
//     cF->cd(multiplicity);
//     TreeHist_Multi_F[7]->GetXaxis()->SetRangeUser(0, 8000);
//     Ref_Hist_Multi_F[7]->GetXaxis()->SetRangeUser(0, 8000);
//     TreeHist_Multi_F[multiplicity]->Scale(Ref_Hist_Multi_F[7]->Integral() / TreeHist_Multi_F[7]->Integral());

//     Ref_Hist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
//     TreeHist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
//     Ref_Hist_Multi_F[multiplicity]->Draw("HIST");
//     TreeHist_Multi_F[multiplicity]->SetLineColor(kRed);
//     TreeHist_Multi_F[multiplicity]->Draw("SAME HIST");

//     TPaveText *pt = new TPaveText(0.7, 0.65, 1, 0.7, "brNDC");
//     pt->SetFillColor(0);
//     pt->AddText(("#chi^{2} = " + to_string(chi2_vec[multiplicity - 1])).c_str());
//     pt->Draw("SAME");
//   }

//   Merged_File->cd();
//   cF->Write();
// }

////////////////////////////////////////////////////////////////////////////////
///////function JUST to plot all the multiplicity FOR ONE GIVEN SIPM////////////
////////////////////////////////////////////////////////////////////////////////
void MakeSiPM_SiPMPlots()
{
  double chi2 = 0;
  double chi2_std = 0;
  chi2_vec.clear();

  for (int multiplicity = 0; multiplicity <= 9; multiplicity++)
  {
    for (int sipm = 1; sipm <= 9; sipm++)
    {
      TreeHist_F[multiplicity][sipm] = (TH1D *)Ref_Hist_F[multiplicity][sipm]->Clone(("SiPM" + to_string(sipm) + "_F_M" + to_string(multiplicity)).c_str());
      TreeHist_F[multiplicity][sipm]->Reset();
    }
  }

  Reader_calib_sipm_F->Restart();
  TTreeReaderValue<double> Tree_SiPM(*Reader_calib_sipm_F, "PlasticScintillator_Deposit_Energy");

  ////////////// conversion exp en kev///////////
  for (int multiplicity = 0; multiplicity <= 9; multiplicity++)
  {
    for (int sipm = 1; sipm <= 9; sipm++)
    {
      Ref_Hist_F[multiplicity][sipm]->Reset();
      TTreeReader *Reader = new TTreeReader(Ref_Tree_F[multiplicity][sipm]);
      TTreeReaderValue<double> value(*Reader, "Channel");
      while (Reader->Next())
      {
        Ref_Hist_F[multiplicity][sipm]->Fill(bestPar_sipm[sipm][0] + bestPar_sipm[sipm][1] * (*value));
      }
    }
  }
  /////////////////////////////////////////////////

  vector<double> vec_e = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double sigma_resolution = 0;
  bool negative_res = false;

  for (int i = 0; i < 10; i++)
  {
    Reader_calib_sipm_F->Restart();
    while (Reader_calib_sipm_F->Next())
    {
      int real_multiplicity = 0;
      vec_e.clear();

      for (int sipm = 1; sipm <= 8; sipm++)
      {
        double energy = (*Tree_SiPM);
        if (sipm == 4)
        {
          continue;
        }
        sigma_resolution = sqrt(pow(bestPar_sipm[sipm][2], 2) + pow(bestPar_sipm[sipm][3] * sqrt(energy), 2) + pow(bestPar_sipm[sipm][4] * pow(energy, 2), 2));
        normal_distribution<> resolution(0, sigma_resolution);

        energy += resolution(gen);

        normal_distribution<> threshold(bestPar_sipm[sipm][5], bestPar_sipm[sipm][6]);
        if (energy > threshold(gen))
        {
          real_multiplicity++;
          vec_e[sipm] = energy;
        }
      }

      for (int multi = 1; multi <= real_multiplicity; multi++)
      {
        for (int sipm = 1; sipm <= 9; sipm++)
        {
          TreeHist_F[multi][sipm]->Fill(vec_e[sipm]);
        }
      }
    }
  }

  for (int sipm = 1; sipm <= 8; sipm++)
  {
    if (sipm == 4)
    {
      continue;
    }
    /////////////////PLOTTING
    TCanvas *cF = new TCanvas(("SiPM" + to_string(sipm) + "_F").c_str(), ("SiPM" + to_string(sipm) + "_F").c_str(), 1920, 1080);
    cF->Divide(3, 3);

    for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
    {
      cF->cd(multiplicity);

      TPaveText *pt = new TPaveText(0.7, 0.65, 1, 0.7, "brNDC");
      pt->SetFillColor(0);
      Ref_Hist_F[7][sipm]->GetXaxis()->SetRangeUser(0, 8000);
      TreeHist_F[7][sipm]->GetXaxis()->SetRangeUser(0, 8000);
      pt->AddText(("#chi^{2} = " + to_string(Ref_Hist_F[multiplicity][sipm]->Chi2Test(TreeHist_F[multiplicity][sipm], " CHI2/NDF"))).c_str());
      TreeHist_F[multiplicity][sipm]->Scale(Ref_Hist_F[7][sipm]->Integral() / TreeHist_F[7][sipm]->Integral());

      Ref_Hist_F[multiplicity][sipm]->GetXaxis()->SetRangeUser(0, 8000);
      TreeHist_F[multiplicity][sipm]->GetXaxis()->SetRangeUser(0, 8000);
      Ref_Hist_F[multiplicity][sipm]->Draw("HIST");
      TreeHist_F[multiplicity][sipm]->SetLineColor(kRed);
      TreeHist_F[multiplicity][sipm]->Draw("SAME HIST");

      pt->Draw("SAME");
    }

    Merged_File->cd();
    cF->Write();
  }
}

void FillHistogramsChi2_conv(const Double_t *bestPar)
{
  for (int index = 0; index < minimizer->NDim(); index++)
  {
    for (auto &p : coefficents)
    {
      if (p.first == minimizer->VariableName(index))
      {
        p.second = bestPar[index];
        // cout << p.first << " : " << p.second << endl;
        break;
      }
    }
  }

  ////////////// conversion exp en kev///////////
  for (int i = 0; i < BETA_SIZE + 1; i++)
  {
    Ref_Hist_Multi_F[i]->Reset();
    TTreeReader *Reader = new TTreeReader(Ref_Tree_Multi_F[i]);
    TTreeReaderValue<double> value(*Reader, "Channel");
    TTreeReaderValue<int> label(*Reader, "Label");
    int counter_sum = 0;
    double value_sum = 0;
    while (Reader->Next())
    {
      Ref_Hist_Multi_F[i]->Fill(coefficents[0].second + coefficents[1].second * (*value));
    }
  }

  TF1 *gauss;
  Merged_File->cd();
  TreeHist_Conv_F[0]->Reset();


  //////////CONVOLUTING RESOLUTION//////////
  for (int i = 1; i <= histF->GetNbinsX(); i++)
  {
    double energy = histF->GetBinCenter(i);
    double sigma_resolution = sqrt(pow(coefficents[2].second, 2) + pow(coefficents[3].second * sqrt(energy), 2) + pow(coefficents[4].second * pow(energy, 2), 2));

    gauss = new TF1("gauss", "gaus", 0, 8000);
    gauss->SetNpx(800);
    gauss->SetParameters(histF->GetBinContent(i) / (sigma_resolution * sqrt(2 * M_PI)), energy, sigma_resolution); /// weighted gaussian

    for (int j = 1; j <= TreeHist_Conv_F[0]->GetNbinsY(); j++)
    {
      TreeHist_Conv_F[0]->SetBinContent(i, j, gauss->Eval(histF->GetBinCenter(j)));
    }
    delete gauss;
  }

  TreeHist_Conv_F_av = (TH2D *)TreeHist_Conv_F[0]->Clone("TreeHist_Conv_F1");

  ////////// CONVOLUTING THRESHOLD//////////
  Threshold = new TF1("Threshold", "gaus", 0, 8000);
  Threshold->SetParameters(1 / (coefficents[6].second * sqrt(2 * M_PI)), coefficents[5].second, coefficents[6].second);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    TreeHist_Conv_F[multiplicity] = (TH2D *)TreeHist_Conv_F[0]->Clone(("TreeHist_Conv_F" + to_string(multiplicity)).c_str());
    TreeHist_Conv_F[multiplicity]->Reset();
    for (int i = 1; i <= TreeHist_Conv_F[0]->GetNbinsX(); i++)
    {
      for (int j = 1; j <= TreeHist_Conv_F[0]->GetNbinsY(); j++)
      {
        TreeHist_Conv_F[0]->SetBinContent(i, j, TreeHist_Conv_F[0]->GetBinContent(i, j) * Threshold->Integral(0, j * 10));
      }
    }
  }


  ////////// CONVOLUTING MULTIPLICITY//////////
  if (h_ratio)
  {
    delete h_ratio;
  }
  h_ratio = new TH1D("h_ratio", "h_ratio", 800, 0, 8000);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    if (h_probability[multiplicity])
    {
      delete h_probability[multiplicity];
    }
    h_probability[multiplicity] = new TH1D(("h_probability" + to_string(multiplicity)).c_str(), ("h_probability" + to_string(multiplicity)).c_str(), 800, 0, 8000);
    TreeHist_Conv_F[multiplicity] = (TH2D *)TreeHist_Conv_F[0]->Clone(("TreeHist_Conv_F" + to_string(multiplicity)).c_str());
    TreeHist_Conv_F[multiplicity]->Reset();

    for (int i = 1; i <= TreeHist_Conv_F[0]->GetNbinsX(); i++)
    {
      TH1D *intav = (TH1D *)TreeHist_Conv_F_av->ProjectionY(("av" + to_string(i)).c_str(), i, i + 1);
      TH1D *intap = (TH1D *)TreeHist_Conv_F[0]->ProjectionY(("ap" + to_string(i)).c_str(), i, i + 1);

      double ratio = 0;
      if (intav->Integral(0, 8000) == 0 || intap->Integral(0, 8000) == 0)
      {
        ratio = 0;
      }
      else
      {
        ratio = intap->Integral(0, 8000) / intav->Integral(0, 8000);
      }

      h_ratio->SetBinContent(i, ratio);

      double sum = 0;
      for (int m = multiplicity; m <= 7; m++)
      {
        sum += TMath::Binomial(7, m) * pow(ratio, m) * pow(1 - ratio, 7 - m);
      }
      h_probability[multiplicity]->SetBinContent(i, sum);

      for (int j = 1; j <= TreeHist_Conv_F[0]->GetNbinsY(); j++)
      {
        TreeHist_Conv_F[multiplicity]->SetBinContent(i, j, TreeHist_Conv_F[0]->GetBinContent(i, j) * sum);
      }
    }
  }
}

inline double Chi2TreeHist_conv(const Double_t *bestPar)
{
  double chi2 = 0;
  double pvalue = 0;


  delete TreeHist_Conv_F[0];
  TreeHist_Conv_F[0] = new TH2D("TreeHist_Conv_F", "TreeHist_Conv_F", 800, 0, 8000, 800, 0, 8000);
  TreeHist_Conv_F[0]->Reset();

  FillHistogramsChi2_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    TreeHist_Conv_F_1D[multiplicity] = (TH1D *)TreeHist_Conv_F[multiplicity]->ProjectionY(("TreeHist_Conv_F_1D" + to_string(multiplicity)).c_str(), 0, 8000);
    
    TreeHist_Conv_F_1D[multiplicity]->GetXaxis()->SetRangeUser(0, 8000);
    Ref_Hist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 8000);
    TreeHist_Conv_F_1D[multiplicity]->Scale(Ref_Hist_Multi_F[multiplicity]->Integral() / TreeHist_Conv_F_1D[multiplicity]->Integral());

  for (int bin = 1; bin <= TreeHist_Conv_F_1D[multiplicity]->GetNbinsX(); bin++)
    {
      if (TreeHist_Conv_F_1D[multiplicity]->GetBinContent(bin) < 1)
      {
        TreeHist_Conv_F_1D[multiplicity]->SetBinContent(bin, 1);
      }
    }

    pvalue = Ref_Hist_Multi_F[multiplicity]->Chi2Test(TreeHist_Conv_F_1D[multiplicity], "UW");
    chi2_multiplicity[multiplicity] = (Ref_Hist_Multi_F[multiplicity]->Chi2Test(TreeHist_Conv_F_1D[multiplicity], "UW CHI2/NDF")); 
    chi2 += chi2_multiplicity[multiplicity];
    // cout << "M" << multiplicity << setprecision(5) << "   p-value = " << pvalue << "  Chi2 = " << setprecision(5) << Ref_Hist_Multi_F[multiplicity]->Chi2Test(TreeHist_Conv_F_1D[multiplicity], "WU CHI2/NDF") << endl;
  }

  chi2 /= 7;
cout << "p-value = " << setw(10) << setprecision(5) << pvalue 
     << "  X² = " << setw(6) << setprecision(5) << chi2 
     << "    " << setw(6) << setprecision(5) << bestPar[0] 
     << "   " << setw(6) << setprecision(5) << bestPar[1] 
     << "   " << setw(6) << setprecision(5) << bestPar[2] 
     << "   " << setw(6) << setprecision(5) << bestPar[3] 
     << "   " << setw(6) << setprecision(5) << bestPar[4] 
     << "   " << setw(6) << setprecision(5) << bestPar[5] 
     << "   " << setw(6) << setprecision(5) << bestPar[6] << endl;  return chi2;
}

inline void MakeSiPM_MultiplicityPlots()
{
  counter = 0;
  int result_counter = 0;
  int run = 0;

  

  int countercd = 0;
  double sigma = 0;
  double value = 0;

  Par = {39, 4.9, 0, 3.65, 2.3e-05, 95.496, 12.0801};

  ROOT::Math::Functor functor(&Chi2TreeHist_conv, 7);
  minimizer->SetFunction(functor);
  minimizer->SetLimitedVariable(0, "Calibration_OffSet", Par[0], 5, 0, 50);
  minimizer->SetLimitedVariable(1, "Calibration", Par[1], 0.1, 4.8, 5.);
  minimizer->SetFixedVariable(2, "Resolution_OffSet", Par[2]);
  minimizer->SetLimitedVariable(3, "Resolution_SQRT", Par[3], 1, 3.5, 3.9);
  minimizer->SetLimitedVariable(4, "Resolution_2", Par[4], 1e-6, 1e-5, 4e-5);
  minimizer->SetLimitedVariable(5, "Threshold", Par[5], 5, 80, 150);
  minimizer->SetLimitedVariable(6, "Threshold_STD", Par[6], 10, 5, 16);
  minimizer->SetPrecision(0.1);
  // minimizer->SetTolerance(0.000000001);
  minimizer->SetMaxFunctionCalls(10000000);
  minimizer->SetMaxIterations(10000000);
  bestCHI2 = 1e6;

  // minimizer->Minimize();
  // const double *bestPar = minimizer->X();
  // minimizer->PrintResults();

  const double *bestPar = Par.data();

  chi2 = Chi2TreeHist_conv(bestPar);

  TCanvas *cF = new TCanvas(("SiPM_Multi_F" + to_string(run)).c_str(), ("SiPM_Multi_F" + to_string(run)).c_str(), 1920, 1080);
  cF->Divide(3, 3);
  
  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    // TreeHist_Conv_F_1D[multiplicity] = (TH1D *)TreeHist_Conv_F[multiplicity]->ProjectionY(("TreeHist_Conv_F_1D" + to_string(multiplicity)).c_str(), 0, 8000);
    cF->cd(multiplicity);
    TreeHist_Conv_F_1D[multiplicity]->GetXaxis()->SetRangeUser(0, 8000);
    Ref_Hist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 8000);
    TreeHist_Conv_F_1D[multiplicity]->Scale(Ref_Hist_Multi_F[7]->Integral() / TreeHist_Conv_F_1D[7]->Integral());
    // TreeHist_Conv_F_1D[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
    // Ref_Hist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
    TreeHist_Conv_F_1D[multiplicity]->GetYaxis()->SetRangeUser(1, -1111);
    Ref_Hist_Multi_F[multiplicity]->GetYaxis()->SetRangeUser(1, -1111);
    TreeHist_Conv_F_1D[multiplicity]->SetLineColor(kRed);
    Ref_Hist_Multi_F[multiplicity]->Draw("HIST");
    TreeHist_Conv_F_1D[multiplicity]->Draw("SAME HIST");

    TPaveText *pt = new TPaveText(0.7, 0.65, 1, 0.7, "brNDC");
    pt->SetFillColor(0);
    pt->AddText(("#chi^{2} = " + to_string(chi2_multiplicity[multiplicity])).c_str());
    pt->Draw("SAME");

    
  }

  cF->cd(8);
  graph_calib = new TGraph();
  graph_calib->SetTitle(("Calibration;Energy[keV]; Channel"));
  int counter_cal = 0;

  for (double e = 1; e <= 6000; e += 10)
  {
    counter_cal++;
    value = (e - bestPar[0]) / bestPar[1];
    graph_calib->SetPoint(counter_cal, e, value);
  }
  graph_calib->Draw("AL");

  cF->cd(9);
  graph_resol = new TGraph();
  graph_resol->SetTitle(("Resolution;Energy[keV]; #sigma [keV]"));
  int counter = 0;
  for (double e = 1; e <= 6000; e += 10)
  {
    counter++;
    sigma = sqrt(pow(bestPar[2], 2) + pow(bestPar[3] * sqrt(e), 2) + pow(bestPar[4] * pow(e, 2), 2));
    graph_resol->SetPoint(counter, e, sigma);
  }
  graph_resol->Draw("AL");
  cF->Write();

  TreeHist_Conv_F_1D[1]->GetXaxis()->SetRangeUser(0, 8000);
  /////////// SAVE STEP CONVOLUTING ////////// 
  gStyle->SetOptStat(0);
  TCanvas *c_sim = new TCanvas(("Simulation" + to_string(run)).c_str(), ("Simulation" + to_string(run)).c_str(), 1920, 1080);
  histF->GetXaxis()->SetRangeUser(0, 8000);
  histF->Scale(TreeHist_Conv_F_1D[1]->Integral() / histF->Integral());
  histF->SetLineColor(kBlack);
  histF->Draw("HIST");
  histF->GetXaxis()->SetTitle("Energy [keV]");
  histF->GetYaxis()->SetTitle("Counts / 10keV");
  histF->GetXaxis()->CenterTitle();
  histF->GetYaxis()->CenterTitle();
  c_sim->Write();

  TCanvas *c_sim_conv = new TCanvas(("Simulation_Convoluted" + to_string(run)).c_str(), ("Simulation_Convoluted" + to_string(run)).c_str(), 1920, 1080);
  TLegend *leg_sim_conv = new TLegend(0.6, 0.7, 0.9, 0.9);
  histF->SetLineColor(kBlack);
  histF->Draw("HIST");
  TH1D* histF_conv = (TH1D *)TreeHist_Conv_F_av->ProjectionY("Simulation_Convoluted", 0, 8000);
  histF_conv->GetXaxis()->SetRangeUser(0, 8000);
  histF_conv->Scale(TreeHist_Conv_F_1D[1]->Integral() / histF_conv->Integral());
  histF_conv->SetLineColor(kRed);
  histF_conv->Draw("SAME HIST");
  leg_sim_conv->AddEntry(histF, "Simulation", "l");
  leg_sim_conv->AddEntry(histF_conv, "Simulation Convoluted", "l");
  leg_sim_conv->Draw("SAME");
  histF->GetXaxis()->SetTitle("Energy [keV]");
  histF->GetYaxis()->SetTitle("Counts / 10keV");
  histF->GetXaxis()->CenterTitle();
  histF->GetYaxis()->CenterTitle();
  c_sim_conv->Write();

  TCanvas *c_sim_conv_th = new TCanvas(("Simulation_Convoluted_Threshold" + to_string(run)).c_str(), ("Simulation_Convoluted_Threshold" + to_string(run)).c_str(), 1920, 1080);
  TLegend *leg_sim_conv_th = new TLegend(0.6, 0.7, 0.9, 0.9);
  histF_conv->SetLineColor(kBlack);
  histF_conv->Draw("HIST");
  TreeHist_Conv_F_1D[1]->Scale(0.96);
  TreeHist_Conv_F_1D[1]->SetLineColor(kRed);
  TreeHist_Conv_F_1D[1]->Draw("SAME HIST");
  leg_sim_conv_th->AddEntry(histF_conv, "Simulation Convoluted", "l");
  leg_sim_conv_th->AddEntry(TreeHist_Conv_F_1D[1], "Simulation Convoluted Threshold", "l");
  leg_sim_conv_th->Draw("SAME");
  histF_conv->GetXaxis()->SetTitle("Energy [keV]");
  histF_conv->GetYaxis()->SetTitle("Counts / 10keV");
  histF_conv->GetXaxis()->CenterTitle();
  histF_conv->GetYaxis()->CenterTitle();
  c_sim_conv_th->Write();

  TCanvas *c_th = new TCanvas(("Threshold" + to_string(run)).c_str(), ("Threshold" + to_string(run)).c_str(), 1920, 1080);
  Threshold->SetNpx(8000);
  Threshold->Draw("L");
  Threshold->DrawIntegral("L");
  Threshold->GetXaxis()->SetTitle("Energy [keV]");
  Threshold->GetYaxis()->SetTitle("Counts");
  Threshold->GetXaxis()->CenterTitle();
  Threshold->GetYaxis()->CenterTitle();
  c_th->Write();

  // TCanvas *c_th_int = new TCanvas(("Threshold_Probability" + to_string(run)).c_str(), ("Threshold_Probability" + to_string(run)).c_str(), 1920, 1080);
  // Threshold->SetNpx(8000);
  
  // Threshold->GetXaxis()->SetTitle("Energy [keV]");
  // Threshold->GetYaxis()->SetTitle("Probability");
  // Threshold->GetXaxis()->CenterTitle();
  // Threshold->GetYaxis()->CenterTitle();
  // c_th_int->Write();

  TCanvas *c_trig_100 = new TCanvas(("Simulation_Convoluted_100keV" + to_string(run)).c_str(), ("Simulation_Convoluted_100keV" + to_string(run)).c_str(), 1920, 1080);
  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  TH1D* hist_trigger_av = (TH1D *)TreeHist_Conv_F_av->ProjectionY("Trigger_Probability_100keV_av", 9, 10);
  TH1D *hist_trigger_ap = (TH1D *)TreeHist_Conv_F[1]->ProjectionY("Trigger_Probability_100keV_ap", 9, 10);
  hist_trigger_av->SetLineColor(kBlack);
  hist_trigger_av->Draw("HIST");
  hist_trigger_ap->SetLineColor(kRed);
  hist_trigger_ap->Draw("SAME HIST");
  leg->AddEntry(hist_trigger_av, "Simulation Convoluted without Threshold", "l");
  leg->AddEntry(hist_trigger_ap, "Simulation Convoluted with Threshold", "l");
  leg->Draw("SAME");
  hist_trigger_av->GetXaxis()->SetTitle("Energy [keV]");
  hist_trigger_av->GetYaxis()->SetTitle("Counts / 10keV");
  hist_trigger_av->GetXaxis()->CenterTitle();
  hist_trigger_av->GetYaxis()->CenterTitle();
  c_trig_100->Write();

  TCanvas *c_trig_1000 = new TCanvas(("Simulation_Convoluted_1000keV" + to_string(run)).c_str(), ("Simulation_Convoluted_1000keV" + to_string(run)).c_str(), 1920, 1080);
  TLegend *leg_1000 = new TLegend(0.6, 0.7, 0.9, 0.9);
  TH1D* hist_trigger_av_1000 = (TH1D *)TreeHist_Conv_F_av->ProjectionY("Trigger_Probability_100keV_av", 99, 100);
  TH1D *hist_trigger_ap_1000 = (TH1D *)TreeHist_Conv_F[1]->ProjectionY("Trigger_Probability_100keV_ap", 99, 100);
  hist_trigger_av->SetLineColor(kBlack);
  hist_trigger_av->Draw("HIST");
  hist_trigger_ap->SetLineColor(kRed);
  hist_trigger_ap->Draw("SAME HIST");
  leg_1000->AddEntry(hist_trigger_av, "Simulation Convoluted without Threshold", "l");
  leg_1000->AddEntry(hist_trigger_ap, "Simulation Convoluted with Threshold", "l");
  leg_1000->Draw("SAME");
  hist_trigger_av->GetXaxis()->SetTitle("Energy [keV]");
  hist_trigger_av->GetYaxis()->SetTitle("Counts / 10keV");
  hist_trigger_av->GetXaxis()->CenterTitle();
  hist_trigger_av->GetYaxis()->CenterTitle();
  c_trig_1000->Write();

  TCanvas *c_trig = new TCanvas(("P_{trig}" + to_string(run)).c_str(), ("Trigger_Probability" + to_string(run)).c_str(), 1920, 1080);
  h_ratio->Draw("HIST");
  h_ratio->GetXaxis()->SetTitle("Energy [keV]");
  h_ratio->GetYaxis()->SetTitle("Probability");
  h_ratio->GetXaxis()->CenterTitle();
  h_ratio->GetYaxis()->CenterTitle();
  c_trig->Write();

  TCanvas *c_prob = new TCanvas(("Trigger_Probability_Multiplicity" + to_string(run)).c_str(), ("P_{trig} with Multiplicity" + to_string(run)).c_str(), 1920, 1080);
  TLegend *leg_prob = new TLegend(0.1, 0.7, 0.48, 0.9);
  leg->SetHeader("Multiplicity");
  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    h_probability[multiplicity]->SetLineColor(multiplicity);
    h_probability[multiplicity]->Draw("SAME HIST");
    leg->AddEntry(h_probability[multiplicity], ("M" + to_string(multiplicity)+"+").c_str(), "l");
  }
  leg_prob->Draw("SAME");
  h_probability[1]->GetXaxis()->SetTitle("Energy [keV]");
  h_probability[1]->GetYaxis()->SetTitle("Counts / 10keV");
  h_probability[1]->GetXaxis()->CenterTitle();
  h_probability[1]->GetYaxis()->CenterTitle();
  c_prob->Write();
}


void MakeSiPM_ParameterPlots()
{
  int run = 0;
  const double *bestPar;

  //////PARMATERS///////
  TH1D* hist_PARAM_0[BETA_SIZE + 1];
  TH1D* hist_PARAM_1[BETA_SIZE + 1];
  TH1D* hist_PARAM_2[BETA_SIZE + 1];
  for (int i = 0; i < BETA_SIZE + 1; i++)
  {
    hist_PARAM_0[i] = (TH1D *)Ref_Hist_Multi_F[i]->Clone(("hist_PARAM0_" + to_string(i)).c_str());
    hist_PARAM_0[i]->Reset();
    hist_PARAM_1[i] = (TH1D *)Ref_Hist_Multi_F[i]->Clone(("hist_PARAM1_" + to_string(i)).c_str());
    hist_PARAM_1[i]->Reset();
    hist_PARAM_2[i] = (TH1D *)Ref_Hist_Multi_F[i]->Clone(("hist_PARAM2_" + to_string(i)).c_str());
    hist_PARAM_2[i]->Reset();
  }


  TCanvas *cF_offset = new TCanvas(("SiPM_OFFSET_F" + to_string(run)).c_str(), ("SiPM_OFFSET_F" + to_string(run)).c_str(), 1920, 1080);
  cF_offset->Divide(3, 3);
  TCanvas *cF_calib = new TCanvas(("SiPM_CALIB_F" + to_string(run)).c_str(), ("SiPM_CALIB_F" + to_string(run)).c_str(), 1920, 1080);
  cF_calib->Divide(3, 3);
  TCanvas *cF_resol_sqrt = new TCanvas(("SiPM_RESOL_SQRT_F" + to_string(run)).c_str(), ("SiPM_RESOL_SQRT_F" + to_string(run)).c_str(), 1920, 1080);
  cF_resol_sqrt->Divide(3, 3);
  TCanvas *cF_resol_2 = new TCanvas(("SiPM_RESOL_2_F" + to_string(run)).c_str(), ("SiPM_RESOL_2_F" + to_string(run)).c_str(), 1920, 1080);
  cF_resol_2->Divide(3, 3);
  TCanvas *cF_th = new TCanvas(("SiPM_TH_F" + to_string(run)).c_str(), ("SiPM_TH_F" + to_string(run)).c_str(), 1920, 1080);
  cF_th->Divide(3, 3);
  TCanvas *cF_th_std = new TCanvas(("SiPM_TH_STD_F" + to_string(run)).c_str(), ("SiPM_TH_STD_F" + to_string(run)).c_str(), 1920, 1080);
  cF_th_std->Divide(3, 3);

  double MAX_OFFSET[BETA_SIZE + 1];
  double MAX_CALIB[BETA_SIZE + 1];
  double MAX_RESOL_SQRT[BETA_SIZE + 1];
  double MAX_RESOL_2[BETA_SIZE + 1];
  double MAX_TH[BETA_SIZE + 1];
  double MAX_TH_STD[BETA_SIZE + 1];

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    hist_PARAM_0[multiplicity] = (TH1D*)Ref_Hist_Multi_F[multiplicity]->Clone(("hist_PARAM0_" + to_string(multiplicity)).c_str());
    hist_PARAM_0[multiplicity]->SetLineColor(kRed);

    cF_offset->cd(multiplicity);
    MAX_OFFSET[multiplicity] = hist_PARAM_0[multiplicity]->GetMaximum();
    hist_PARAM_0[multiplicity]->Draw("SAME HIST");

    cF_calib->cd(multiplicity);
    MAX_CALIB[multiplicity] = hist_PARAM_0[multiplicity]->GetMaximum();
    hist_PARAM_0[multiplicity]->Draw("SAME HIST");

    hist_PARAM_0[multiplicity] = (TH1D*)TreeHist_Conv_F_1D[multiplicity]->Clone(("hist_PARAM0_" + to_string(multiplicity)).c_str());
    cF_resol_sqrt->cd(multiplicity);
    MAX_RESOL_SQRT[multiplicity] = hist_PARAM_0[multiplicity]->GetMaximum();
    hist_PARAM_0[multiplicity]->Draw("SAME HIST");

    cF_resol_2->cd(multiplicity);
    MAX_RESOL_2[multiplicity] = hist_PARAM_0[multiplicity]->GetMaximum();
    hist_PARAM_0[multiplicity]->Draw("SAME HIST");

    cF_th->cd(multiplicity);
    MAX_TH[multiplicity] = hist_PARAM_0[multiplicity]->GetMaximum();
    hist_PARAM_0[multiplicity]->Draw("SAME HIST");

    cF_th_std->cd(multiplicity);
    MAX_TH_STD[multiplicity] = hist_PARAM_0[multiplicity]->GetMaximum();
    hist_PARAM_0[multiplicity]->Draw("SAME HIST");
  }

  ////OFFSET////
  Par[0] = 0;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_offset->cd(multiplicity);
    hist_PARAM_1[multiplicity] = (TH1D *)Ref_Hist_Multi_F[multiplicity]->Clone(("hist_PARAM1_" + to_string(multiplicity)).c_str());
    hist_PARAM_1[multiplicity]->SetLineColor(8);
    hist_PARAM_1[multiplicity]->Draw("SAME HIST");
  }

  Par[0] = 80;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_offset->cd(multiplicity);
    hist_PARAM_2[multiplicity] = (TH1D *)Ref_Hist_Multi_F[multiplicity]->Clone(("hist_PARAM2_" + to_string(multiplicity)).c_str());
    hist_PARAM_2[multiplicity]->SetLineColor(9);
    hist_PARAM_2[multiplicity]->Draw("SAME HIST");
  }

  Par[0] = 39;
  cF_offset->Write();

  ////CALIBRATION////
  Par[1] = 4.8;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_calib->cd(multiplicity);
    hist_PARAM_1[multiplicity] = (TH1D *)Ref_Hist_Multi_F[multiplicity]->Clone(("hist_PARAM1_" + to_string(multiplicity)).c_str());
    hist_PARAM_1[multiplicity]->SetLineColor(8);
    hist_PARAM_1[multiplicity]->Draw("SAME HIST");
  }

  Par[1] = 5;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_calib->cd(multiplicity);
    hist_PARAM_2[multiplicity] = (TH1D *)Ref_Hist_Multi_F[multiplicity]->Clone(("hist_PARAM2_" + to_string(multiplicity)).c_str());
    hist_PARAM_2[multiplicity]->SetLineColor(9);
    hist_PARAM_2[multiplicity]->Draw("SAME HIST");
  }
  Par[1] = 4.9;
  cF_calib->Write();

  ////RESOLUTION SQRT////
  Par[3] = 2.6;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_resol_sqrt->cd(multiplicity);
    hist_PARAM_1[multiplicity] = (TH1D *)TreeHist_Conv_F_1D[multiplicity]->Clone(("hist_PARAM1_" + to_string(multiplicity)).c_str());
    hist_PARAM_1[multiplicity]->SetLineColor(8);
    hist_PARAM_1[multiplicity]->Scale(MAX_RESOL_SQRT[multiplicity] / hist_PARAM_1[multiplicity]->GetMaximum());
    hist_PARAM_1[multiplicity]->Draw("SAME HIST");
  }

  Par[3] = 4.6;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_resol_sqrt->cd(multiplicity);
    hist_PARAM_2[multiplicity] = (TH1D *)TreeHist_Conv_F_1D[multiplicity]->Clone(("hist_PARAM2_" + to_string(multiplicity)).c_str());
    hist_PARAM_2[multiplicity]->SetLineColor(9);
    hist_PARAM_2[multiplicity]->Scale(MAX_RESOL_SQRT[multiplicity] / hist_PARAM_2[multiplicity]->GetMaximum());
    hist_PARAM_2[multiplicity]->Draw("SAME HIST");
  }

  Par[3] = 3.65;
  cF_resol_sqrt->Write();

  ////RESOLUTION 2////
  Par[4] = 0.1e-5;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_resol_2->cd(multiplicity);
    hist_PARAM_1[multiplicity] = (TH1D *)TreeHist_Conv_F_1D[multiplicity]->Clone(("hist_PARAM1_" + to_string(multiplicity)).c_str());
    hist_PARAM_1[multiplicity]->SetLineColor(8);
    // hist_PARAM_1[multiplicity]->Scale(MAX_RESOL_2[multiplicity] / hist_PARAM_1[multiplicity]->GetMaximum());
    hist_PARAM_1[multiplicity]->Draw("SAME HIST");
  }

  Par[4] = 4e-5;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_resol_2->cd(multiplicity);
    hist_PARAM_2[multiplicity] = (TH1D *)TreeHist_Conv_F_1D[multiplicity]->Clone(("hist_PARAM2_" + to_string(multiplicity)).c_str());
    hist_PARAM_2[multiplicity]->SetLineColor(9);
    // hist_PARAM_2[multiplicity]->Scale(MAX_RESOL_2[multiplicity] / hist_PARAM_2[multiplicity]->GetMaximum());
    hist_PARAM_2[multiplicity]->Draw("SAME HIST");
  }

  Par[4] = 2.3e-05;
  cF_resol_2->Write();

  ////THRESHOLD////
  Par[5] = 80;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_th->cd(multiplicity);
    hist_PARAM_1[multiplicity] = (TH1D *)TreeHist_Conv_F_1D[multiplicity]->Clone(("hist_PARAM1_" + to_string(multiplicity)).c_str());
    hist_PARAM_1[multiplicity]->SetLineColor(8);
    hist_PARAM_1[multiplicity]->Scale(MAX_TH[multiplicity] / hist_PARAM_1[multiplicity]->GetMaximum());
    hist_PARAM_1[multiplicity]->Draw("SAME HIST");
  }

  Par[5] = 150;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_th->cd(multiplicity);
    hist_PARAM_2[multiplicity] = (TH1D *)TreeHist_Conv_F_1D[multiplicity]->Clone(("hist_PARAM2_" + to_string(multiplicity)).c_str());
    hist_PARAM_2[multiplicity]->SetLineColor(9);
    hist_PARAM_2[multiplicity]->Scale(MAX_TH[multiplicity] / hist_PARAM_2[multiplicity]->GetMaximum());
    hist_PARAM_2[multiplicity]->Draw("SAME HIST");
  }

  Par[5] = 95.496;
  cF_th->Write();

  ////THRESHOLD STD////
  Par[6] = 16;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_th_std->cd(multiplicity);
    hist_PARAM_1[multiplicity] = (TH1D *)TreeHist_Conv_F_1D[multiplicity]->Clone(("hist_PARAM1_" + to_string(multiplicity)).c_str());
    hist_PARAM_1[multiplicity]->SetLineColor(8);
    hist_PARAM_1[multiplicity]->Scale(MAX_TH_STD[multiplicity] / hist_PARAM_1[multiplicity]->GetMaximum());
    hist_PARAM_1[multiplicity]->Draw("SAME HIST");
  }

  Par[6] = 8.0801;
  bestPar = Par.data();
  chi2 = Chi2TreeHist_conv(bestPar);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF_th_std->cd(multiplicity);
    hist_PARAM_2[multiplicity] = (TH1D *)TreeHist_Conv_F_1D[multiplicity]->Clone(("hist_PARAM2_" + to_string(multiplicity)).c_str());
    hist_PARAM_2[multiplicity]->SetLineColor(9);
    hist_PARAM_2[multiplicity]->Scale(MAX_TH_STD[multiplicity] / hist_PARAM_2[multiplicity]->GetMaximum());
    hist_PARAM_2[multiplicity]->Draw("SAME HIST");
  }

  Par[6] = 12.0801;
  cF_th_std->Write();


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

vector<int> LabelSimToExp(int label_sim)
{
  vector<int> label_exp(4);

  if (label_sim > 0)
  {
    label_exp[0] = 110 - label_sim % 10;
    label_exp[1] = 120 - label_sim % 10;
    label_exp[2] = 130 - label_sim % 10;
    label_exp[3] = 140 - label_sim % 10;
  }
  else
  {
    label_exp[0] = 150 + label_sim % 10;
    label_exp[1] = 160 + label_sim % 10;
    label_exp[2] = 170 + label_sim % 10;
    label_exp[3] = 180 + label_sim % 10;
  }

  for (int &label : label_exp)
  {
    for (int i = 0; i < SIGNAL_MAX; i++)
    {
      if (detectorInfo[i] == label)
      {
        label = detectorCoder[i];
        break;
      }
    }
  }

  return label_exp;
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

  file.close();

  // FOR CHECK
  //  int index = 0;
  //  for (auto &name : detectorName)
  //  {
  //    cout << name << "       " << SiliconCalib_Fermi[index] << endl;
  //    index++;
  //  }

  return error;
}

int InitCalibSiPM()
{
  int error = 0;
  ///////FOR SIPM
  clock_t start = clock(), Current;
  TFile *Calibration_File = new TFile("../../../../../../mnt/hgfs/shared-2/32Ar_a1_b0_1.root");
  TTree *tree = (TTree *)Calibration_File->Get("Tree");
  TTreeReader *Reader = new TTreeReader(tree);
  TTreeReaderArray<int> DetectorCode(*Reader, "Silicon_Detector_Code");
  TTreeReaderValue<double> PlasticEnergy(*Reader, "PlasticScintillator_Deposit_Energy");
  TTreeReaderArray<double> SiliconEnergy(*Reader, "Silicon_Detector_Deposit_Energy");

  TFile *f = new TFile("f.root", "RECREATE");
  double Energy_SiPM = 0;
  TTree *Filtered_Tree_F = new TTree("Filtered_Tree_F", "Filtered_Tree_F");
  Filtered_Tree_F->Branch("PlasticScintillator_Deposit_Energy", &Energy_SiPM, "PlasticScintillator_Deposit_Energy/D");
  TTree *Filtered_Tree_GT = new TTree("Filtered_Tree_GT", "Filtered_Tree_GT");
  Filtered_Tree_GT->Branch("PlasticScintillator_Deposit_Energy", &Energy_SiPM, "PlasticScintillator_Deposit_Energy/D");

  bool fermi = false;
  int TotalEntries = Reader->GetEntries();
  while (Reader->Next())
  {

    ULong64_t cEntry = Reader->GetCurrentEntry();
    if (cEntry % 100 == 0 && cEntry > 0)
    {
      ProgressBar(cEntry, TotalEntries, start, Current, "<SAM> Reading Simulation File | ");
    }

    bool interstrip = false;
    bool fermi = false;
    bool gt = false;
    for (int strip_A : DetectorCode)
    {
      for (int strip_B : DetectorCode)
      {
        if ((strip_A / 10 == strip_B / 10) && (strip_A != strip_B))
        {
          interstrip = true;
          break;
        }
      }
    }

    if (!interstrip && *PlasticEnergy > 0)
    {
      for (int i = 0; i < DetectorCode.GetSize(); i++)
      {
        for (int label : LabelSimToExp(DetectorCode[i]))
        {
          if (Is_F(SiliconEnergy[i], label))
          { // histF->Scale(histGT->Integral() / histF->Integral());

            fermi = true;
            break;
          }
          if (Is_GT(SiliconEnergy[i], label))
          {
            gt = true;
            break;
          }
        }
      }

      if (fermi)
      {
        Energy_SiPM = (*PlasticEnergy) / 1000;
        Filtered_Tree_F->Fill();
        histF->Fill(Energy_SiPM);
      }
      else if (gt)
      {
        Energy_SiPM = (*PlasticEnergy) / 1000;
        Filtered_Tree_GT->Fill();
        histGT->Fill(Energy_SiPM);
      }
    }
  }
  Reader_calib_sipm_F = new TTreeReader(Filtered_Tree_F);
  Reader_calib_sipm_GT = new TTreeReader(Filtered_Tree_GT);
  TCanvas *cFGT = new TCanvas("cFGT", "cFGT", 800, 600);
  cFGT->cd();
  // histF->Scale(histGT->Integral() / histF->Integral());
  histGT->SetLineColor(kRed);
  histGT->Draw("HIST");
  histF->Draw("SAME HIST");
  cFGT->Write();
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