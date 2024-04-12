#ifndef MERGER_HH
#define MERGER_HH
//
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
TH1D *Ref_Hist_Surface_F[BETA_SIZE + 1];
TH1D *TreeHist_F[BETA_SIZE + 1][BETA_SIZE + 1];
TH1D *TreeHist_Multi_F[BETA_SIZE + 1];
TH1D *TreeHist_Surface_F[BETA_SIZE + 1];
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
TTree* Ref_Tree_F[BETA_SIZE + 1][BETA_SIZE + 1];
TTree* Ref_Tree_Multi_F[BETA_SIZE + 1];
TTree* Ref_Tree_Surface_F[BETA_SIZE + 1];
vector<array<double, 7>> bestPar_sipm;
int current_sipm;
TH1D *checker_si[SIGNAL_MAX];
TF1 *linear = new TF1("linear", "[0]*x");
TH1D *histF = new TH1D("histF", "histF", 1000, 0, 10000);
TH1D *histGT = new TH1D("histGT", "histGT", 1000, 0, 10000);
TGraph *graph_res[BETA_SIZE+1];
TGraph *graph_cal[BETA_SIZE+1];
TF1 *dark;
TFile* dark_file;

int Verbose = 0;

bool Calibration = false;

int counter;

double bin_time;
int n_bin_time;
double start_time;
double stop_time;

void generate_combinations(std::vector<pair<int, double>>& elements, int combination_size, int start, std::vector<pair<int, double>>& current_combination, double Channel, double Label, TTree* Tree) {
    if (current_combination.size() == combination_size) {
        for (auto i : current_combination) {
          Label = i.first;
          Channel = i.second;
          Tree->Fill();
        }
    } else {
        for (int i = start; i < elements.size(); ++i) {
            current_combination.push_back(elements[i]);
            generate_combinations(elements, combination_size, i + 1, current_combination, Channel, Label, Tree);
            current_combination.pop_back();
        }
    }
}

void generate_combinations_sim(std::vector<double>& elements, int combination_size, int start, std::vector<double>& current_combination, TH1D* Hist) {
    if (current_combination.size() == combination_size) {
      double sum = 0;
        for (auto i : current_combination) {
          sum += i;
        }
        Hist->Fill(sum/combination_size);
        
    } else {
        for (int i = start; i < elements.size(); ++i) {
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

// void FillHistogramsChi2(const Double_t *bestPar)
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


//     Reader_calib_sipm_F->Restart();
//     TTreeReaderValue<double> Tree_SiPM(*Reader_calib_sipm_F, "PlasticScintillator_Deposit_Energy");

//     ////////////// conversion exp en kev///////////
//     for (int i = 0; i < BETA_SIZE+1; i++)
//     {
//       Ref_Hist_F[i]->Reset();
//       TTreeReader *Reader = new TTreeReader(Ref_Tree_F[i]);
//       TTreeReaderValue<double> value(*Reader, "Channel");
//       while (Reader->Next())
//       {
//         Ref_Hist_F[i]->Fill(coefficents[0].second + coefficents[1].second * (*value));
//       }
//     }
//     /////////////////////////////////////////////////

//     vector<double> vec_e;
//     double e_sum;
//     double sigma_resolution = 0;
//     double sigma_resolution_sum = 0;
//     bool negative_res = false;

//     for (int i = 0; i < 10; i++)
//     {
//       Reader_calib_sipm_F->Restart();
//       while (Reader_calib_sipm_F->Next())
//       {
        
//         double energy = (*Tree_SiPM);
//         double energy_sum = (*Tree_SiPM);
//         int real_multiplicity = 0;
//         vec_e.clear();
//         e_sum = 0;

//         for (int sipm = 1; sipm <= 7; sipm++)
//         {
//           sigma_resolution = sqrt(pow(coefficents[2].second, 2) + pow(coefficents[3].second * sqrt(energy), 2) + pow(coefficents[4].second * pow(energy, 2), 2));
//           sigma_resolution_sum = sqrt(pow(coefficents[2].second, 2) + pow(coefficents[3].second * sqrt(energy), 2));

//           normal_distribution<> resolution(0, sigma_resolution);
//           normal_distribution<> resolution_sum(0, sigma_resolution_sum);

//           energy += resolution(gen);
//           energy_sum += resolution_sum(gen);

//           normal_distribution<> threshold(coefficents[5].second, coefficents[6].second);
//           if (energy > threshold(gen))
//           {
//             real_multiplicity++;
//             vec_e.push_back(energy);
//             e_sum += energy_sum;
//           }
//         }

//         for (int multi = 1; multi <= real_multiplicity; multi++)
//         {
//           for (int sipm = 1; sipm <= real_multiplicity; sipm++)
//           {
//             TreeHist_F[multi]->Fill(vec_e[sipm-1]);
//           }
//         }
//         if (real_multiplicity == 7)
//         {
//           TreeHist_F[0]->Fill(e_sum / 7);
//         }
//       }
//     }

//   }
//   // if (negative_res)
//   // {
//   //   for (int multi = 1; multi <= 7; multi++)
//   //       {
//   //           TreeHist[multi]->Reset();

//   //       }
//   // }


// inline double Chi2TreeHist(const Double_t *bestPar)
// {
//   double chi2 = 0;
//   double chi2_std = 0;
//   chi2_vec.clear();

  

//   for (int multiplicity = 0; multiplicity <= 7; multiplicity++)
//   {
//     TreeHist_F[multiplicity] = (TH1D *)Ref_Hist_F[multiplicity]->Clone(("SiPM_F" + to_string(multiplicity)).c_str());
//     TreeHist_F[multiplicity]->Reset();
//   }

//   FillHistogramsChi2(bestPar);

//   int min = 250;
//   int max = 2000;

//   for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
//   {
//     // TreeHist_F[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
//     // Ref_Hist_F[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
//     // chi2_vec.push_back(Ref_Hist_F[multiplicity]->Chi2Test(TreeHist_F[multiplicity], " CHI2/NDF"));

//     // Ref_Hist_F[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
//     // TreeHist_F[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
//     // chi2_vec.push_back(Ref_Hist_F[multiplicity]->Chi2Test(TreeHist_F[multiplicity], "CHI2/NDF"));
//     // cout << Ref_Hist_F[multiplicity]->Chi2Test(TreeHist_F[multiplicity], "CHI2/NDF") << endl;

//     Ref_Hist_F[multiplicity]->GetXaxis()->SetRangeUser(0, 7000);
//     TreeHist_F[multiplicity]->GetXaxis()->SetRangeUser(0, 7000);
//     chi2_vec.push_back(Ref_Hist_F[multiplicity]->Chi2Test(TreeHist_F[multiplicity], "CHI2/NDF"));
//     cout << Ref_Hist_F[multiplicity]->Chi2Test(TreeHist_F[multiplicity], "CHI2/NDF") << endl;

//     // Ref_Hist_F[multiplicity]->GetXaxis()->SetRangeUser(1000, 1300);
//     // TreeHist_F[multiplicity]->GetXaxis()->SetRangeUser(1000, 1300);
//     // chi2_vec.push_back(Ref_Hist_F[multiplicity]->Chi2Test(TreeHist_F[multiplicity], " CHI2/NDF"));
//   }

//   double sum = accumulate(chi2_vec.begin(), chi2_vec.end(), 0.0);
//   chi2 = sum / chi2_vec.size();

//   double sq_sum = inner_product(chi2_vec.begin(), chi2_vec.end(), chi2_vec.begin(), 0.0);
//   chi2_std = sqrt(sq_sum / chi2_vec.size() - chi2 * chi2);

//   cout << chi2 << " +/- " << chi2_std << "    {";
//   for (size_t i = 0; i < 7; ++i)
//   {
//     cout << bestPar[i] << ", ";
//   }
//   cout << "};";
//   if (bestCHI2 > chi2)
//   {
//     bestCHI2 = chi2;
//     cout << "    BEST" << endl;
//   }
//   else
//   {
//     cout << endl;
//   }
//   return chi2;
// }

// void MakeSiPMCalibration(int run)
// {
//   counter = 0;
//   int result_counter = 0;

//   //////////////////// MINIMIZER ////////////////////
//   ////// First

//   ROOT::Math::Functor functor(&Chi2TreeHist, 7);
//   minimizer->SetFunction(functor);
//   minimizer->SetFixedVariable(0, "Calibration_OffSet", 9);
//   minimizer->SetLimitedVariable(1, "Calibration", 4.5, 0.1, 4., 5.5);
//   minimizer->SetFixedVariable(2, "Resolution_OffSet", 0);
//   minimizer->SetLimitedVariable(3, "Resolution_SQRT", 0, 1, 0, 5);
//   minimizer->SetLimitedVariable(4, "Resolution_2", 0, 1e6, 0, 5e-5);
//   minimizer->SetFixedVariable(5, "Threshold", 81.7);
//   minimizer->SetFixedVariable(6, "Threshold_STD", 18);
//   minimizer->SetPrecision(0.1);
//   minimizer->SetTolerance(0.1);
//   minimizer->SetMaxFunctionCalls(10000000);
//   minimizer->SetMaxIterations(10000000);

//   // minimizer->Minimize();
//   // const double *bestPar = minimizer->X();
//   // minimizer->PrintResults();

  

//   // double bestPar[6] = {0.21575, 1.96249, 1.62856e-06, 0., 13.9911, 2.74464};          /// BEST    Positive offset
//   // double bestPar[7] = {0.5, 0.203913, 0, 1.85754, 6.62335e-05, 13.3651, 2.12369};       //BEST    Negative offset propop, res, res2, offset, th, th_sigma

//   // double bestPar[7] = {-3, 0.200901, 0., 1.9503, 5e-5, 14., 3.26826, }; //// BEST FOR MULTIPLICITY*
//   double bestPar[7] = {9, 5.019, 0, 3.2545, 7e-6, 81.6739, 17.9915, }; //// BEST FOR MULTIPLICITY


//   double chi2 = Chi2TreeHist(bestPar);

//   TCanvas *cF = new TCanvas(("SiPM_F" + to_string(run)).c_str(), ("SiPM_F" + to_string(run)).c_str(), 1920, 1080);
//   cF->Divide(3, 3);
//   TCanvas *cGT = new TCanvas(("SiPM_GT" + to_string(run)).c_str(), ("SiPM_GT" + to_string(run)).c_str(), 1920, 1080);
//   cGT->Divide(3, 3);
//   for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
//   {
//     cF->cd(multiplicity);
//     TreeHist_F[7]->GetXaxis()->SetRangeUser(0, 8000);
//     Ref_Hist_F[7]->GetXaxis()->SetRangeUser(0, 8000);
//     TreeHist_F[multiplicity]->Scale(Ref_Hist_F[7]->Integral() / TreeHist_F[7]->Integral());

//     Ref_Hist_F[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
//     TreeHist_F[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
//     Ref_Hist_F[multiplicity]->Draw("HIST");
//     TreeHist_F[multiplicity]->SetLineColor(kRed);
//     TreeHist_F[multiplicity]->Draw("SAME HIST");

//     TPaveText *pt = new TPaveText(0.7, 0.65, 1, 0.7, "brNDC");
//     pt->SetFillColor(0);
//     pt->AddText(("#chi^{2} = " + to_string(chi2_vec[multiplicity - 1])).c_str());
//     pt->Draw("SAME");

//   }
//   // cF->cd(8);
//   // TreeHist_F[7]->GetXaxis()->SetRangeUser(0, 6000);
//   // Ref_Hist_F[7]->GetXaxis()->SetRangeUser(0, 6000);
//   // Ref_Hist_F[0]->Rebin(2);
//   // Ref_Hist_F[0]->Scale(Ref_Hist_F[7]->Integral() / TreeHist_F[7]->Integral());
//   // Ref_Hist_F[0]->SetLineColor(kBlack);
//   // Ref_Hist_F[0]->Draw("HIST");
//   // TreeHist_F[0]->Scale(Ref_Hist_F[7]->Integral() / TreeHist_F[7]->Integral() / 28);
//   // TreeHist_F[0]->SetLineColor(kRed);
//   // TreeHist_F[0]->Draw("SAME HIST");

//   cF->cd(9);
//   graph_res = new TGraph();
//   graph_res->SetTitle("Resolution;Energy[keV]; #sigma");
//   graph_sqrt = new TGraph();
//   int counter = 0;

//   for (double e = 1; e <= 6000; e += 10)
//   {
//     counter++;
//     double sigma = sqrt(pow(bestPar[2], 2) + pow(bestPar[3] * sqrt(e), 2) + pow(bestPar[4] * pow(e, 2), 2));
//     graph_res->SetPoint(counter, e, sigma);

//     graph_sqrt->SetPoint(counter, e, sqrt(pow(bestPar[2], 2) + pow(bestPar[3] * sqrt(e), 2)));
//   }
//   graph_res->SetLineColor(kRed);
//   graph_res->Draw("AL");
//   graph_sqrt->SetLineColor(kBlack);
//   graph_sqrt->Draw("SAME");

//   Merged_File->cd();
//   cF->Write();
//   cGT->Write();

// }
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

  vector<double> vec_e;
  double e_sum;
  double sigma_resolution = 0;
  double sigma_resolution_sum = 0;
  bool negative_res = false;

  for (int i = 0; i < 10; i++)
  {
    Reader_calib_sipm_F->Restart();
    while (Reader_calib_sipm_F->Next())
    {
      double energy = (*Tree_SiPM);
      int real_multiplicity = 0;
      vec_e.clear();
      e_sum = 0;
      sigma_resolution = sqrt(pow(coefficents[2].second, 2) + pow(coefficents[3].second * sqrt(energy), 2) + pow(coefficents[4].second * pow(energy, 2), 2));

      normal_distribution<> resolution(0, sigma_resolution);
      energy += resolution(gen);
      normal_distribution<> threshold(coefficents[5].second, coefficents[6].second);

      if (energy > threshold(gen))
      {
        TreeHist_F[1][current_sipm]->Fill(energy);
      }
    }
  }
}

inline double Chi2TreeHist(const Double_t *bestPar)
{
  double chi2 = 0;
  double chi2_std = 0;
  chi2_vec.clear();

  TreeHist_F[1][current_sipm] = (TH1D *)Ref_Hist_F[1][current_sipm]->Clone(("SiPM_FF" + to_string(current_sipm)).c_str());
  TreeHist_F[1][current_sipm]->Reset();

  FillHistogramsChi2(bestPar);

  int min = 250;
  int max = 2000;

  Ref_Hist_F[1][current_sipm]->GetXaxis()->SetRangeUser(0, 8000);
  TreeHist_F[1][current_sipm]->GetXaxis()->SetRangeUser(0, 8000);
  chi2_vec.push_back(Ref_Hist_F[1][current_sipm]->Chi2Test(TreeHist_F[1][current_sipm], "CHI2/NDF"));
  // cout << Ref_Hist_F[multiplicity]->Chi2Test(TreeHist_F[multiplicity], "CHI2/NDF") << endl;

  double sum = accumulate(chi2_vec.begin(), chi2_vec.end(), 0.0);
  chi2 = sum / chi2_vec.size();

  double sq_sum = inner_product(chi2_vec.begin(), chi2_vec.end(), chi2_vec.begin(), 0.0);
  chi2_std = sqrt(sq_sum / chi2_vec.size() - chi2 * chi2);

  cout << chi2 << " +/- " << chi2_std << "    {";
  for (size_t i = 0; i < 7; ++i)
  {
    cout << bestPar[i] << ", ";
  }
  cout << "};";
  if (bestCHI2 > chi2)
  {
    bestCHI2 = chi2;
    cout << "    BEST" << endl;
  }
  else
  {
    cout << endl;
  }
  return chi2;
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
  // bestPar_sipm.push_back({6.06931, 4.96007, 0, 3.82433, 2.24847e-05, 81.1015, 19.4501, }); 
  // bestPar_sipm.push_back({2.58732, 5.0408, 0, 2.90152, 2.9204e-05, 87.252, 16.7987, });
  // bestPar_sipm.push_back({ 0, 0,0,0,0,0,0});
  // bestPar_sipm.push_back({2.7675, 4.99999, 0, 2.9588, 2.78769e-05, 83.4966, 17.721, });
  // bestPar_sipm.push_back({14.9584, 4.83285, 0, 2.11295, 3.7311e-05, 74.5177, 12.5693, });
  // bestPar_sipm.push_back({8.32008, 5.00363, 0, 4.99757, 3.31788e-05, 83.5593, 15.6704, });
  // bestPar_sipm.push_back({9.31388, 4.9, 0, 2.0, 2.21952e-05, 81.7, 18, });

  bestPar_sipm.push_back({ 0, 0,0,0,0,0,0});
  bestPar_sipm.push_back({0, 5.01, 0, 4.719, 1.89548e-05, 75.4259, 13.6421, });
  bestPar_sipm.push_back({0, 5.0, 0, 5.2043, 2.00847e-05, 75.7, 13.0109, }); 
  bestPar_sipm.push_back({0, 5.152, 0, 5.40899, 2.259526e-05, 80.4515, 13.9441, });
  bestPar_sipm.push_back({ 0, 0,0,0,0,0,0});
  bestPar_sipm.push_back({0, 5.05118, 0, 3.87573, 3.08645e-05, 79.2214, 18.3628, });
  bestPar_sipm.push_back({0, 4.93709, 0, 3.73354, 3.31023e-05, 60.7, 10.0289, });
  bestPar_sipm.push_back({0, 5.01339, 0, 3.43591, 3.57933e-05, 79.4809, 14.7005, });
  bestPar_sipm.push_back({0, 4.95, 0, 3.75, 2.01322e-05, 70.7679, 15.9544, });

  

  for (int sipm = 0; sipm <= 9; sipm++)
  {
    if (sipm == 4 || sipm == 0 || sipm == 9)
    {
      bestPar_sipm.push_back({ 0, 0,0,0,0,0,0});
      continue;
    }
    countercd++;
    cout << "SiPM : " << sipm << endl;
    current_sipm = sipm;

    ROOT::Math::Functor functor(&Chi2TreeHist, 7);
    minimizer->SetFunction(functor);
    minimizer->SetFixedVariable(0, "Calibration_OffSet", 0);
    minimizer->SetLimitedVariable(1, "Calibration", 5, 0.1, 4.8, 5.2);
    minimizer->SetFixedVariable(2, "Resolution_OffSet", 0);
    minimizer->SetLimitedVariable(3, "Resolution_SQRT", 3.5, 1, 2.5, 5.5);
    minimizer->SetLimitedVariable(4, "Resolution_2", 1e-5, 1e-7, 1e-7, 5e-5);
    minimizer->SetLimitedVariable(5, "Threshold", 81.7, 5, 60, 90);
    minimizer->SetLimitedVariable(6, "Threshold_STD", 18, 5, 10, 30);
    minimizer->SetPrecision(0.000001);
    minimizer->SetTolerance(0.000001);
    minimizer->SetMaxFunctionCalls(10000000);
    minimizer->SetMaxIterations(10000000);
    bestCHI2 = 1e6;

    // minimizer->Minimize();
    // const double *bestPar = minimizer->X();
    // minimizer->PrintResults();
    
    // bestPar_sipm.push_back({ bestPar[0], bestPar[1], bestPar[2], bestPar[3], bestPar[4], bestPar[5], bestPar[6]});

    const double *bestPar = bestPar_sipm[sipm].data();

    chi2_sipm[sipm] = Chi2TreeHist(bestPar);

    cF->cd(countercd);
    TreeHist_F[1][sipm]->GetXaxis()->SetRangeUser(0, 8000);
    Ref_Hist_F[1][sipm]->GetXaxis()->SetRangeUser(0, 8000);
    TreeHist_F[1][sipm]->Scale(Ref_Hist_F[1][sipm]->Integral() / TreeHist_F[1][sipm]->Integral());

    Ref_Hist_F[1][sipm]->GetXaxis()->SetRangeUser(0, 2000);
    TreeHist_F[1][sipm]->GetXaxis()->SetRangeUser(0, 2000);
    Ref_Hist_F[1][sipm]->Draw("HIST");
    TreeHist_F[1][sipm]->SetLineColor(kRed);
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
      value = (e-bestPar_sipm[sipm][0])/bestPar_sipm[sipm][1];
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
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MakeSiPM_MultiplicityPlots()
{
  double chi2 = 0;
  double chi2_std = 0;
  chi2_vec.clear();

  
  for (int multiplicity = 0; multiplicity <= 7; multiplicity++)
  {
    TreeHist_Multi_F[multiplicity] = (TH1D *)Ref_Hist_Multi_F[multiplicity]->Clone(("SiPM_MF" + to_string(multiplicity)).c_str());
    TreeHist_Multi_F[multiplicity]->Reset();
  }

  Reader_calib_sipm_F->Restart();
  TTreeReaderValue<double> Tree_SiPM(*Reader_calib_sipm_F, "PlasticScintillator_Deposit_Energy");

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
      if (i == 0) /////mean of sipm
      {
        counter_sum++;
        value_sum += (bestPar_sipm[*label][0] + bestPar_sipm[*label][1] * (*value));
        
        if (counter_sum == 7)
        {
          Ref_Hist_Multi_F[i]->Fill(value_sum / 7);
          value_sum = 0;
          counter_sum = 0;
        }
      }

      else ///multiplicity
      {
        Ref_Hist_Multi_F[i]->Fill(bestPar_sipm[*label][0] + bestPar_sipm[*label][1] * (*value));
      }
    }
  }
  /////////////////////////////////////////////////
  

  vector<double> vec_e;
  double energy_sum = 0;
  double sigma_resolution = 0;
  double sigma_resolution_SUM = 0;
  bool negative_res = false;

  for (int i = 0; i < 10; i++)
  {
    Reader_calib_sipm_F->Restart();
    while (Reader_calib_sipm_F->Next())
    {
      int real_multiplicity = 0;
      vec_e.clear();
      double energy_sum_final = 0;

      for (int sipm = 1; sipm <= 8; sipm++)
      {
        double energy = (*Tree_SiPM);
        energy_sum = (*Tree_SiPM);
        if (sipm == 4)
        {
          continue;
        }
        sigma_resolution = sqrt(pow(bestPar_sipm[sipm][2], 2) + pow(bestPar_sipm[sipm][3] * sqrt(energy), 2) + pow(bestPar_sipm[sipm][4] * pow(energy, 2), 2));
        normal_distribution<> resolution(0, sigma_resolution);

        sigma_resolution_SUM = sqrt(pow(bestPar_sipm[sipm][2], 2) + pow(bestPar_sipm[sipm][3] * sqrt(energy), 2));
        normal_distribution<> resolution_SUM(0, sigma_resolution_SUM);

        energy += resolution(gen);
        energy_sum += resolution_SUM(gen);

        normal_distribution<> threshold(bestPar_sipm[sipm][5], bestPar_sipm[sipm][6]);
        if (energy > threshold(gen))
        {
          real_multiplicity++;
          vec_e.push_back(energy);
          energy_sum_final += energy_sum;
        }
      }
      
      for (int multi = 1; multi <= real_multiplicity; multi++)
      {
        for (int sipm = 1; sipm <= real_multiplicity; sipm++)
        {
          TreeHist_Multi_F[multi]->Fill(vec_e[sipm - 1]);
        }
      }

      if (real_multiplicity == 7)
      {
        TreeHist_Multi_F[0]->Fill(energy_sum_final / 7);
      }
    }
  }

  int min = 250;
  int max = 2000;

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    Ref_Hist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 7000);
    TreeHist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 7000);
    chi2_vec.push_back(Ref_Hist_Multi_F[multiplicity]->Chi2Test(TreeHist_Multi_F[multiplicity], "CHI2/NDF"));
    cout << Ref_Hist_Multi_F[multiplicity]->Chi2Test(TreeHist_Multi_F[multiplicity], "CHI2/NDF") << endl;

  }

  double sum = accumulate(chi2_vec.begin(), chi2_vec.end(), 0.0);
  chi2 = sum / chi2_vec.size();

  double sq_sum = inner_product(chi2_vec.begin(), chi2_vec.end(), chi2_vec.begin(), 0.0);
  chi2_std = sqrt(sq_sum / chi2_vec.size() - chi2 * chi2);

  cout << chi2 << " +/- " << chi2_std << endl;;
  

  /////////////////PLOTTING
  TCanvas *cF = new TCanvas("SiPM_F", "SiPM_F", 1920, 1080);
  cF->Divide(3, 3);

  for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
  {
    cF->cd(multiplicity);
    TreeHist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 8000);
    Ref_Hist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 8000);
    TreeHist_Multi_F[multiplicity]->Scale(Ref_Hist_Multi_F[multiplicity]->Integral() / TreeHist_Multi_F[multiplicity]->Integral());

    Ref_Hist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
    TreeHist_Multi_F[multiplicity]->GetXaxis()->SetRangeUser(0, 2000);
    Ref_Hist_Multi_F[multiplicity]->Draw("HIST");
    TreeHist_Multi_F[multiplicity]->SetLineColor(kRed);
    TreeHist_Multi_F[multiplicity]->Draw("SAME HIST");

    TPaveText *pt = new TPaveText(0.7, 0.65, 1, 0.7, "brNDC");
    pt->SetFillColor(0);
    pt->AddText(("#chi^{2} = " + to_string(chi2_vec[multiplicity - 1])).c_str());
    pt->Draw("SAME");
  }

  Merged_File->cd();
  cF->Write();
}

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
      TreeHist_F[multiplicity][sipm] = (TH1D *)Ref_Hist_F[multiplicity][sipm]->Clone(("SiPM"+to_string(sipm)+"_F_M" + to_string(multiplicity)).c_str());
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
  

  vector<double> vec_e = {0,0,0,0,0,0,0,0,0,0};
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
    TCanvas *cF = new TCanvas(("SiPM"+to_string(sipm)+"_F").c_str(), ("SiPM"+to_string(sipm)+"_F").c_str(), 1920, 1080);
    cF->Divide(3, 3);

    for (int multiplicity = 1; multiplicity <= 7; multiplicity++)
    {
      cF->cd(multiplicity);
      
        TreeHist_F[1][sipm]->GetXaxis()->SetRangeUser(0, 8000);
        Ref_Hist_F[1][sipm]->GetXaxis()->SetRangeUser(0, 8000);
        TPaveText *pt = new TPaveText(0.7, 0.65, 1, 0.7, "brNDC");
      pt->SetFillColor(0);
      Ref_Hist_F[multiplicity][sipm]->GetXaxis()->SetRangeUser(0, 8000);
      TreeHist_F[multiplicity][sipm]->GetXaxis()->SetRangeUser(0, 8000);
      pt->AddText(("#chi^{2} = " + to_string(Ref_Hist_F[multiplicity][sipm]->Chi2Test(TreeHist_F[multiplicity][sipm], " CHI2/NDF"))).c_str());
        TreeHist_F[multiplicity][sipm]->Scale(Ref_Hist_F[multiplicity][sipm]->Integral() / TreeHist_F[multiplicity][sipm]->Integral());
      


      Ref_Hist_F[multiplicity][sipm]->GetXaxis()->SetRangeUser(0, 2000);
      TreeHist_F[multiplicity][sipm]->GetXaxis()->SetRangeUser(0, 2000);
      Ref_Hist_F[multiplicity][sipm]->Draw("HIST");
      TreeHist_F[multiplicity][sipm]->SetLineColor(kRed);
      TreeHist_F[multiplicity][sipm]->Draw("SAME HIST");

      
      pt->Draw("SAME");
    }

    Merged_File->cd();
    cF->Write();
  }
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
          {  // histF->Scale(histGT->Integral() / histF->Integral());

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
  TCanvas* cFGT = new TCanvas("cFGT", "cFGT", 800, 600);
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