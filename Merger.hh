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
string baseFileName;
double detectorCalib[SIGNAL_MAX];
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

TH2D *HSiliconRun_Channel_Unmatched[SIGNAL_MAX];
TH2D *HSiliconRun_Channel_Matched[SIGNAL_MAX];
TH2D *HSiliconRun[SIGNAL_MAX];

TH1D *HSilicon_Channel_Unmatched_temp[SIGNAL_MAX];
TH1D *HSilicon_Channel_Matched_temp[SIGNAL_MAX];
TH1D *HSilicon_temp[SIGNAL_MAX];
TH1D *to_save_unmatched;
TH1D *to_save_matched;




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

vector<string> listFilesInDirectory(const string &directoryPath)
{
  vector<string> files;
  fs::path path(directoryPath);
  if (fs::exists(path) && fs::is_directory(path))
  {
    for (auto &entry : fs::directory_iterator(path))
    {
      if (fs::is_regular_file(entry))
      {
        files.push_back(entry.path().filename().string());
      }
    }
  }
  return files;
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
  double t_max = file_Time[file_Time.size()-1].second;
  n_bin_time = (int)(t_max - t_min) / bin_time;

  cout<<t_min<<" "<<t_max<<" "<<n_bin_time<<endl;
  //////////////////////////////////////////

  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      HSiliconRun_Channel_Unmatched[i] = new TH2D(("HSiliconRun_Channel_Unmatched_" + detectorName[i]).c_str(), ("HSiliconRun_Channel_Unmatched_" + detectorName[i]).c_str(), n_bin_time, t_min, t_max,  eSiliN, eSiliMin, eSiliMax);
      HSiliconRun_Channel_Unmatched[i]->GetXaxis()->SetTitle("Time (h)");
      HSiliconRun_Channel_Unmatched[i]->GetYaxis()->SetTitle("Channel");
      HSiliconRun_Channel_Unmatched[i]->GetXaxis()->CenterTitle();
      HSiliconRun_Channel_Unmatched[i]->GetYaxis()->CenterTitle();

      HSiliconRun_Channel_Matched[i] = new TH2D(("HSiliconRun_Channel_Matched_" + detectorName[i]).c_str(), ("HSiliconRun_Channel_Matched_" + detectorName[i]).c_str(), n_bin_time, t_min, t_max,  eSiliN, eSiliMin, eSiliMax);
      HSiliconRun_Channel_Matched[i]->GetXaxis()->SetTitle("Runs");
      HSiliconRun_Channel_Matched[i]->GetYaxis()->SetTitle("Channel");
      HSiliconRun_Channel_Matched[i]->GetXaxis()->CenterTitle();
      HSiliconRun_Channel_Matched[i]->GetYaxis()->CenterTitle();

      HSiliconRun[i] = new TH2D(("HSiliconRun"+ detectorName[i]).c_str(), ("HSiliconRun"+ detectorName[i]).c_str(), n_bin_time, t_min, t_max, 10000, 0, 10000);
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

      HSilicon[i] = new TH1D(("HStrip_" + detectorName[i]).c_str(), ("HStrip_" + detectorName[i]).c_str(), 10000, 0, 10000);
      HSilicon[i]->GetXaxis()->SetTitle("Energy [keV]");
      HSilicon[i]->GetYaxis()->SetTitle("Counts");
      HSilicon[i]->GetXaxis()->CenterTitle();
      HSilicon[i]->GetYaxis()->CenterTitle();

      HSilicon_coinc[i] = new TH1D(("HStrip_coinc_" + detectorName[i]).c_str(), ("HStrip_coinc_" + detectorName[i]).c_str(), 10000, 0, 10000);
      HSilicon_coinc[i]->GetXaxis()->SetTitle("Energy [keV]");
      HSilicon_coinc[i]->GetYaxis()->SetTitle("Counts");
      HSilicon_coinc[i]->GetXaxis()->CenterTitle();
      HSilicon_coinc[i]->GetYaxis()->CenterTitle();

      HSilicon_no_coinc[i] = new TH1D(("HStrip_no_coinc_" + detectorName[i]).c_str(), ("HStrip_no_coinc_" + detectorName[i]).c_str(), 10000, 0, 10000);
      HSilicon_no_coinc[i]->GetXaxis()->SetTitle("Energy [keV]");
      HSilicon_no_coinc[i]->GetYaxis()->SetTitle("Counts");
      HSilicon_no_coinc[i]->GetXaxis()->CenterTitle();
      HSilicon_no_coinc[i]->GetYaxis()->CenterTitle();
    }

    if (IsDetectorBetaHigh(i))
    {
      HSiPM[GetDetectorChannel(i)] = new TH1D(("HSiPM_M"+to_string(GetDetectorChannel(i))).c_str(), ("HSiPM_M"+to_string(GetDetectorChannel(i))).c_str(), eLowN/10, eLowMin, eLowMax*5);
      HSiPM[GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
      HSiPM[GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
      HSiPM[GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
      HSiPM[GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

      HSiPM_F[GetDetectorChannel(i)] = new TH1D(("HSiPM_F_M"+to_string(GetDetectorChannel(i))).c_str(), ("HSiPM_F_M"+to_string(GetDetectorChannel(i))).c_str(), eLowN/10, eLowMin, eLowMax*5);
      HSiPM_F[GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
      HSiPM_F[GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
      HSiPM_F[GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
      HSiPM_F[GetDetectorChannel(i)]->GetYaxis()->CenterTitle();

      HSiPM_GT[GetDetectorChannel(i)] = new TH1D(("HSiPM_GT_M"+to_string(GetDetectorChannel(i))).c_str(), ("HSiPM_GT_M"+to_string(GetDetectorChannel(i))).c_str(), eLowN/10, eLowMin, eLowMax*5);
      HSiPM_GT[GetDetectorChannel(i)]->GetXaxis()->SetTitle("Channel");
      HSiPM_GT[GetDetectorChannel(i)]->GetYaxis()->SetTitle("Counts");
      HSiPM_GT[GetDetectorChannel(i)]->GetXaxis()->CenterTitle();
      HSiPM_GT[GetDetectorChannel(i)]->GetYaxis()->CenterTitle();
    }
  }
  
}
void WriteRunHistograms(TFile* file, int i, double unmatched, double matched)
{
  pair<string, string> string_time = GetTime(file);
  start_time = Convert_DatetoTime(string_time.first, file_string_Time[0].first);
  stop_time = Convert_DatetoTime(string_time.second, file_string_Time[0].first);

  for (double j = start_time; j < stop_time; j += bin_time)
  {
    HSiliconRun_Channel_Unmatched[i]->Fill(j, unmatched);
    HSiliconRun_Channel_Matched[i]->Fill(j, matched);
  }
}

void WriteHistograms()
{
  Merged_File->cd();
  TDirectory *dir_HSilicon_Channel_Unmatched = Merged_File->mkdir("Silicon_Unmatched");
  TDirectory *dir_HSilicon_Channel_Matched = Merged_File->mkdir("Silicon_Matched");
  TDirectory *dir_HSilicon = Merged_File->mkdir("Silicon_Calibrated");
  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorSiliStrip(i))
    {
      dir_HSilicon_Channel_Unmatched->cd();
      HSilicon_Channel_Unmatched[i]->Write();
      HSiliconRun_Channel_Unmatched[i]->Write();

      dir_HSilicon_Channel_Matched->cd();
      HSilicon_Channel_Matched[i]->Write();
      HSiliconRun_Channel_Matched[i]->Write();
      
      dir_HSilicon->cd();
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
  for (int i = 0; i < BETA_SIZE + 1; i++)
  {
    HSiPM[GetDetectorChannel(i)]->Write();
    HSiPM_F[GetDetectorChannel(i)]->Write();
    HSiPM_GT[GetDetectorChannel(i)]->Write();
  }
}

void InitFiles()
{
  vector<string> Grouped_Files = listFilesInDirectory("./Grouped/");
  vector<string> Matched_Files = listFilesInDirectory("./Matched/");

  auto processFiles = [&](const vector<string> &files)
  {
    for (const auto &file : files)
    {
      boost::filesystem::path filePath(file);
      string base_name = filePath.stem().string();

      string suffix1 = "_matched";
      string suffix2 = "_grouped";
      if (base_name.rfind(suffix1) != string::npos)
      {
        base_name = base_name.substr(0, base_name.size() - suffix1.size());
      }
      else if (base_name.rfind(suffix2) != string::npos)
      {
        base_name = base_name.substr(0, base_name.size() - suffix2.size());
      }

      if (fileMap.find(base_name) == fileMap.end())
      {
        fileMap[base_name] = vector<string>();
      }

      fileMap[base_name].push_back(file);
    }
  };

  processFiles(Grouped_Files);
  processFiles(Matched_Files);

  for (auto it = fileMap.begin(); it != fileMap.end();)
  {
    if (it->second.size() != 2)
    {
      it = fileMap.erase(it);
    }
    else
    {
      ++it;
    }
  }

  // Print the map
  // for (const auto &pair : fileMap)
  // {
  //   cout << pair.first << ":\n";
  //   for (const auto &file : pair.second)
  //   {
  //     cout << "  " << file << "\n";
  //   }
  // }
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
      HSilicon_Channel_Matched[i]->GetXaxis()->SetRangeUser(40000, 42000);
      detectorCalib[i] /= HSilicon_Channel_Matched[i]->GetMean();
    }
    
  }
}

int InitCalib()
{
  int error = 0;
  string CalibDir = "../Calib/";
  string CalibFileName = "Sim_32Ar_MylarAl_F.txt";
  ifstream file(CalibDir + CalibFileName);

  for (auto &value : detectorCalib)
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
      cout<<DetName<<endl;
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
            detectorCalib[i] = E;
          }
        }
      }
      else if (DetName.find("Down") != string::npos)
      {
        for (int i = 0; i < SIGNAL_MAX; i++)
        {
          if (GetDetector(i) > 4 and GetDetectorChannel(i) == strip)
          {
            detectorCalib[i] = E;
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

  for (int i = 0; i < SIGNAL_MAX; i++)
  {
    if (IsDetectorBetaHigh(i))
    {
      detectorCalib[i] = 1390./24000000;
    }
    if (IsDetectorBetaLow(i))
    {
      detectorCalib[i] = 6000./9000000;
    }
  }

  // FOR CHECK
  // int index = 0;
  // for (auto &name : detectorName)
  // {
  //   cout << name << "       " << detectorCalib[index] << endl;
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


