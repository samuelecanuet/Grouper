#ifndef DETECTORS_HH
#define DETECTORS_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream> 
#include <gsl/gsl_statistics.h>
#include "TFile.h"

#include "/home/local1/Documents/lib/GTools1.0/include/GString.hh"


#define FDATA_MAX 80  ///< Maximum coder label
#define SIGNAL_MAX 80 ///< Maximum number of signals
#define BETA_NUM 2    ///< Number of beta (SiPM) detector groups
#define SILI_NUM 8    ///< Number of silicon detector groups
#define BETA_SIZE 9   ///< Number of channels in beta (SiPM) detector groups
#define SILI_SIZE 6   ///< Number of channels in silicon detector groups

#define BETA_HI 1
#define BETA_LO 2

using namespace std;

string detectorFileName; ///< Detectors definition file name
size_t detectorNum;      ///< Number of defined detectors channels
string detectorName[SIGNAL_MAX];
int detectorCoder[SIGNAL_MAX]; ///< Coder label of all signals
int detectorInfo[SIGNAL_MAX];  ///< Detector information of all signals (type identifier)

int detBeta[BETA_NUM][BETA_SIZE]; ///< Detector number of beta signals
int detSili[SILI_NUM][SILI_SIZE]; ///< Detector number of silicon signals

int coderDetector[FDATA_MAX]; ///< Detector number of all coders


/// Silicon ///
double winSiliMin = -50;
double winSiliMax = 50;
int winSiliN = (abs(winSiliMin) + winSiliMax) / 8;
double eSiliMin = 0;
double eSiliMax = 100000;
int eSiliN = 10000;

/// SiPM ///
/// High
double winHighMin = 100;
double winHighMax = 250;
int winHighN = (abs(winHighMin) + winHighMax) / 2;
double eHighMin = 0;
double eHighMax = 4000000;
int eHighN = 4000;

/// Low
double winLowMin = 100;
double winLowMax = 250;
int winLowN = (abs(winLowMin) + winLowMax) / 2;
double eLowMin = 0;
double eLowMax = 4000000;
int eLowN = 4000;

/// TOTAL WINDOW ///
double tabMIN[3] = {winSiliMin, winHighMin, winLowMin};
double winTotalMin = gsl_stats_min(tabMIN, 1, 3);
double tabMAX[3] = {winSiliMax, winHighMax, winLowMax};
double winTotalMax = gsl_stats_max(tabMAX, 1, 3);


inline bool IsDetectorBeta(int det)
{
  bool res = false;
  if (det >= 0)
    res = ((detectorInfo[det] / 100) == 0);
  return (res);
}

inline bool IsDetectorBetaLow(int det)
{
  bool res = false;
  res = ((detectorInfo[det] / 20) == 1);
  return (res);
}

inline bool IsDetectorBetaHigh(int det)
{
  bool res = false;
  res = (!IsDetectorBetaLow(det) && IsDetectorBeta(det));
  return (res);
}

inline bool IsDetectorSili(int det)
{
  bool res = false;
  if (det >= 0)
    res = ((detectorInfo[det] / 100) == 1);
  return (res);
}

inline int GetDetectorChannel(int det)
{
  return (detectorInfo[det] % 10);
}

inline int GetDetector(int det)
{
  return (detectorInfo[det] / 10) % 10;
}

//----------------------------------------------------------------------
inline bool IsDetectorSiliBack(int det)
{
  return (IsDetectorSili(det) && (GetDetectorChannel(det) == 6));
}

inline bool IsDetectorSiliStrip(int det)
{
  return (IsDetectorSili(det) && (GetDetectorChannel(det) != 6));
}

inline bool IsSameSiliDetector(int detstrip1, int detstrip2)
{
  bool res = false;
    res = (GetDetector(detstrip1) == GetDetector(detstrip2));
  return res;
}

inline int InitDetectors(const string &fname)
{
  int error = 0;

  FILE *fp = fopen(fname.c_str(), "r");
  if (fp != NULL)
  {
    detectorNum = 0;
    while (!feof(fp))
    {
      char fline[512];
      fgets(fline, 511, fp);

      if (!feof(fp))
      {
        string defline = GString(fline).NoEndSpace();
        if (defline.length() > 0)
        {
          if (defline.at(0) != '#') // not a comment
          {
            istringstream iss(defline);
            int coder;
            string name;
            iss >> coder >> name;

            if (!iss.fail())
            {
              detectorCoder[detectorNum] = coder;
              detectorName[detectorNum] = name;

              // set detector groups
              if (name.substr(0, 4) == "Beta")
              {
                detectorInfo[detectorNum] = 0;
                int d = -1;
                sscanf(&name.c_str()[6], "%d", &d);
                if ((d > 0) && (d <= 9))
                {
                  if (name.substr(4, 2) == "Hi")
                  {
                    detBeta[0][d - 1] = detectorNum;
                    detectorInfo[detectorNum] += 10 + d;
                    // cerr << "** BetaHi " << d << " = " << detectorInfo[detectorNum] << endl;
                  }
                  else if (name.substr(4, 2) == "Lo")
                  {
                    detBeta[1][d - 1] = detectorNum;
                    detectorInfo[detectorNum] += 20 + d;
                    // cerr << "** BetaLo " << d << " = " << detectorInfo[detectorNum] << endl;
                  }
                  else
                    GLogWarning("Bad beta detector: " + defline);
                }
                else
                  GLogWarning("Bad beta detector number: " + defline);
              }
              else if (name.substr(0, 2) == "Si")
              {
                detectorInfo[detectorNum] = 100;
                name.at(3) = ' ';
                if (name.at(4) == 'R')
                  name.at(4) = '6';
                int s = -1;
                int d = -1;
                sscanf(&name.c_str()[2], "%d %d", &s, &d);
                if ((s > 0) && (s <= 8) && (d > 0) && (d <= 6))
                {
                  detSili[s - 1][d - 1] = detectorNum;
                  detectorInfo[detectorNum] += 10 * s + d;
                }
                else
                  GLogWarning("Bad silicon detector number: " + defline);
              }
              else
                GLogWarning("Error in definition line: " + defline);

              detectorNum++;
            }
            else
              GLogWarning("Error in definition line: " + defline);
          }
        }
      }
    } // read loop

    fclose(fp);
    GLogMessage("<WIS> Number of defined detectors: " + GGetString((int)detectorNum));

    // update coder -> detector conversion table
    for (size_t i = 0; i < detectorNum; ++i)
    {
      // cerr << detectorName[i] << " -> " << GGetString(detectorInfo[i],3) << endl;
      if (detectorCoder[i] >= 0)
        coderDetector[detectorCoder[i]] = i;
    }
  }
  else
    GLogWarning("Error opening definition file: " + fname);

  return (error);
}

void WriteTime(TFile* from_File, TFile* to_file)
{
  TObject* start = from_File->Get("Start_time");
  TObject* stop = from_File->Get("Stop_time");

  to_file->cd();
  TNamed("Start_time", start->GetTitle()).Write();
  TNamed("Stop_time", stop->GetTitle()).Write();
}

pair<string, string> GetTime(TFile* File)
{
  TObject* start = File->Get("Start_time");
  TObject* stop = File->Get("Stop_time");

  string str_start = start->GetTitle();
  string str_stop = stop->GetTitle(); 

  return make_pair(str_start, str_stop);
}

#endif