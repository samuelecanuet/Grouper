#ifndef ANALSYER_HH
#define ANALSYER_HH

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
#include "TProfile.h"
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


string Catcher;
vector<double> Convolution_Parameters;
TFile* Simulation_File[3];
double a[3];

vector<double> Extract_Parameters(TFile* Merged_File)
{
    vector<double> Parameters;
    return Parameters;
}

void FillHistograms_EXP()
{

}

void FillHistograms_SIM()
{

}

#endif