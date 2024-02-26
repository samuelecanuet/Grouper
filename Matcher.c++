#include "Matcher.hh"
#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

using namespace std;
using namespace ROOT::Math;

int main(int argc, char *argv[])
{
    if (argc == 1)
    {
        cout << "Default Runs Used : 034" << endl;
        RunCal = 34;
    }
    else
    {
        RunCal = atoi(argv[1]);
    }

    ///////////INITIALISATION///////////
    TFile *OutputFile = new TFile(("./Matched/run_0"+to_string(RunCal)+"_32Ar_matched.root").c_str(), "RECREATE");
    TFile *Ref1234_File = new TFile(("./Cleaned/run_0"+to_string(RunRef_1234)+"_32Ar_cleaned.root").c_str(), "READ");
    TFile *Ref5678_File = new TFile(("./Cleaned/run_0"+to_string(RunRef_5678)+"_32Ar_cleaned.root").c_str(), "READ");
    TFile *Cal_File = new TFile(("./Cleaned/run_0"+to_string(RunCal)+"_32Ar_cleaned.root").c_str(), "READ");

    TGraphErrors* Result = new TGraphErrors();
    Result->SetDrawOption("*AP");
    Result->SetTitle("Proportionality;Detector;Proportionality");
    int result_counter = 0;
    InitDetectors("../Analysis/ReadFaster/Data/detectors.dat");

    ///////////MINIMIZER///////////
    string dir;
    string det;
    int rebin;
    for (size_t i = 0; i < detectorNum; ++i)
    {
        if (IsDetectorSili(i))
        {
            det = "Silicon";
            range_low = 20000;
            range_high = 90000;
            rebin = 1;
            if (GetDetector(i) <= 4 )
            {
                Ref_File = Ref1234_File;
            }
            else
            {
                Ref_File = Ref5678_File;
            }

            if (IsDetectorSiliStrip(i))
            {
                
                dir = "Strip";
            }
            else
            {
                dir = "Rear";
            }
        }
        else if (IsDetectorBeta(i))
        {
            det = "SiPM";
            Ref_File = Ref5678_File;
            range_low = 0;
            range_high = 4000000;
            rebin = 1;
            guess_high = 1.2;
            guess_low = 0.8;
            if (IsDetectorBetaHigh(i))
            {
                dir = "SiPMHigh";
            }
            else
            {
                dir = "SiPMLow";
                rebin = 8;
            }
        }
            Ref_Hist = (TH1I *)Ref_File->Get((dir+"_Channel"+"/H"+dir+"_Channel_C_"+detectorName[i]).c_str());
            Ref_Hist->Rebin(rebin);
            Tree = (TTree *)Cal_File->Get(("Tree_"+det+"_"+detectorName[i]).c_str());
            Reader = new TTreeReader (Tree);

            graph = new TGraph();
            graph->SetTitle("Gain_Matching;Proportionality; #chi^{2}");
            counter=0;
            
            //////////////////// MINIMIZER ////////////////////
            ////// First
            Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
            ROOT::Math::Functor functor(&Chi2TreeHist, 1);
            minimizer->SetFunction(functor);
            minimizer->SetLimitedVariable(0, "Proportionality", 1., .05, guess_low, guess_high);
            minimizer->SetPrecision(0.05);
            minimizer->Minimize();
            const double *bestPar = minimizer->X();
            
            ////// Second
            minimizer->SetLimitedVariable(0, "Proportionality", bestPar[0], .05, guess_low, guess_high);
            minimizer->SetPrecision(1e-5);
            minimizer->Minimize();
            bestPar = minimizer->X();

            ////// Third
            minimizer->SetLimitedVariable(0, "Proportionality", bestPar[0], .05, guess_low, guess_high);
            minimizer->SetPrecision(1e-10);
            minimizer->Minimize();
            bestPar = minimizer->X();
            double error = minimizer->Errors()[0];
            std::cout << detectorName[i] << "   CHI2 : "<< chi2 <<"    Best Proportionality: " << bestPar[0] << " +/- " << error << std::endl;



            //////////////////// SAVE ////////////////////
            TCanvas *c = new TCanvas(detectorName[i].c_str(), detectorName[i].c_str(), 800, 600);
            c->Divide(1, 2);
            c->cd(1);
            Ref_Hist->Draw("HIST");

            Reader->Restart();
            TTreeReaderValue<int> Channel(*Reader, "Channel");

            TH1I *TreeHist = (TH1I *)Ref_Hist->Clone((detectorName[i] + "_" + to_string(bestPar[0])).c_str());
            TreeHist->Reset();
            

            while (Reader->Next())
            {
                TreeHist->Fill((*Channel) * bestPar[0]);
            }

            TreeHist->SetLineColor(kRed);
            TreeHist->Draw("SAME");

            c->cd(2);
            graph->Draw("*AP");
            OutputFile->cd();
            c->Write();

            result_counter++;
            Result->SetPoint(result_counter, i, bestPar[0]);
            Result->SetPointError(result_counter, 0, error);

            
        }

    Result->SetName("10000");
    Result->Write();

    return 0;
}
