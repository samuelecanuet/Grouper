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
    if (RunCal <= 40)
    {
        RunRef_1234 = 16;
        RunRef_5678 = 30;
    }
    else
    {
        RunRef_1234 = 00;
        RunRef_5678 = 00;
    }
    cout<<"<SAM> Current Run : "<<RunCal<<endl;
    TFile *OutputFile = new TFile(("./Matched/run_0" + to_string(RunCal) + "_32Ar_matched.root").c_str(), "RECREATE");
    TFile *Ref1234_File = new TFile(("./Cleaned/run_0" + to_string(RunRef_1234) + "_32Ar_cleaned.root").c_str(), "READ");
    TFile *Ref5678_File = new TFile(("./Cleaned/run_0" + to_string(RunRef_5678) + "_32Ar_cleaned.root").c_str(), "READ");
    TFile *Cal_File = new TFile(("./Cleaned/run_0" + to_string(RunCal) + "_32Ar_cleaned.root").c_str(), "READ");


    //////////////MATCHING SILICON////////////
    cout<<"<SAM> Matching Silicons"<<endl;
    TDirectory *dir_Chi2Silicon = OutputFile->mkdir("Silicon_RunMatching");
    dir_Chi2Silicon->cd();
    TGraphErrors *Result_Si = new TGraphErrors();
    Result_Si->SetDrawOption("*AP");
    Result_Si->SetName("Silicon_RunMatching");
    Result_Si->SetTitle("Silicon_RunMatching;Detector Label;Silicon_RunMatching");
    int result_counter = 0;
    InitDetectors("../Analysis/ReadFaster/Data/detectors.dat");

    ///////////MINIMIZER///////////
    string dir;
    string det;
    for (size_t i = 0; i < detectorNum; ++i)
    {
        if (IsDetectorSili(i))
        {
            det = "Silicon";
            range_low = 20000;
            range_high = 90000;
            if (GetDetector(i) <= 4)
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

            Ref_Hist = (TH1D *)Ref_File->Get((dir + "_Channel" + "/H" + dir + "_Channel_C_" + detectorName[i]).c_str());
            Tree = (TTree *)Cal_File->Get(("Tree_" + det + "_" + detectorName[i]).c_str());
            Reader = new TTreeReader(Tree);

            graph = new TGraph();
            graph->SetTitle("Gain_Matching;Proportionality; #chi^{2}");
            counter = 0;

            //////////////////// MINIMIZER ////////////////////
            ////// First
            Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
            ROOT::Math::Functor functor(&Chi2TreeHist, 1);
            minimizer->SetFunction(functor);
            minimizer->SetLimitedVariable(0, "Proportionality", 1., .1, 0.8, 1.2);
            minimizer->SetPrecision(0.05);
            minimizer->Minimize();
            const double *bestPar = minimizer->X();

            ////// Second
            minimizer->SetLimitedVariable(0, "Proportionality", bestPar[0], .05, 0.8, 1.2);
            minimizer->SetPrecision(1e-5);
            minimizer->Minimize();
            bestPar = minimizer->X();

            ////// Third
            minimizer->SetLimitedVariable(0, "Proportionality", bestPar[0], .005, 0.8, 1.2);
            minimizer->SetPrecision(1e-10);
            minimizer->Minimize();
            bestPar = minimizer->X();
            double error = minimizer->Errors()[0];
            std::cout << detectorName[i] << "\t CHI2 : " << chi2 << "\t\t Best Proportionality : " << bestPar[0] << " +/- " << error << std::endl;

            //////////////////// SAVE ////////////////////
            TCanvas *c = new TCanvas(detectorName[i].c_str(), detectorName[i].c_str(), 800, 600);
            c->Divide(1, 2);
            c->cd(1);
            Ref_Hist->Draw("HIST");

            Reader->Restart();
            TTreeReaderValue<double> Channel(*Reader, "Channel");

            TH1D *TreeHist = (TH1D *)Ref_Hist->Clone((detectorName[i] + "_" + to_string(bestPar[0])).c_str());
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
            dir_Chi2Silicon->cd();
            c->Write();

            result_counter++;
            Result_Si->SetPoint(result_counter, i, bestPar[0]);
            SiliconMatching[i] = bestPar[0];
            Result_Si->SetPointError(result_counter, 0, error);
        }
    }
    dir_Chi2Silicon->cd();
    Result_Si->Write();

    /////// Coef between Low and High /////
    cout<<"<SAM> Fitting gain coefficients Low/High"<<endl;
    dir_Chi2SiPMLowHigh = OutputFile->mkdir("SiPM_LowHigh_Matching");
    dir_Chi2SiPMLowHigh->cd();
    TGraphErrors *Result_SiPM_LowHigh = new TGraphErrors();
    Result_SiPM_LowHigh->SetDrawOption("*AP");
    Result_SiPM_LowHigh->SetName("SiPM_LowHigh_Matching");
    Result_SiPM_LowHigh->SetTitle("SiPM_LowHigh_Matching;Detector Label;SiPM_LowHigh_Matching");
    result_counter=0;
    for (int i = 1; i <= BETA_SIZE; i++)
    {
        TCanvas *canvas = (TCanvas *)Cal_File->Get(("SiPM_Multiplicities/SiPM_" + to_string(GetDetectorChannel(i))).c_str());
        TGraph *graph = (TGraph *)canvas->GetPrimitive(("SiPM_" + to_string(GetDetectorChannel(i))).c_str());

        double par = 0;
        double err = 0;
        if (graph->GetN() != 0)
        {
            TF1 *f1 = new TF1("f1", "[0]*x", 60000, 180000);
            f1->SetParameter(0, 10);
            graph->Fit(f1, "QR");

            par = f1->GetParameter(0);
            err = f1->GetParError(0);

            OutputFile->cd();
            TCanvas *c = new TCanvas(("SiPM_" + to_string(GetDetectorChannel(i))).c_str(), ("SiPM_" + to_string(GetDetectorChannel(i))).c_str(), 800, 600);
            c->cd();
            graph->Draw("*AP");
            f1->Draw("SAME");
            dir_Chi2Silicon->cd();
            c->Write();
            delete c;
        }
        else
        {
            par = 0;
            err = 0;

            OutputFile->cd();
            TCanvas *c = new TCanvas(("SiPM_" + to_string(GetDetectorChannel(i))).c_str(), ("SiPM_" + to_string(GetDetectorChannel(i))).c_str(), 800, 600);
            c->cd();
            dir_Chi2SiPMLowHigh->cd();
            c->Write();
            delete c;
        }

        

        result_counter++;
        Result_SiPM_LowHigh->SetPoint(result_counter, i, par);
        SiPMLowHighMatching[i] = par;
        Result_SiPM_LowHigh->SetPointError(result_counter, 0, err);
    }

    dir_Chi2SiPMLowHigh->cd();
    Result_SiPM_LowHigh->Write();

    OutputFile->cd();

    
    /////// Merging Low and High, apply matching on Silicon ///////
    cout<<"<SAM> Merging Low/High"<<endl;
    InitHistograms();
    Tree_Read = (TTree *)Cal_File->Get("Tree");
    Reader = new TTreeReader(Tree_Read);
    TTreeReaderArray<Signal> *Tree_Silicon = new TTreeReaderArray<Signal>(*Reader, "Tree_Silicon");
    TTreeReaderArray<Signal> *Tree_SiPMHigh = new TTreeReaderArray<Signal>(*Reader, "Tree_SiPMHigh");
    TTreeReaderArray<Signal> *Tree_SiPMLow = new TTreeReaderArray<Signal>(*Reader, "Tree_SiPMLow");

    Tree_Write = new TTree("Tree_LowHighMatched", "Tree_LowHighMatched");
    vector<Signal> Tree_Silicon_W;
    vector<Signal> Tree_SiPM_W;
    Tree_Write->Branch("Tree_Silicon", &Tree_Silicon_W);
    Tree_Write->Branch("Tree_SiPM", &Tree_SiPM_W);

    while (Reader->Next())
    {
        Tree_Silicon_W.clear();
        Tree_SiPM_W.clear();
        //////MATCHING SILICON
        for (auto &signal : *Tree_Silicon)
        {
            signal.Channel = SiliconMatching[signal.Label] * signal.Channel;
            Tree_Silicon_W.push_back(signal);
        }
        

        //////MATCHING AND MERGING SiPM Low and High
        bool doHigh = true;

        for (auto &signal : *Tree_SiPMHigh)
        {
            if (1 / SiPMLowHighMatching[GetDetectorChannel(signal.Label)] * signal.Channel > 180000)
            {
                doHigh = false;
            }
        }

        if (doHigh)
        {
            for (auto &signal : *Tree_SiPMHigh)
            {
                signal.Channel = 1 / SiPMLowHighMatching[GetDetectorChannel(signal.Label)] * signal.Channel;
                Tree_SiPM_W.push_back(signal);

                HSiPMHigh_Channel[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
                HSiPM_Channel[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
            }
        }
        else
        {
            for (auto &signal : *Tree_SiPMLow)
            {
                Tree_SiPM_W.push_back(signal);
                HSiPMLow_Channel[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
                HSiPM_Channel[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
            }
        }
        Tree_Write->Fill();

        /////FOR SAVING BEFORE MERGING//////
        for (auto &signal : *Tree_SiPMHigh)
        {
            HSiPMHigh_Channel_all[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
        }

        for (auto &signal : *Tree_SiPMLow)
        {
            HSiPMLow_Channel_all[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
        }
        ///////////////////////////////////
    }

    Tree_Write->Write();
    

    ////// Matching all SiPMs/////
    cout<<"<SAM> Matching SiPMs"<<endl;
    dir_Chi2SiPMs = OutputFile->mkdir("SiPMs_Matching");
    dir_Chi2SiPMs->cd();
    TGraphErrors *Result_SiPMs = new TGraphErrors();
    Result_SiPMs->SetDrawOption("*AP");
    Result_SiPMs->SetName("SiPMs_Matching");
    Result_SiPMs->SetTitle("SiPMs_Matching;Detector Label;SiPMs_Matching");
    result_counter=0;
    Ref_Hist = HSiPM_Channel[1];  //////////////////ref is 1
    range_low = 400000;
    range_high = 1000000;
    Tree = (TTree *)OutputFile->Get("Tree_LowHighMatched");
    Reader = new TTreeReader(Tree_Write);
    for (size_t i = 1; i <= BETA_SIZE; ++i)
    {
        graph = new TGraph();
        graph->SetTitle("Gain_Matching;Proportionality; #chi^{2}");
        counter = 0;

        current_SiPM = i;

        //////////////////// MINIMIZER ////////////////////
            ////// First
            Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
            ROOT::Math::Functor functor(&Chi2SiPM, 1);
            minimizer->SetFunction(functor);
            minimizer->SetLimitedVariable(0, "Proportionality", 1., .05, 0.6, 1.2);
            minimizer->SetPrecision(0.05);
            minimizer->Minimize();
            const double *bestPar = minimizer->X();

            ////// Second
            minimizer->SetLimitedVariable(0, "Proportionality", bestPar[0], .05, 0.6, 1.2);
            minimizer->SetPrecision(1e-5);
            minimizer->Minimize();
            bestPar = minimizer->X();

            ////// Third
            minimizer->SetLimitedVariable(0, "Proportionality", bestPar[0], .005, 0.6, 1.2);
            minimizer->SetPrecision(1e-10);
            minimizer->Minimize();
            bestPar = minimizer->X();
            double error = minimizer->Errors()[0];
            std::cout << "SiPM "+to_string(i) << "\t CHI2 : " << chi2 << "\t\t Best Proportionality : " << bestPar[0] << " +/- " << error << std::endl;

            //////////////////// SAVE ////////////////////
            TCanvas *c = new TCanvas(("HSiPM_Matched" + to_string(i) + "_Channel").c_str(), ("HSiPM_Matched" + to_string(i) + "_Channel").c_str(), 800, 600);
            c->Divide(1, 3);
            c->cd(1);
            Ref_Hist->Draw("HIST");

            Reader->Restart();
            TTreeReaderArray<Signal> *SiPM = new TTreeReaderArray<Signal>(*Reader, "Tree_SiPM");
            TH1D *TreeHist = (TH1D *)Ref_Hist->Clone((detectorName[i] + "_" + to_string(bestPar[0])).c_str());
            TreeHist->Reset();
            while (Reader->Next())
            {
                for (size_t j = 0; j < SiPM->GetSize(); j++)
                {
                    if (current_SiPM == GetDetectorChannel((*SiPM)[j].Label))
                        TreeHist->Fill((*SiPM)[j].Channel * bestPar[0]);
                }
            }

            TreeHist->SetLineColor(kRed);
            TreeHist->GetXaxis()->UnZoom();
            TreeHist->Draw("SAME");

            c->cd(2);
            TH1D *h = (TH1D *)Ref_Hist->Clone("Diff");
            h->Add(TreeHist, -1.0);
            h->Draw("EP");


            c->cd(3);
            graph->Draw("*AP");
            OutputFile->cd();
            dir_Chi2SiPMs->cd();
            c->Write();

            result_counter++;
            Result_SiPMs->SetPoint(result_counter, i, bestPar[0]);
            SiPMsMatching[i] = bestPar[0];
            Result_SiPMs->SetPointError(result_counter, 0, error);
    }

    dir_Chi2SiPMs->cd();
    Result_SiPMs->Write();


    ///// Apply Matching between SiPMs and Merge them ///////
    cout<<"<SAM> Merging SiPMs"<<endl;
    Reader->Restart();
    Tree_Silicon = new TTreeReaderArray<Signal>(*Reader, "Tree_Silicon");
    TTreeReaderArray<Signal> *Tree_SiPM = new TTreeReaderArray<Signal>(*Reader, "Tree_SiPM");

    Tree_Final = new TTree("Tree", "Tree");
    vector<Signal> Tree_SiPM_final = {Signal()};
    vector<Signal> Tree_Silicon_final;
    Tree_Final->Branch("Tree_Silicon", &Tree_Silicon_final);
    Tree_Final->Branch("Tree_SiPM", &Tree_SiPM_final);

    while (Reader->Next())
    {
        Tree_Silicon_final = vector<Signal>(Tree_Silicon->begin(), Tree_Silicon->end());
        if (Tree_SiPM->GetSize() != 0)
        {
            Tree_SiPM_final[0] = Signal();
            //////MERGING SiPMs
            for (auto &signal : *Tree_SiPM)
            {
                Tree_SiPM_final[0] += SiPMsMatching[GetDetectorChannel(signal.Label)]*signal;
            }
            HSiPM->Fill(Tree_SiPM_final[0].Channel);
        }  
        Tree_Final->Fill();      
    }

    WriteHistograms();

    OutputFile->cd();
    Tree_Final->Write();

    WriteTime(Cal_File, OutputFile);
    OutputFile->Close();
    Cal_File->Close();

    return 0;
}
