#include "Matcher.hh"
#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include <ctime>

using namespace std;
using namespace ROOT::Math;

int main(int argc, char *argv[])
{

    clock_t start = clock(), Current;
    string filename = "run_069_90Sr";
    TFile *OutputFile = new TFile(("./../../../../../../../mnt/hgfs/shared-2/"+filename+"_matched.root").c_str(), "RECREATE");
    TFile *Cal_File = new TFile(("./../../../../../../../mnt/hgfs/shared-2/"+filename+"_grouped.root").c_str(), "READ");
    InitDetectors("../../Analysis/ReadFaster/Data/detectors.dat");
    /////// Coef between Low and High /////
    cout<<"<SAM> Fitting gain coefficients Low/High"<<endl;
    dir_Chi2SiPMLowHigh = OutputFile->mkdir("SiPM_LowHigh_Matching");
    dir_Chi2SiPMLowHigh->cd();
    TGraphErrors *Result_SiPM_LowHigh = new TGraphErrors();
    Result_SiPM_LowHigh->SetDrawOption("*AP");
    Result_SiPM_LowHigh->SetName("SiPM_LowHigh_Matching");
    Result_SiPM_LowHigh->SetTitle("SiPM_LowHigh_Matching;SiPM Channel;SiPM_LowHigh_Matching");
    int result_counter=0;
    for (int i = 1; i <= BETA_SIZE; i++)
    {
        TCanvas *canvas = (TCanvas *)Cal_File->Get(("SiPM_Multiplicities/SiPM_" + to_string(i)).c_str());
        TGraph *graph = (TGraph *)canvas->GetPrimitive(("SiPM_" + to_string(i)).c_str());

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
            TCanvas *c = new TCanvas(("SiPM_" + to_string(i)).c_str(), ("SiPM_" + to_string(i)).c_str(), 800, 600);
            c->cd();
            graph->Draw("*AP");
            f1->Draw("SAME");
            dir_Chi2SiPMLowHigh->cd();
            c->Write();
            delete c;
        }
        else
        {
            par = 0;
            err = 0;

            OutputFile->cd();
            TCanvas *c = new TCanvas(("SiPM_" + to_string(i)).c_str(), ("SiPM_" + to_string(i)).c_str(), 800, 600);
            c->cd();
            dir_Chi2SiPMLowHigh->cd();
            c->Write();
            delete c;
        }

        
        Result_SiPM_LowHigh->SetPoint(result_counter, i, par);
        SiPMLowHighMatching[i] = par;
        Result_SiPM_LowHigh->SetPointError(result_counter, 0, err);
        result_counter++;
    }

    dir_Chi2SiPMLowHigh->cd();
    Result_SiPM_LowHigh->Write();

    OutputFile->cd();

    /////// Merging Low and High, Reconstruct SIPM losses, apply matching on Silicon ///////
    cout<<"<SAM> Merging and Reconstructing Low/High"<<endl;
    InitHistograms();
    dir_SiPMsReconstruction = OutputFile->mkdir("SiPMs_Reconstruction");
    Tree_Read = (TTree *)Cal_File->Get("Tree");
    Reader = new TTreeReader(Tree_Read);
    TTreeReaderArray<Signal> *Tree_SiPMHigh = new TTreeReaderArray<Signal>(*Reader, "Tree_SiPMHigh");
    TTreeReaderArray<Signal> *Tree_SiPMLow = new TTreeReaderArray<Signal>(*Reader, "Tree_SiPMLow");

    Tree_Write = new TTree("Tree_LowHighMatched", "Tree_LowHighMatched");
    vector<Signal> Tree_SiPM_W;
    Tree_Write->Branch("Tree_SiPM", &Tree_SiPM_W);

    TProfile *profile[BETA_SIZE+1];
    
    for (int i=1; i <= BETA_SIZE; i++)
    {
        profile[i] = (TProfile *)Cal_File->Get(("SiPM_Multiplicities/PSiPM_" + to_string(i)).c_str());
    }

    vector<int> vec_det;
        vector<Signal> new_SiPMHigh;
        vector<Signal> new_SiPMLow;
    
    while (Reader->Next() && Reader->GetCurrentEntry() < 3000000)
    {

        ProgressBar(Reader->GetCurrentEntry(), Reader->GetEntries(), start, Current, "");

        Tree_SiPM_W.clear();
        vec_det.clear();
        new_SiPMHigh.clear();
        new_SiPMLow.clear();


        /////FOR SAVING BEFORE MERGING//////
        for (auto &signal : *Tree_SiPMHigh)
        {
            HSiPMHigh_Channel_all[GetDetectorChannel(signal.Label)]->Fill( 1 / SiPMLowHighMatching[GetDetectorChannel(signal.Label)] * signal.Channel);
        }

        for (auto &signal : *Tree_SiPMLow)
        {
            HSiPMLow_Channel_all[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
        }
        ///////////////////////////////////

        // RECONSTRUCTING SiPM LOSSES ////////////////////////////////////////////////////////////////////
        
        // cout << "Event : " << Reader->GetCurrentEntry() << " (Si : "<< (*Tree_Silicon)[0].Channel << " )"<<endl;
        // for (auto &ith : *Tree_SiPMHigh)
        // {
        //     cout<<ith.Label<< " : " <<ith.Channel<<endl;
        // }
        // for (auto &itl : *Tree_SiPMLow)
        // {
        //     cout<<itl.Label<< " : " <<itl.Channel<<endl;
        // }

        for (auto &ith : *Tree_SiPMHigh)
        {
            for (auto &itl : *Tree_SiPMLow)
            {
                if (GetDetectorChannel(ith.Label) == GetDetectorChannel(itl.Label))
                {
                    // cout<<GetDetectorChannel(ith.Label)<< " : " <<ith.Channel<< "  " <<itl.Channel<< setprecision(15)<< "    TIME_DIFF : "<< ith.Time-itl.Time<<endl;
                    vec_det.push_back(GetDetectorChannel(ith.Label));
                }
            }
        }


        for (auto &ith : *Tree_SiPMHigh)
        {
            if (std::find(vec_det.begin(), vec_det.end(), GetDetectorChannel(ith.Label)) == vec_det.end())
            {
                HSiPMHigh_False->Fill(ith.Channel);
                // HSiPMHigh_Time_False->Fill(ith.Time - (*Tree_Silicon)[0].Time);
                HSiPMHigh_SiPM_False->Fill(GetDetectorChannel(ith.Label));


                //////////RECONSTRUCTING LOW SiPM//////////
                Signal new_low;
                if (ith.Channel < 2000000)
                {
                    for (int i = 1; i <= profile[GetDetectorChannel(ith.Label)]->GetNbinsX(); i++)
                    {
                        if (profile[GetDetectorChannel(ith.Label)]->GetBinContent(i) > ith.Channel)
                        {
                            new_low = Signal(ith.Label, ith.Time, profile[GetDetectorChannel(ith.Label)]->GetBinCenter(i));
                            break;
                        }
                    }
                }
                else
                {
                    double sum = 0.;
                    if (Tree_SiPMLow->GetSize() != 0)
                    {
                        for (auto &itl : *Tree_SiPMLow)
                        {
                            sum += itl.Channel;
                        }
                        new_low = Signal(ith.Label, ith.Time, sum/Tree_SiPMLow->GetSize());  
                    }
                    else
                    {
                        new_low = Signal(ith.Label, ith.Time, 0);
                    }
                }

                new_SiPMLow.push_back(new_low);

                HSiPMLow_New->Fill(new_low.Channel);
                // HSiPMLow_Time_New->Fill(new_low.Time - (*Tree_Silicon)[0].Time);
                HSiPMLow_SiPM_New->Fill(GetDetectorChannel(new_low.Label));
            }
        }


        for (auto &itl : *Tree_SiPMLow)
        {
            if (std::find(vec_det.begin(), vec_det.end(), GetDetectorChannel(itl.Label)) == vec_det.end())
            {
                HSiPMLow_False->Fill(itl.Channel);
                // HSiPMLow_Time_False->Fill(itl.Time - (*Tree_Silicon)[0].Time);
                HSiPMLow_SiPM_False->Fill(GetDetectorChannel(itl.Label));

                //////////RECONSTRUCTING HIGH SiPM//////////
                Signal new_high = Signal(HighLowLabelConversion(itl.Label), itl.Time, profile[GetDetectorChannel(itl.Label)]->GetBinContent(profile[GetDetectorChannel(itl.Label)]->FindBin(itl.Channel)));
                if (new_high.Channel == 0 && itl.Channel > 400000)
                {
                    new_high.Channel = 400000;
                }
                new_SiPMHigh.push_back(new_high);
                
                HSiPMHigh_New->Fill(new_high.Channel);
                // HSiPMHigh_Time_New->Fill(new_high.Time - (*Tree_Silicon)[0].Time);
                HSiPMHigh_SiPM_New->Fill(GetDetectorChannel(new_high.Label));
            }
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        //////MATCHING AND MERGING SiPM Low and High
        bool doHigh = true;

        for (auto &signal : *Tree_SiPMHigh)
        {
            if (1 / SiPMLowHighMatching[GetDetectorChannel(signal.Label)] * signal.Channel > 180000)
            {
                doHigh = false;
            }
        }

        for (auto signal : new_SiPMHigh)
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

            // for (auto &signal : new_SiPMHigh)
            // {
            //     signal.Channel = 1 / SiPMLowHighMatching[GetDetectorChannel(signal.Label)] * signal.Channel;
            //     Tree_SiPM_W.push_back(signal);

            //     HSiPMHigh_Channel[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
            //     HSiPM_Channel[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
            // }

        }
        else
        {
            for (auto &signal : *Tree_SiPMLow)
            {
                Tree_SiPM_W.push_back(signal);
                HSiPMLow_Channel[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
                HSiPM_Channel[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
            }

            // for (auto &signal : new_SiPMLow)
            // {
            //     Tree_SiPM_W.push_back(signal);
            //     HSiPMLow_Channel[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
            //     HSiPM_Channel[GetDetectorChannel(signal.Label)]->Fill(signal.Channel);
            // }
        }
        Tree_Write->Fill();

        
    }

    Tree_Write->Write();

    // ////// Matching all SiPMs/////
    cout << "<SAM> Matching SiPMs" << endl;
    dir_Chi2SiPMs = OutputFile->mkdir("SiPMs_Matching");
    dir_Chi2SiPMs->cd();
    TGraphErrors *Result_SiPMs = new TGraphErrors();
    Result_SiPMs->SetDrawOption("*AP");
    Result_SiPMs->SetName("SiPMs_Matching");
    Result_SiPMs->SetTitle("SiPMs_Matching;Detector Label;SiPMs_Matching");
    result_counter = 0;
    Ref_Hist = HSiPM_Channel[1]; //////////////////ref is 1
    range_low = 400000;
    range_high = 1000000;
    Tree = (TTree *)OutputFile->Get("Tree_LowHighMatched");
    Reader = new TTreeReader(Tree_Write);

    for (size_t i = 1; i <= BETA_SIZE; ++i)
    {

    // graph = new TGraph();
    // graph->SetTitle("Gain_Matching;Proportionality; #chi^{2}");
    // counter = 0;

    // current_SiPM = i;

    // //////////////////// MINIMIZER ////////////////////
    // ////// First
    // Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    // ROOT::Math::Functor functor(&Chi2SiPM, 1);
    // minimizer->SetFunction(functor);
    // minimizer->SetLimitedVariable(0, "Proportionality", 1., .05, 0.6, 1.2);
    // minimizer->SetPrecision(0.05);
    // minimizer->Minimize();
    // const double *bestPar = minimizer->X();

    // ////// Second
    // minimizer->SetLimitedVariable(0, "Proportionality", bestPar[0], .05, 0.6, 1.2);
    // minimizer->SetPrecision(1e-5);
    // minimizer->Minimize();
    // bestPar = minimizer->X();

    // ////// Third
    // minimizer->SetLimitedVariable(0, "Proportionality", bestPar[0], .005, 0.6, 1.2);
    // minimizer->SetPrecision(1e-10);
    // minimizer->Minimize();
    // bestPar = minimizer->X();
    // double error = minimizer->Errors()[0];
    // std::cout << "SiPM " + to_string(i) << "\t CHI2 : " << chi2 << "\t\t Best Proportionality : " << bestPar[0] << " +/- " << error << std::endl;

    // //////////////////// SAVE ////////////////////
    // TCanvas *c = new TCanvas(("HSiPM_Matched" + to_string(i) + "_Channel").c_str(), ("HSiPM_Matched" + to_string(i) + "_Channel").c_str(), 800, 600);
    // c->Divide(1, 3);
    // c->cd(1);
    // Ref_Hist->GetXaxis()->UnZoom();
    // Ref_Hist->Draw("HIST");

    // Reader->Restart();
    // TTreeReaderArray<Signal> *SiPM = new TTreeReaderArray<Signal>(*Reader, "Tree_SiPM");
    // TH1D *TreeHist = (TH1D *)Ref_Hist->Clone((detectorName[i] + "_" + to_string(bestPar[0])).c_str());
    // TreeHist->Reset();
    // while (Reader->Next())
    // {
    //     for (size_t j = 0; j < SiPM->GetSize(); j++)
    //     {
    //         if (current_SiPM == GetDetectorChannel((*SiPM)[j].Label))
    //             TreeHist->Fill((*SiPM)[j].Channel * bestPar[0]);
    //     }
    // }

    // TreeHist->SetLineColor(kRed);
    // TreeHist->GetXaxis()->UnZoom();
    // TreeHist->Draw("SAME");

    // c->cd(2);
    // TH1D *h = (TH1D *)Ref_Hist->Clone("Diff");
    // h->Add(TreeHist, -1.0);
    // h->Draw("EP");

    // c->cd(3);
    // graph->Draw("*AP");
    // OutputFile->cd();
    // dir_Chi2SiPMs->cd();
    // c->Write();

    // result_counter++;
    // Result_SiPMs->SetPoint(result_counter, i, bestPar[0]);
    // SiPMsMatching[i] = bestPar[0];
    // Result_SiPMs->SetPointError(result_counter, 0, error);
    // delete c;
    }
    //     dir_Chi2SiPMs->cd();
    // Result_SiPMs->Write();

    SiPMsMatching[1] = 1;
    SiPMsMatching[2] = 0.94;
    SiPMsMatching[3] = 1.01;
    SiPMsMatching[4] = 1.;
    SiPMsMatching[5] = 1.03;
    SiPMsMatching[6] = 0.81;
    SiPMsMatching[7] = 1;
    SiPMsMatching[8] = 0.85;
    SiPMsMatching[9] = 1;


    ////other technique to maching sipms///////////////////////////////////////////////////////////
    // TProfile *graph1[BETA_SIZE+1];
    // int counters[BETA_SIZE+1];
    // for (int i = 0; i <= BETA_SIZE; i++)
    // {
    //     graph1[i] = new TProfile(("PSiPM_" + to_string(i)).c_str(), ("PSiPM_" + to_string(i)).c_str(), 500, 0, eLowMax, 0, eLowMax);
    // }

    // TTreeReaderArray<Signal> *SiPM = new TTreeReaderArray<Signal>(*Reader, "Tree_SiPM");
    // // int TotalEntries = Reader->GetEntries();
    // // string Prefix = "SiPMs Matching";
    // // clock_t start = clock(), Current;
    // while (Reader->Next())
    // {
    //     Signal Ref;
    //     // ULong64_t cEntry = Reader->GetCurrentEntry();
    //         // if ((cEntry % 100 == 0 && cEntry > 0) || (cEntry == TotalEntries-1))
    //         // {
    //         //     ProgressBar(cEntry, TotalEntries, start, Current, Prefix);
    //         // }
    //     if (SiPM->GetSize() == 7)
    //     {
    //         for (int i = 0; i < 7; i++)
    //         {
    //             if (GetDetectorChannel((*SiPM)[i].Label) == 1)
    //             {
    //                 Ref = (*SiPM)[i];
    //                 break;
    //             }
    //         }
    //         for (int i = 0; i < 7; i++)
    //         {
    //             if (GetDetectorChannel((*SiPM)[i].Label) != 1 && Ref.isValid)
    //             {
    //                 graph1[GetDetectorChannel((*SiPM)[i].Label)]->Fill((*SiPM)[i].Channel, Ref.Channel);
    //                 counters[GetDetectorChannel((*SiPM)[i].Label)]++;
    //             }
    //         }
    //     }
    // }

    // double new_matching[BETA_SIZE+1];
    // TCanvas *c = new TCanvas("SiPMs_Matching1", "SiPMs_Matching1", 800, 600);
    // c->Divide(3, 3);
    
    // for (size_t i = 1; i <= BETA_SIZE; ++i)
    // {
    //     if (graph1[i]->GetEntries() == 0)
    //     {
    //         new_matching[i] = 0;//
    //         continue;
    //     }
    //     c->cd(i);
    //     graph1[i]->Draw();
    //     TF1 *fit = new TF1("fit", "[0]*x");
    //     graph1[i]->Fit(fit);
    //     new_matching[i] = fit->GetParameter(0);
        
    // }
    // OutputFile->cd();
    // c->Write();
    // delete c;

    //////////////////// SAVE ////////////////////
    // for (size_t i = 1; i <= BETA_SIZE; ++i)
    // {
    //     TCanvas *c8 = new TCanvas(("HSiP_Matched" + to_string(i) + "_Channel").c_str(), ("HSiP_Matched" + to_string(i) + "_Channel").c_str(), 800, 600);
    //     c8->cd();
    //     current_SiPM = i;
    //     Ref_Hist->GetXaxis()->UnZoom();
    //     Ref_Hist->Draw("HIST");

    //     Reader->Restart();
    //     TTreeReaderArray<Signal> *SiPM = new TTreeReaderArray<Signal>(*Reader, "Tree_SiPM");
    //     TH1D *TreeHist = (TH1D *)Ref_Hist->Clone((detectorName[i] + "_" + to_string(new_matching[i])).c_str());
    //     TreeHist->Reset();
    //     while (Reader->Next())
    //     {
    //         for (size_t j = 0; j < SiPM->GetSize(); j++)
    //         {
    //             if (current_SiPM == GetDetectorChannel((*SiPM)[j].Label))
    //                 TreeHist->Fill((*SiPM)[j].Channel * new_matching[current_SiPM]);
    //         }
    //     }

    //     cout<< Ref_Hist->Chi2Test(TreeHist, "CHI2/NDF")<<endl;
    //     TreeHist->SetLineColor(kRed);
    //     TreeHist->GetXaxis()->UnZoom();
    //     TreeHist->Draw("SAME");
    //     OutputFile->cd();
    //     c8->Write();
    //     delete c8;
    // }
    ////////////////////////////////////////////////////////////////////////////////////////

    ///// Apply Matching between SiPMs and Merge them ///////
    cout << "<SAM> Merging SiPMs" << endl;
    Reader->Restart();
    TTreeReaderArray<Signal> *Tree_SiPM = new TTreeReaderArray<Signal>(*Reader, "Tree_SiPM");

    Tree_Final = new TTree("Tree", "Tree");
    vector<Signal> Tree_SiPM_final;
    Tree_Final->Branch("Tree_SiPM", &Tree_SiPM_final);

    while (Reader->Next())
    {
        Tree_SiPM_final.clear();

        ////////////////CUMULATIVE
        // Tree_Silicon_final = vector<Signal>(Tree_Silicon->begin(), Tree_Silicon->end());

        if (Tree_SiPM->GetSize() != 0)
        {
            for (auto &signal : *Tree_SiPM)
            {
                signal.Channel = SiPMsMatching[GetDetectorChannel(signal.Label)] * signal.Channel;
                Tree_SiPM_final.push_back(signal);
                HSiPM->Fill(signal.Channel);
            }
        }
        else
        {
            Tree_SiPM_final = vector<Signal>();
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
