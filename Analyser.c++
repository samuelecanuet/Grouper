#include "Analyser.hh"


int main(int argc, char *argv[])
{
    if (argc == 2)
    {
        if (argv[1] == "Al")
        {
            Catcher = "Al";
        }
        else
        {
            Catcher = "Mylar";
        }
    }
    else
    {
        Catcher = "Al";
    }



    ////// TAKING PARAMETERS FROM CONVOLUTION OF SIPMs //////
    TFile* Merged_File = new TFile(("Merged/runs_32Ar_merged.root").c_str(), "READ");
    Convolution_Parameters = Extract_Parameters(Merged_File);


    ////// TAKING DATA FROM SIMUMLATION //////
    a = {-1, 0, 1};
    Simulation_File[0] = new TFile(("../../../../../../mnt/hgfs/shared-2/32Ar_a1_b0.root").c_str(), "READ");
    Simulation_File[1] = new TFile(("../../../../../../mnt/hgfs/shared-2/32Ar_a1_b0.root").c_str(), "READ");
    Simulation_File[2] = new TFile(("../../../../../../mnt/hgfs/shared-2/32Ar_a1_b0.root").c_str(), "READ");
    FillHistograms_SIM();
    
    ////// TAKING DATA FROM EXPERIMENT //////
    FillHistograms_EXP();


    ////// COMPUTE ESHIFT FOR SIMULATION //////
    Compute_Eshift();




}

