#define analyzer_cxx
#include "analyzer.h"
#include <TH2.h>
#include <TStyle.h>


void analyzer::Begin(TTree * /*tree*/)
{
    // Job config
    TString option  = GetOption();
    TObjArray *args = (TObjArray*)option.Tokenize(" ");
    string suffix   = (string)((TObjString*)args->At(0))->GetString();

    histoFile = new TFile(("histos/test_" + suffix + ".root").c_str(), "RECREATE");
    histoFile->mkdir(suffix.c_str(), suffix.c_str());
    //histoFile->cd(suffix.c_str);

    h1_Met  = new TH1D("h1_Met", "MET", 20, 0., 100.);
}

Bool_t analyzer::Process(Long64_t entry)
{
    GetEntry(entry);

    h1_Met->Fill(met);

    if (jetMult > 2 && HT < 60)
        cout << HT << endl;

    //!! Z-veto !!//
    /*
    if (
            leptons.size() == 2 
            && leptons[0].Charge() == leptons[1].Charge()
            && leptons[0].Type() == leptons[1].Type()
            && fabs((leptons[0] + leptons[1]).M() - 91.2) < 15. 
       ) 
        return true;
    else if (leptons.size() == 3 && (zTagged || (dileptonMassOS > 40 && fabs(trileptonMass - 91.2) < 7.5))) 
        return true;


    //!! Require at least one b-jet !!//
    if (bJetsM.size() + jets.size() <= 1) return true;

    if (jetMult > 1 && HT < 60)
        cout << HT << endl;


    //!! MET cut !!//
    if (leptons.size() == 2){
        if (leptons[0].Charge() == leptons[1].Charge()) 
            if (recoMET->Mod() < metCut[0])
                return true;
    } else if (leptons.size() == 3) {
        if (recoMET->Mod() < metCut[1]) 
            return true;
    }

    //!! HT cut !!//
    if (leptons.size() == 2){
        if (leptons[0].Charge() == leptons[1].Charge()) 
            if (sqrt(HT) < htCut[0])
                return true;
    } else if (leptons.size() == 3) {
        if (sqrt(HT) < htCut[1]) 
            return true;
    }
    */


    return kTRUE;
}

void analyzer::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

    histoFile->Write();
    histoFile->Close();


}
