#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>

#include <iostream>

using namespace std;

void merge2_convert_sweights(int year, uint nsubs)
{

  auto fin = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/%i_data_beforsel.root",year),"READ");
  auto ntuple = (TTree*)fin->Get("ntuple");

  double val_sw = 0.;
  double runN = 0.;
  Long64_t eventN = 0;
  ntuple->SetBranchAddress("nsig_sw",&val_sw);
  ntuple->SetBranchAddress("runN",&runN);
  ntuple->SetBranchAddress("eventN",&eventN);
  int nev = ntuple->GetEntries();

  auto fout = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/sWeightsConversion/%i_data_beforsel_posWei_mod%i.root",year,nsubs),"RECREATE");
  fout->cd();
  auto tout = new TTree("ntuple_posWei","ntuple");
  double nsig_sw_pos = 0.;
  tout->Branch("runN", &runN, "runN/D");
  tout->Branch("eventN", &eventN, "eventN/L");
  tout->Branch("nsig_sw_pos", &nsig_sw_pos, "nsig_sw_pos/D");

  vector<TFile*> fin_sub (nsubs);
  vector<TTree*> ntuple_sub (nsubs);
  for (uint subs=0; subs<nsubs; ++subs) {
    fin_sub[subs] = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/sWeightsConversion/%i_data_beforsel_%imod%i_v3.root",year,subs,nsubs),"READ");
    ntuple_sub[subs] = (TTree*)fin_sub[subs]->Get(Form("ntuple_%imod%i",subs,nsubs));
    ntuple_sub[subs]->SetBranchAddress("nsig_sw_pos",&nsig_sw_pos);

    if (nev != ntuple_sub[subs]->GetEntries()) {
      cout<<"Error: wrong ntuple dimension in sub: "<<subs<<endl;
      return;
    }

  }

  for (int i=0; i<nev; ++i) {
    ntuple->GetEntry(i);
    ntuple_sub[i%nsubs]->GetEntry(i);
    if (nsig_sw_pos<0)
      cout<<"Error: negative weight "<<nsig_sw_pos<<" for event "<<i<<endl;
    tout->Fill();
  }

  for (uint subs=0; subs<nsubs; ++subs) {
    fin_sub[subs]->Close();
  }

  fout->cd();
  tout->Write();
  fout->Close();
  
  fin->Close();

}

int main(int argc, char** argv)
{

  int year = 2016;
  uint nsubs = 1;

  if ( argc > 1 ) nsubs = atoi(argv[1]);
  if ( argc > 2 ) year = atoi(argv[2]);

  merge2_convert_sweights(year,nsubs);

  return 0;

}
