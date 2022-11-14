#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>

#include <iostream>

using namespace std;

void plot_convert_sweights(int year, int divi)
{

  auto fin = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/%i_data_beforsel.root",year),"READ");
  // auto fin2 = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/sWeightsConversion/%i_data_beforsel_posWei_mod%i.root",year,nsubs),"READ");
  auto ntuple = (TTree*)fin->Get("ntuple");
  ntuple->SetName("ntuData");
  ntuple->AddFriend("ntuple_posWei",Form("/eos/user/a/aboletti/BdToKstarMuMu/sWeightsConversion/%i_data_beforsel_posWei_div%i.root",year,divi));

  auto finMC = new TFile(Form("/eos/cms/store/group/phys_bphys/fiorendi/p5prime/ntuples/preselection/scale_median_and_add_vars/MC_JPSI_%i_preBDT_scale_add_vars_median_no_mass_window.root",year),"READ");
  auto ntupleMC = (TTree*)finMC->Get("ntuple");

  // vector<string> vars = {"bVtxCL", "bLBS", "bLBSE" ,"bCosAlphaBS", "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt","kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk1DCABSE","kstTrk2DCABS","kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta","sum_isopt_04"};
  vector<string> vars = {"bVtxCL", "bLBS", "bLBSE" ,"bCosAlphaBS", "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt","kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk2DCABS","sum_isopt_04",
    "kstTrk1DCABSE","kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta"};
  vector<double> vals(vars.size(),0.);
  double val_sw = 0.;
  double val_sw_pos = 0.;
  Long64_t eventN = 0;

  vector<double> var_min(vars.size(),999.);
  vector<double> var_max(vars.size(),-999.);
  vector<double> x2sum (vars.size(),0.);
  vector<double> xsum (vars.size(),0.);
  vector<uint> xcnt (vars.size(),0);

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  ntuple->SetBranchStatus("*",0);
  ntupleMC->SetBranchStatus("*",0);
  for (uint i=0; i<vars.size(); ++i) {
    ntuple->SetBranchStatus(vars[i].c_str(),1);
    ntuple->SetBranchAddress(vars[i].c_str(),&vals[i]);
    ntupleMC->SetBranchStatus(vars[i].c_str(),1);
    ntupleMC->SetBranchAddress(vars[i].c_str(),&vals[i]);
  }
  ntuple->SetBranchStatus("nsig_sw",1);
  ntuple->SetBranchStatus("nsig_sw_pos",1);
  ntuple->SetBranchStatus("eventN",1);
  ntuple->SetBranchAddress("nsig_sw",&val_sw);
  ntuple->SetBranchAddress("nsig_sw_pos",&val_sw_pos);
  ntuple->SetBranchAddress("eventN",&eventN);

  Long64_t nev = ntuple->GetEntries();
  for (int i=0; i<nev; i++) {
    if (i%divi != 0) continue;
    ntuple->GetEntry(i);
    if (eventN%2 != 0) continue;
    for (uint iVar=0; iVar<vars.size(); ++iVar)
      if (vals[iVar]==vals[iVar]) {
	if (var_min[iVar]>vals[iVar]) var_min[iVar] = vals[iVar];
	if (var_max[iVar]<vals[iVar]) var_max[iVar] = vals[iVar];
	xcnt[iVar]++;
	xsum[iVar] += vals[iVar];
	x2sum[iVar] += vals[iVar]*vals[iVar];
      }
  }

  vector<TH1D*> hori(0);
  vector<TH1D*> hunw(0);
  vector<TH1D*> hxgb(0);
  vector<TH1D*> hpos(0);
  vector<TH1D*> hMCu(0);

  vector<TH1D*> runw(0);
  vector<TH1D*> rxgb(0);
  vector<TH1D*> rpos(0);
  vector<TH1D*> rMCu(0);

  for (uint iVar=0; iVar<vars.size(); ++iVar) {
    // cout<<xsum[iVar]<<"\t"<<xcnt[iVar]<<"\t"<<x2sum[iVar]<<"\t"<<xsum[iVar]/xcnt[iVar]<<"\t"<<sqrt(x2sum[iVar]/xcnt[iVar] - xsum[iVar]*xsum[iVar]/xcnt[iVar]/xcnt[iVar])<<endl;
    double varlow  = max( var_min[iVar], xsum[iVar]/xcnt[iVar] - 4 * sqrt(x2sum[iVar]/xcnt[iVar] - xsum[iVar]*xsum[iVar]/xcnt[iVar]/xcnt[iVar]) );
    double varhigh = min( var_max[iVar], xsum[iVar]/xcnt[iVar] + 4 * sqrt(x2sum[iVar]/xcnt[iVar] - xsum[iVar]*xsum[iVar]/xcnt[iVar]/xcnt[iVar]) );
    hori.push_back( new TH1D( Form("hori%i",iVar), Form("%s distribution comparison;%s;Events",vars[iVar].c_str(),vars[iVar].c_str()), 100, varlow, varhigh ) );
    hunw.push_back( new TH1D( Form("hunw%i",iVar), Form("%s distribution comparison;%s;Events",vars[iVar].c_str(),vars[iVar].c_str()), 100, varlow, varhigh ) );
    hxgb.push_back( new TH1D( Form("hxgb%i",iVar), Form("%s distribution comparison;%s;Events",vars[iVar].c_str(),vars[iVar].c_str()), 100, varlow, varhigh ) );
    hpos.push_back( new TH1D( Form("hpos%i",iVar), Form("%s distribution comparison;%s;Events",vars[iVar].c_str(),vars[iVar].c_str()), 100, varlow, varhigh ) );
    hMCu.push_back( new TH1D( Form("hMCu%i",iVar), Form("%s distribution comparison;%s;Events",vars[iVar].c_str(),vars[iVar].c_str()), 100, varlow, varhigh ) );
    runw.push_back( new TH1D( Form("runw%i",iVar), Form(";%s;Ratio",vars[iVar].c_str()), 100, varlow, varhigh ) );
    rxgb.push_back( new TH1D( Form("rxgb%i",iVar), Form(";%s;Ratio",vars[iVar].c_str()), 100, varlow, varhigh ) );
    rpos.push_back( new TH1D( Form("rpos%i",iVar), Form(";%s;Ratio",vars[iVar].c_str()), 100, varlow, varhigh ) );
    rMCu.push_back( new TH1D( Form("rMCu%i",iVar), Form(";%s;Ratio",vars[iVar].c_str()), 100, varlow, varhigh ) );
  }
  
  for (int i=0; i<nev; i++) {
    ntuple->GetEntry(i);
    for (uint iVar=0; iVar<vars.size(); ++iVar)
      if (vals[iVar]==vals[iVar]) {
	hori[iVar]->Fill(vals[iVar],val_sw);
	hxgb[iVar]->Fill(vals[iVar],val_sw>0?val_sw:0);
	hunw[iVar]->Fill(vals[iVar]);
      }
    if (i%divi != 0) continue;
    if (eventN%2 != 0) continue;
    for (uint iVar=0; iVar<vars.size(); ++iVar)
      if (vals[iVar]==vals[iVar])
	hpos[iVar]->Fill(vals[iVar],val_sw_pos);
  }

  for (int i=0; i<ntupleMC->GetEntries(); i++) {
    ntupleMC->GetEntry(i);
    for (uint iVar=0; iVar<vars.size(); ++iVar)
      if (vals[iVar]==vals[iVar])
	hMCu[iVar]->Fill(vals[iVar]);
  }

  TLine line;
  line.SetLineColor(13);
  line.SetLineStyle(2);

  for (uint iVar=0; iVar<vars.size(); ++iVar) {

    auto can = new TCanvas(Form("can%i",iVar),Form("can%i",iVar),1000,1000);
    can->cd();
    hunw[iVar]->Scale(hori[iVar]->Integral()/hunw[iVar]->Integral());
    hxgb[iVar]->Scale(hori[iVar]->Integral()/hxgb[iVar]->Integral());
    hpos[iVar]->Scale(hori[iVar]->Integral()/hpos[iVar]->Integral());
    hMCu[iVar]->Scale(hori[iVar]->Integral()/hMCu[iVar]->Integral());
    hpos[iVar]->SetLineColor(2);
    hori[iVar]->SetLineColor(1);
    hunw[iVar]->SetLineColor(4);
    hxgb[iVar]->SetLineColor(8);
    hMCu[iVar]->SetLineColor(4);
    hori[iVar]->SetMinimum(0);
    hori[iVar]->Draw();
    // hunw[iVar]->Draw("same");
    hxgb[iVar]->Draw("same");
    hpos[iVar]->Draw("same");
    hMCu[iVar]->Draw("same");

    TLegend leg(0.3,0.2);
    leg.AddEntry(hori[iVar],"data sPlot","lep");
    leg.AddEntry(hxgb[iVar],"data sPlot (pos w only)","lep");
    leg.AddEntry(hpos[iVar],"data converted sWeights","lep");
    leg.AddEntry(hMCu[iVar],"MC unweighted","lep");
    leg.Draw();

    can->SaveAs(Form("plot_d/%s_dist_div%i_v2.pdf",vars[iVar].c_str(),divi));

    for (int iBin=1; iBin<=hori[iVar]->GetNbinsX(); ++iBin) {
      double ref = hori[iVar]->GetBinContent(iBin);
      double refe2 = hori[iVar]->GetBinError(iBin)*hori[iVar]->GetBinError(iBin)/ref/ref;
      rpos[iVar]->SetBinContent(iBin,hpos[iVar]->GetBinContent(iBin)/ref);
      rxgb[iVar]->SetBinContent(iBin,hxgb[iVar]->GetBinContent(iBin)/ref);
      runw[iVar]->SetBinContent(iBin,hunw[iVar]->GetBinContent(iBin)/ref);
      rMCu[iVar]->SetBinContent(iBin,hMCu[iVar]->GetBinContent(iBin)/ref);
      rpos[iVar]->SetBinError(iBin,rpos[iVar]->GetBinContent(iBin) * sqrt(hpos[iVar]->GetBinError(iBin)*hpos[iVar]->GetBinError(iBin)/hpos[iVar]->GetBinContent(iBin)/hpos[iVar]->GetBinContent(iBin) + refe2));
      rxgb[iVar]->SetBinError(iBin,rxgb[iVar]->GetBinContent(iBin) * sqrt(hxgb[iVar]->GetBinError(iBin)*hxgb[iVar]->GetBinError(iBin)/hxgb[iVar]->GetBinContent(iBin)/hxgb[iVar]->GetBinContent(iBin) + refe2));
      runw[iVar]->SetBinError(iBin,runw[iVar]->GetBinContent(iBin) * sqrt(hunw[iVar]->GetBinError(iBin)*hunw[iVar]->GetBinError(iBin)/hunw[iVar]->GetBinContent(iBin)/hunw[iVar]->GetBinContent(iBin) + refe2));
      rMCu[iVar]->SetBinError(iBin,rMCu[iVar]->GetBinContent(iBin) * sqrt(hMCu[iVar]->GetBinError(iBin)*hMCu[iVar]->GetBinError(iBin)/hMCu[iVar]->GetBinContent(iBin)/hMCu[iVar]->GetBinContent(iBin) + refe2));
    }

    auto can2 = new TCanvas(Form("can2%i",iVar),Form("can2%i",iVar),1000,1000);
    can2->cd();
    rpos[iVar]->SetLineColor(2);
    runw[iVar]->SetLineColor(4);
    rxgb[iVar]->SetLineColor(8);
    rMCu[iVar]->SetLineColor(4);
    rpos[iVar]->Draw();
    line.DrawLine(rpos[iVar]->GetXaxis()->GetXmin(),1,rpos[iVar]->GetXaxis()->GetXmax(),1);
    rxgb[iVar]->Draw("same");
    // runw[iVar]->Draw("same");
    rMCu[iVar]->Draw("same");
    rpos[iVar]->GetYaxis()->SetRangeUser(0.01,1.99);
    can2->SaveAs(Form("plot_d/%s_distratio_div%i_v2.pdf",vars[iVar].c_str(),divi));

  }

  // for (uint i=0; i<vars.size(); ++i) {
  //   // auto hori = new TH1D(Form("hori%i",i),Form("%s distribution;%s;Events",vars[i].c_str(),vars[i].c_str()));
  //   auto can = new TCanvas(Form("can%i",i),Form("can%i",i));
  //   can->cd();
  //   ntuple->Draw(Form("%s>>hori%i",vars[i].c_str(),i),"nsig_sw");
  //   ntuple->Draw(Form("%s>>hpos%i",vars[i].c_str(),i),"nsig_sw_pos","same");
  //   ((TH1F*)gDirectory->Get(Form("hpos%i",i)))->SetLineColor(2);
  //   can->SaveAs(Form("plot_d/%s_dist_mod%i.pdf",vars[i].c_str(),nsubs));
  // }
  
  fin->Close();

}

int main(int argc, char** argv)
{

  int nsubs = 1;
  int year = 2016;

  if ( argc > 1 ) nsubs = atoi(argv[1]);
  if ( argc > 2 ) year = atoi(argv[2]);

  plot_convert_sweights(year,nsubs);

  return 0;

}
