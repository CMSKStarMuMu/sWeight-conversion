#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>

#include <iostream>

using namespace std;

vector<double> rms(0);

vector<double> val_sw_pos(0);
vector<bool> is_sw_mod(0);
double val_sw_ori = 0.;

float getDist(vector<double>* v1, vector<double>* v2)
{
  float ret = 0;
  for (uint i=0; i<v1->size(); ++i)
    ret += (v1->at(i)-v2->at(i))*(v1->at(i)-v2->at(i))/rms[i]/rms[i];
  return sqrt(ret);
}

float processList(vector<uint> locmap, TTree* ntuple, vector<double>* vals, vector<double>* vals_ori)
{
  vector<int> candidates(0);
  vector<float> dists(0);

  // for (int j=0; j<nev; ++j) {
  float sum = 0.;
  for (uint ij=0; ij<locmap.size(); ++ij) {
    int j = locmap[ij];

    if (val_sw_pos[j]<=0) continue;

    ntuple->GetEntry(j);

    float dist = getDist(vals,vals_ori);
    if ( dists.size() > 0 && dist > dists.back() && sum + val_sw_pos[j] > -1*val_sw_ori )
      continue;
    // if (i==155) cout<<j<<" == "<<(dists.size()>0 ? dist-dists.back() : 999)<<" "<<sum + val_sw_pos[j] + val_sw_ori<<endl;

    uint k;
    for (k=0; k<candidates.size(); ++k)
      if (dist < dists[k]) break;
    dists.insert(dists.begin()+k,dist);
    candidates.insert(candidates.begin()+k,j);


    sum = 0.;
    for (k=0; k<candidates.size(); ++k) {
      if (sum>-1*val_sw_ori) {
	candidates.erase(candidates.begin()+k);
	dists.erase(dists.begin()+k);
	--k;
      } else
	sum += val_sw_pos[candidates[k]];
    }

  }

  // cout<<i<<" "<<val_sw_ori<<" - ";
  // for (k=0; k<candidates.size(); ++k)
  //   cout<<candidates[k]<<","<<dists[k]<<","<<val_sw_pos[candidates[k]]<<"\t";
  // cout<<endl;
      
  sum = -1*val_sw_ori;
  for (uint k=0; k<candidates.size(); ++k) {
    if ( sum > val_sw_pos[candidates[k]] ) {
      sum -= val_sw_pos[candidates[k]];
      val_sw_pos[candidates[k]] = 0;
    } else {
      val_sw_pos[candidates[k]] -= sum;
      sum = 0;
      break;
    }
  }

  // cout<<i<<" "<<val_sw_ori<<" - ";
  // for (k=0; k<candidates.size(); ++k)
  //   cout<<candidates[k]<<","<<dists[k]<<","<<val_sw_pos[candidates[k]]<<"\t";
  // cout<<endl;

  return sum;

}

void merge_convert_sweights(int year, uint subs, uint nsubs)
{

  auto fin = new TFile(Form("/eos/user/x/xuqin/workdir/B0KstMuMu/reweight/Tree/final/XGBV5/%i/%i_data_beforsel.root",year,year),"READ");
  // auto fin = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/%i_data_beforsel.root",year),"READ");
  auto ntuple = (TTree*)fin->Get("ntuple");

  // vector<string> vars = {"bVtxCL", "bLBS", "bLBSE" ,"bCosAlphaBS", "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt","kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk1DCABSE","kstTrk2DCABS","kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta","sum_isopt_04"};
  vector<string> vars = {"bVtxCL", "bLBS", "bLBSE" ,"bCosAlphaBS", "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt","kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk2DCABS","sum_isopt_04",
    "kstTrk1DCABSE","kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta"};
  vector<double> vals(vars.size(),0.);
  vector<double> vals_ori(0);
  
  string var_sw = "nsig_sw";
  double val_sw = 0.;

  vector<double> x2sum (vars.size(),0.);
  vector<double> xsum (vars.size(),0.);
  vector<uint> xcnt (vars.size(),0);
  
  for (uint i=0; i<vars.size(); ++i)
    ntuple->SetBranchAddress(vars[i].c_str(),&vals[i]);
  ntuple->SetBranchAddress(var_sw.c_str(),&val_sw);

  TStopwatch timer;
  timer.Start(true);

  uint nPosWei = 0;
  uint nNegWei = 0;
  uint nZerWei = 0;

  Long64_t nev = ntuple->GetEntries();
  for (int i=0; i<nev; ++i) {
    ntuple->GetEntry(i);
    val_sw_pos.push_back(val_sw);
    is_sw_mod.push_back(false);
    if (val_sw>0) nPosWei++;
    else if (val_sw<0) nNegWei++;
    else nZerWei++;
    for (uint iVar=0; iVar<vars.size(); ++iVar)
      if (vals[iVar]==vals[iVar]) {
	xcnt[iVar]++;
	xsum[iVar] += vals[iVar];
	x2sum[iVar] += vals[iVar]*vals[iVar];
      } //else
	//cout<<"NaN "<<vars[iVar]<<" value spotted: ev "<<i<<endl;
  }

  timer.Stop();
  cout<<"Vector and sums filling time: "<<timer.RealTime()<<endl;
  cout<<"Original number of weights: pos: "<<nPosWei<<" zero: "<<nZerWei<<" neg: "<<nNegWei<<endl;
  timer.Start(true);

  for (uint i=0; i<vars.size(); ++i)
    rms.push_back( x2sum[i]/xcnt[i] - xsum[i]*xsum[i]/xcnt[i]/xcnt[i] );

  for (int subs=0; subs<8192; ++subs) {
    string filename = Form("/eos/user/a/aboletti/BdToKstarMuMu/sWeightsConversion/%i_data_beforsel_o%i_v2.root",year,subs);
    if ( access( filename.c_str(), F_OK ) == -1 ) continue;
    auto fin_sub = new TFile(filename.c_str(),"READ");
    if (!fin_sub) continue;	// TEMP
    auto tin_sub = (TTree*)fin_sub->Get(Form("ntuple_%i",subs));
    if (!tin_sub) continue;	// TEMP
    double val_sw2 = 0;
    int index = 0;
    tin_sub->SetBranchAddress("nsig_sw_pos",&val_sw2);
    tin_sub->SetBranchAddress("index",&index);
    // tin_sub->SetBranchAddress("index",(double*)&index);
    for (int i=0; i<tin_sub->GetEntries(); ++i) {
      tin_sub->GetEntry(i);
      // ntuple->GetEntry(index);
      // if (val_sw_pos[index]>0 && val_sw2!=val_sw_pos[index])
      // 	cout<<i<<"\t"<<index<<"\t"<<val_sw<<"\t"<<val_sw_pos[index]<<"\t"<<val_sw2<<endl;
      if (is_sw_mod[index]) {
	ntuple->GetEntry(index);
	cout<<"ERROR! Weight of event "<<index<<" ("<<i<<"th in subs "<<subs<<") already modified: "<<val_sw<<" -> "<<val_sw_pos[index]<<" -> "<<val_sw2<<endl;
	continue;
      }
      if (fabs(val_sw2) - fabs(val_sw_pos[index]) > 1e-6 ) {
	cout<<"ERROR! Weight of event "<<index<<" ("<<i<<"th in subs "<<subs<<") did not decrease: "<<val_sw_pos[index]<<" -> "<<val_sw2<<" diff = "<<fabs(val_sw2)-fabs(val_sw_pos[index])<<endl;
	continue;
      }
      val_sw_pos[index] = val_sw2;
      is_sw_mod[index] = true;
    }
    fin_sub->Close();
  }

  timer.Stop();
  cout<<"Weight merging time: "<<timer.RealTime()<<endl;

  nPosWei = nNegWei = nZerWei = 0;
  for (int i=0; i<nev; ++i) {
    if (val_sw_pos[i]>0) nPosWei++;
    else if (val_sw_pos[i]<0) nNegWei++;
    else nZerWei++;
  }
  cout<<"Processed number of weights: pos: "<<nPosWei<<" zero: "<<nZerWei<<" neg: "<<nNegWei<<endl;
  timer.Start(true);

  vector<uint> map;
  for (int i=0; i<nev; ++i) {
    if (i%nsubs != subs) continue;
    ntuple->GetEntry(i);
    if (val_sw<=0) continue;
    bool isProblematic = false;
    for (uint iVar=0; iVar<vars.size(); ++iVar)
      if (vals[iVar]!=vals[iVar]) {
	isProblematic = true;
	break;
      }
    if (isProblematic) continue;
    map.push_back(i);
  }

  timer.Stop();
  cout<<"Map filling time: "<<timer.RealTime()<<endl;
  timer.Start(true);

  for (int i=0; i<nev; ++i) {

    // if (i%(nev/10000)==0) cout<<i*100000/nev<<"%"<<endl;
    // if (i%(nev/10000)==0) {
    //   cout<<i*100000/nev<<"% time: "<<timer.RealTime()<<endl;
    //   timer.Start(false);
    // }
    if (i%(nev/10)==0) cout<<i*100/nev<<"%"<<endl;

    if (i%nsubs != subs) continue;
    ntuple->GetEntry(i);
    if (val_sw_pos[i]>=0)
      continue;

    bool isProblematic = false;
    for (uint iVar=0; iVar<vars.size(); ++iVar)
      if (vals[iVar]!=vals[iVar]) {
	isProblematic = true;
	break;
      }
    if (isProblematic) continue;

    cout<<i<<"\t"<<map.size()<<endl;

    vals_ori = vals;

    float leftover = processList(map, ntuple, &vals, &vals_ori);

    if (leftover>0) cout<<"ERROR, insufficient numer of events!"<<endl;
    val_sw_pos[i] = -1*leftover;

  }

  timer.Stop();
  cout<<"Main loop time: "<<timer.RealTime()<<endl;

  auto fout = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/sWeightsConversion/%i_data_beforsel_%imod%i.root",year,subs,nsubs),"RECREATE");
  fout->cd();
  auto tout = new TTree(Form("ntuple_%imod%i",subs,nsubs),"ntuple");
  double nsig_sw_pos = 0.;
  tout->Branch("nsig_sw_pos", &nsig_sw_pos, "nsig_sw_pos/D");
  for (int i=0; i<nev; ++i) {
    nsig_sw_pos = val_sw_pos[i];
    tout->Fill();
  }
  tout->Write();
  fout->Close();
  
  // ntu_out->Write("", TObject::kOverwrite);
  // delete &map;
  fin->Close();
  // delete fin;

}

int main(int argc, char** argv)
{

  int year = 2016;
  uint subs = 0;
  uint nsubs = 1;

  if ( argc > 1 ) subs = atoi(argv[1]);
  if ( argc > 2 ) nsubs = atoi(argv[2]);
  if ( argc > 3 ) year = atoi(argv[3]);

  if (subs >= nsubs) {
    cout<<"Invalid configuration: "<<subs<<" % "<<nsubs<<endl;
    return 1;
  }

  merge_convert_sweights(year,subs,nsubs);

  return 0;

}
