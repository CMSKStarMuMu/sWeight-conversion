#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>

#include <iostream>

using namespace std;

vector<double> rms(0);

vector<double> val_sw_pos(0);
vector<bool> is_sw_mod(0);

float getDist(vector<double>* v1, vector<double>* v2)
{
  float ret = 0;
  for (uint i=0; i<v1->size(); ++i)
    ret += (v1->at(i)-v2->at(i))*(v1->at(i)-v2->at(i))/rms[i]/rms[i];
  return sqrt(ret);
}

float processList(vector<uint> locmap, TTree* ntuple, vector<double>* vals, vector<double>* vals_ori, double val_sw_ori)
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
      // cout<<"=====DEBUG1: "<<k<<" "<<candidates[k]<<" -> "<<val_sw_pos[candidates[k]]<<" "<<sum<<" "<<val_sw_ori<<endl;
      is_sw_mod[candidates[k]] = true;
    } else {
      // cout<<"=====DEBUG2: "<<k<<" "<<candidates[k]<<" -> "<<val_sw_pos[candidates[k]]<<" "<<sum<<" "<<val_sw_ori<<endl;
      val_sw_pos[candidates[k]] -= sum;
      is_sw_mod[candidates[k]] = true;
      sum = 0;
      // cout<<"=====DEBUG3: "<<k<<" "<<candidates[k]<<" -> "<<val_sw_pos[candidates[k]]<<" "<<sum<<" "<<val_sw_ori<<endl;
      break;
    }
  }

  // cout<<i<<" "<<val_sw_ori<<" - ";
  // for (k=0; k<candidates.size(); ++k)
  //   cout<<candidates[k]<<","<<dists[k]<<","<<val_sw_pos[candidates[k]]<<"\t";
  // cout<<endl;

  return sum;

}

void convert_sweights(int year, int subs)
{

  auto fin = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/%i_data_beforsel.root",year),"READ");
  auto ntuple = (TTree*)fin->Get("ntuple");

  // vector<string> vars = {"bVtxCL", "bLBS", "bLBSE" ,"bCosAlphaBS", "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt","kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk1DCABSE","kstTrk2DCABS","kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta","sum_isopt_04"};
  vector<string> vars = {"bVtxCL", "bLBS", "bLBSE" ,"bCosAlphaBS", "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt","kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk2DCABS","sum_isopt_04",
    "kstTrk1DCABSE","kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta"};
  vector<double> vals(vars.size(),0.);
  vector<double> vals_ori(0);
  vector<uint> mapIdxVar = {0,1,2,3,4,5,10,11,12,13,14,16,18};
  
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

  Long64_t nev = ntuple->GetEntries();
  for (int i=0; i<nev; i+=10) {
    ntuple->GetEntry(i);
    for (uint iVar=0; iVar<vars.size(); ++iVar)
      if (vals[iVar]==vals[iVar]) {
	xcnt[iVar]++;
	xsum[iVar] += vals[iVar];
	x2sum[iVar] += vals[iVar]*vals[iVar];
      } //else
	//cout<<"NaN "<<vars[iVar]<<" value spotted: ev "<<i<<endl;
  }

  timer.Stop();
  cout<<"Sums filling time: "<<timer.RealTime()<<endl;

  vector<double> mean(0);
  for (uint i=0; i<vars.size(); ++i) {
    mean.push_back( xsum[i]/xcnt[i] );
    rms.push_back( x2sum[i]/xcnt[i] - xsum[i]*xsum[i]/xcnt[i]/xcnt[i] );
    // cout<<vars[i]<<" "<<mean[i]<<" "<<rms[i]<<" "<<xsum[i]<<" "<<x2sum[i]<<endl;
  }

  timer.Start(true);

  // ulong nmap2 = pow(4,mapIdxVar.size());
  int nmap = pow(2,vars.size());
  if (subs>-1) nmap = nmap/8192;
  cout<<"Total size = "<<nmap<<endl;
  auto map = new vector<uint> [nmap];
  // vector<uint> pcnt (8192,0);
  // vector<uint> ncnt (8192,0);
  for (int i=0; i<nev; ++i) {
    ntuple->GetEntry(i);
    val_sw_pos.push_back(val_sw);
    is_sw_mod.push_back(false);
    if (val_sw<=0) continue;
    bool isProblematic = false;
    for (uint iVar=0; iVar<vars.size(); ++iVar)
      if (vals[iVar]!=vals[iVar]) {
	isProblematic = true;
	break;
      }
    if (isProblematic) continue;
    // ulong idx2 = 0;
    uint idx = 0;
    // for (uint iMapVar=0; iMapVar<mapIdxVar.size(); ++iMapVar) {
    for (uint iVar=0; iVar<vars.size(); ++iVar) {
      // uint iVar = mapIdxVar[iMapVar];
      if (vals[iVar]>mean[iVar])
	idx = 1<<iVar | idx;
      // if (vals[iVar]>mean[iVar]+rms[iVar])
      // 	idx2 = 3<<(2*iMapVar) | idx2;
      // else if (vals[iVar]>mean[iVar])
      // 	idx2 = 2<<(2*iMapVar) | idx2;
      // else if (vals[iVar]>mean[iVar]-rms[iVar])
      // 	idx2 = 1<<(2*iMapVar) | idx2;
      // cout<<" - "<<idx2;
    }
    // cout<<endl;
    // if (val_sw<=0) {ncnt[idx%8192]++; continue;}
    // pcnt[idx%8192]++;
    if (subs<0 || subs==idx%8192)
      map[subs<0 ? idx : idx/8192].push_back(i);
  }

  cout<<"Map size"<<endl;
  for (int j=0; j<nmap/8; ++j) {
    for (int i=0; i<8; ++i)
      cout<<map[j*8+i].size()<<"\t";
    cout<<endl;
  }
  // cout<<"Splitting size: neg pos"<<endl;
  // for (int i=0; i<8192; ++i)
  //   cout<<i<<"\t"<<ncnt[i]<<"\t"<<pcnt[i]<<endl;

  timer.Stop();
  cout<<"Vector and map filling time: "<<timer.RealTime()<<endl;

  // TStopwatch timer2, timer3, timer4;
  timer.Start(true);

  for (int i=0; i<nev; ++i) {

    // if (i%(nev/1000)==0) cout<<i*10000/nev<<"%"<<endl;
    // if (i%(nev/100000)==0) {
    //   cout<<i*1000000/nev<<"% time: "<<timer.RealTime()<<endl;
    //   timer.Start(false);
    // }
    // if (i%(nev/10)==0) cout<<i*100/nev<<"%"<<endl;

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

    uint idx = 0;
    for (uint iVar=0; iVar<vars.size(); ++iVar) {
      if (vals[iVar]>mean[iVar])
	idx = 1<<iVar | idx;
    }
    // cout<<idx%8192<<endl;
    if (subs>-1 && subs!=idx%8192)
      continue;
    // ulong idx2 = 0;
    // for (uint iMapVar=0; iMapVar<mapIdxVar.size(); ++iMapVar) {
    //   uint iVar = mapIdxVar[iMapVar];
    //   if (vals[iVar]>mean[iVar]+rms[iVar])
    // 	idx2 = 3<<(2*iMapVar) | idx2;
    //   else if (vals[iVar]>mean[iVar])
    // 	idx2 = 2<<(2*iMapVar) | idx2;
    //   else if (vals[iVar]>mean[iVar]-rms[iVar])
    // 	idx2 = 1<<(2*iMapVar) | idx2;
    // }
    if (subs>-1) idx = idx/8192;
    cout<<i<<"\t"<<idx<<"\t"<<map[idx].size()<<endl;

    vals_ori = vals;

    float leftover = processList(map[idx], ntuple, &vals, &vals_ori, val_sw_pos[i]);

    if (leftover>0) cout<<"ERROR, insufficient numer of events!"<<endl;
    val_sw_pos[i] = -1*leftover;
    is_sw_mod[i] = true;

    // for (uint iVar=0; iVar<vars.size(); ++iVar) {
    //   auto plusmap = map[1<<iVar ^ idx];
    //   locmap.insert(locmap.end(),plusmap.begin(),plusmap.end());
    // }

  }

  timer.Stop();
  cout<<"Main loop time: "<<timer.RealTime()<<endl;
  // cout<<"Main loop time: "<<timer.RealTime()<<" ("<<timer2.RealTime()<<" - "<<timer3.RealTime()<<" - "<<timer4.RealTime()<<")"<<endl;

  auto fout = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/sWeightsConversion/%i_data_beforsel_o%i_v2.root",year,subs),"RECREATE");
  fout->cd();
  auto tout = new TTree(Form("ntuple_%i",subs),"ntuple");
  double nsig_sw_pos = 0.;
  int iEv = 0;
  tout->Branch("index", &iEv, "index/I");
  tout->Branch("nsig_sw_pos", &nsig_sw_pos, "nsig_sw_pos/D");
  for (iEv=0; iEv<nev; ++iEv)
    if (is_sw_mod[iEv]) {
      nsig_sw_pos = val_sw_pos[iEv];
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

  int subs = -1;
  int year = 2016;

  if ( argc > 1 ) subs = atoi(argv[1]);
  if ( argc > 2 ) year = atoi(argv[2]);

  if (subs<-1 || subs>=8192) {
    cout<<"subs parameter needs to be in the range [0,8192), or -1 for full sample"<<endl;
    return 1;
  }

  convert_sweights(year,subs);

  return 0;

}
