#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>

#include <iostream>

using namespace std;

vector<double> rms2(0);

vector<double> val_sw_pos(0);

// vector<TStopwatch*> sw(5);

float getDist(vector<double>* v1, vector<double>* v2)
{
  float ret = 0;
  for (uint i=0; i<v1->size(); ++i)
    ret += (v1->at(i)-v2->at(i))*(v1->at(i)-v2->at(i))/rms2[i];
  return sqrt(ret);
}

float processList(vector<uint> locmap, TTree* ntuple, vector<double>* vals, vector<double>* vals_ori, double val_sw_ori)
{
  vector<int> candidates(0);
  vector<float> dists(0);

  // for (int j=0; j<nev; ++j) {
  float sum = 0.;
  int j = 0;
  for (uint ij=0; ij<locmap.size(); ++ij) {
    j = locmap[ij];

    if (val_sw_pos[j]<=0) continue;

    // sw[0]->Start(false);
    ntuple->GetEntry(j);
    // sw[0]->Stop();

    // sw[1]->Start(false);
    float dist = getDist(vals,vals_ori);
    // sw[1]->Stop();
    // sw[2]->Start(false);
    if ( dists.size() > 0 && dist > dists.back() && sum + val_sw_pos[j] > -1*val_sw_ori ) {
      // sw[2]->Stop();
      continue;
    }
    // sw[2]->Stop();
    // if (i==155) cout<<j<<" == "<<(dists.size()>0 ? dist-dists.back() : 999)<<" "<<sum + val_sw_pos[j] + val_sw_ori<<endl;

    // sw[3]->Start(false);
    uint k;
    for (k=0; k<candidates.size(); ++k)
      if (dist < dists[k]) break;
    dists.insert(dists.begin()+k,dist);
    candidates.insert(candidates.begin()+k,j);

    // sw[3]->Stop();
    // sw[4]->Start(false);

    sum = 0.;
    for (k=0; k<candidates.size(); ++k) {
      if (sum>-1*val_sw_ori) {
	candidates.erase(candidates.begin()+k);
	dists.erase(dists.begin()+k);
	--k;
      } else
	sum += val_sw_pos[candidates[k]];
    }
    // sw[4]->Stop();

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
    } else {
      // cout<<"=====DEBUG2: "<<k<<" "<<candidates[k]<<" -> "<<val_sw_pos[candidates[k]]<<" "<<sum<<" "<<val_sw_ori<<endl;
      val_sw_pos[candidates[k]] -= sum;
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

void convert_sweights(int year, int divi)
{

  auto fin = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/%i_data_beforsel.root",year),"READ");
  auto ntuple = (TTree*)fin->Get("ntuple");

  // vector<string> vars = {"bVtxCL", "bLBS", "bLBSE" ,"bCosAlphaBS", "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt","kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk1DCABSE","kstTrk2DCABS","kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta","sum_isopt_04"};
  vector<string> vars = {"bVtxCL", "bLBS", "bLBSE" ,"bCosAlphaBS", "bDCABS","bDCABSE","kstTrk1Pt", "kstTrk2Pt","kstTrk1Eta", "kstTrk2Eta","kstTrk1DCABS","kstTrk2DCABS","sum_isopt_04",
    "kstTrk1DCABSE","kstTrk2DCABSE","mu1Pt","mu2Pt","mu1Eta","mu2Eta"};
  vector<double> vals(vars.size(),0.);
  vector<double> vals_ori(0);
  
  double val_sw = 0.;
  Long64_t eventN = 0;

  vector<double> x2sum (vars.size(),0.);
  vector<double> xsum (vars.size(),0.);
  vector<uint> xcnt (vars.size(),0);
  
  ntuple->SetBranchStatus("*",0);
  for (uint i=0; i<vars.size(); ++i) {
    ntuple->SetBranchStatus(vars[i].c_str(),1);
    ntuple->SetBranchAddress(vars[i].c_str(),&vals[i]);
  }
  ntuple->SetBranchStatus("nsig_sw",1);
  ntuple->SetBranchStatus("eventN",1);
  ntuple->SetBranchAddress("nsig_sw",&val_sw);
  ntuple->SetBranchAddress("eventN",&eventN);

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
    rms2.push_back( x2sum[i]/xcnt[i] - xsum[i]*xsum[i]/xcnt[i]/xcnt[i] );
    // cout<<vars[i]<<" "<<mean[i]<<" "<<rms2[i]<<" "<<xsum[i]<<" "<<x2sum[i]<<endl;
  }

  timer.Start(true);

  vector<uint> map (0);
  // vector<uint> pcnt (8192,0);
  // vector<uint> ncnt (8192,0);

  for (int i=0; i<nev; ++i) {
    // Reject a fraction of the entries
    if (i%divi != 0) {
      val_sw_pos.push_back(0);
      continue;
    }

    ntuple->GetEntry(i);

    // Reject odd entries
    if (eventN%2 != 0) {
      val_sw_pos.push_back(0);
      continue;
    }

    // Reject events with NaN values
    bool isProblematic = false;
    for (uint iVar=0; iVar<vars.size(); ++iVar)
      if (vals[iVar]!=vals[iVar]) {
	isProblematic = true;
	break;
      }
    if (isProblematic) {
      val_sw_pos.push_back(0);
      continue;
    }

    val_sw_pos.push_back(val_sw);
    if (val_sw>0) map.push_back(i);
    // if (val_sw<=0) {ncnt[idx%8192]++; continue;}
    // pcnt[idx%8192]++;
  }

  // cout<<"Splitting size: neg pos"<<endl;
  // for (int i=0; i<8192; ++i)
  //   cout<<i<<"\t"<<ncnt[i]<<"\t"<<pcnt[i]<<endl;

  timer.Stop();
  cout<<"Vector and map filling time: "<<timer.RealTime()<<endl;

  ntuple->SetBranchStatus("nsig_sw",0);
  ntuple->SetBranchStatus("eventN",0);

  // for (uint iSW=0; iSW<sw.size(); ++iSW)
  //   sw[iSW] = new TStopwatch();
  timer.Start(true);

  for (int i=0; i<nev; ++i) {

    // if (i%(nev/10000)==0) cout<<i*100000/nev<<"%"<<endl;
    // if (i%(nev/100000)==0) {
    //   cout<<i*1000000/nev<<"% time: "<<timer.RealTime()<<endl;
    //   timer.Start(false);
    // }
    // if (i%(nev/10)==0) cout<<i*100/nev<<"%"<<endl;

    if (val_sw_pos[i]>=0)
      continue;

    ntuple->GetEntry(i);

    cout<<i<<"\t"<<map.size()<<endl;

    vals_ori = vals;

    float leftover = processList(map, ntuple, &vals, &vals_ori, val_sw_pos[i]);

    if (leftover>0) cout<<"ERROR, insufficient numer of events!"<<endl;
    val_sw_pos[i] = -1*leftover;

  }

  timer.Stop();
  cout<<"Main loop time: "<<timer.RealTime()<<endl;
  // for (vector<TStopwatch*>::iterator isw = sw.begin(); isw!=sw.end(); ++isw)
  //   cout<<(*isw)->RealTime()<<" ";
  // cout<<endl;

  auto fout = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/sWeightsConversion/%i_data_beforsel_posWei_div%i.root",year,divi),"RECREATE");
  fout->cd();
  auto tout = new TTree("ntuple_posWei","ntuple");
  double nsig_sw_pos = 0.;
  double runN = 0.;
  tout->Branch("runN", &runN, "runN/D");
  tout->Branch("eventN", &eventN, "eventN/L");
  tout->Branch("nsig_sw_pos", &nsig_sw_pos, "nsig_sw_pos/D");

  ntuple->SetBranchStatus("*",0);
  ntuple->SetBranchStatus("runN",1);
  ntuple->SetBranchStatus("eventN",1);
  ntuple->SetBranchAddress("runN",&runN);

  for (int i=0; i<nev; ++i) {
    ntuple->GetEntry(i);
    nsig_sw_pos = val_sw_pos[i];
    // if (nsig_sw_pos<0)
    //   cout<<"Error: negative weight "<<nsig_sw_pos<<" for event "<<i<<endl;
    tout->Fill();
  }
  tout->Write();
  fout->Close();
  
  fin->Close();

}

int main(int argc, char** argv)
{

  int divi = 1;
  int year = 2016;

  if ( argc > 1 ) divi = atoi(argv[1]);
  if ( argc > 2 ) year = atoi(argv[2]);

  if (divi<1) {
    cout<<"divi parameter needs to be larger than 0"<<endl;
    return 1;
  }

  convert_sweights(year,divi);

  return 0;

}
