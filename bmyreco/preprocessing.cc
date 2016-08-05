#include "BRawMasterHeader.h"
#include "BExtractedImpulse.h"
#include "BExtractedImpulseTel.h"
#include "BGeom.h"
#include "BGeomApply.h"
#include "BGeomTel2015.h"
#include "BGeomTel2014.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TMath.h"
#include <iostream>
#include <vector>

const float clight=3e8*pow(1.33,-1);


std::vector<float> getHitLocation(BGeomTel2015 geom, int chanID)
{
  std::vector<float> point;
  point.push_back((geom.At(chanID))->GetX());
  point.push_back((geom.At(chanID))->GetY());
  point.push_back((geom.At(chanID))->GetZ());
  return point;
}

//return time of propagation from arbitrary point A on trajectory to coordinates M
 
float getTimeEstimate(std::vector<float> A, std::vector<float> s, std::vector<float> M)
{
  std::vector<float> AM;
  AM.push_back(M[0]-A[0]);
  AM.push_back(M[1]-A[1]);
  AM.push_back(M[2]-A[2]);
  float modAM=sqrt(pow(AM[0],2)+pow(AM[1],2)+pow(AM[2],2));
  float modS=sqrt(pow(s[0],2)+pow(s[1],2)+pow(s[2],2));
  //  std::cout<<"AM: "<<modAM<<"  s: "<<modS<<std::endl;

  float cosAlpha=(AM[0]*s[0]+AM[1]*s[1]+AM[2]*s[2])/(modAM*modS);
  float sinAlpha=sqrt(1-pow(cosAlpha,2));

  //  std::cout<<"cos, sin: "<<cosAlpha<<"   "<<sinAlpha<<std::endl;

  float dist=modAM*(cosAlpha - 2.47*sinAlpha);
  
  //  std::cout<<"dist:  "<<dist<<std::endl;

  float time=1e9*dist/(clight);
    
  return time;
}

//filter hits based on amplitude
//order them in time
std::vector<int> filterHits(BExtractedImpulseTel* pulses)
{
  int impulse_n=pulses->GetNimpulse();
  std::vector<int> selectByAmpl_id; 
  for (int j=0; j< impulse_n; j++){
    //Geometry for 120 modules is available. 118 and 119 are noisy. Take hits with integrated charge > 200
    if (pulses->GetNch(j)<118&&pulses->GetQ(j)>200) selectByAmpl_id.push_back(j); 
  }

  std::vector<int> orderT_id;
  orderT_id.resize(selectByAmpl_id.size());

  for (Int_t j = 0; j < selectByAmpl_id.size(); j++){
    int nBefore=0;
    for (int k=0; k<selectByAmpl_id.size(); k++) {
      if (pulses->GetT(selectByAmpl_id[j])>pulses->GetT(selectByAmpl_id[k])) nBefore++;
    }
    orderT_id[nBefore]=selectByAmpl_id[j];
  }
  return orderT_id;
}

std::vector<int> makePrelimTriplet(BGeomTel2015 geom, BExtractedImpulseTel* pulses, std::vector<int> orderT_id)
{

  //  std::cout<<"clight: "<<clight<<std::endl;
  //count pairs with deltaT< c*deltaR
  //make "triplets" from such pairs
  
      //count good pairs
  int npair=0;
  
  std::vector<int> triplet_id;

  std::vector<bool> isInTriplet_chan;
    isInTriplet_chan.resize(120);
    for (int p=0; p<120; p++){
      isInTriplet_chan[p]=false;
    }      

  for (int j=0; j< orderT_id.size(); j++){
    //      std::cout<<j<<"   "<<pulses->GetNch(orderT_id[j])<<"   "<<pulses->GetQ(orderT_id[j])<<"   "<<pulses->GetT(orderT_id[j])<<std::endl;
    //coordinates vector
    int chanID_first=pulses->GetNch(orderT_id[j]);
    std::vector <float> A;
    A=getHitLocation(geom, chanID_first);
    
    //plot pairwise speed
    for (int k=j+1; k < orderT_id.size(); k++){
      int chanID_sec=pulses->GetNch(orderT_id[k]);
      if (chanID_sec==chanID_first) continue;
      std::vector <float> B;
      B=getHitLocation(geom, chanID_sec);
      
      float dR=sqrt(pow(A[0]-B[0],2)+pow(A[1]-B[1],2)+pow(A[2]-B[2],2));
      float dT=-pulses->GetT(orderT_id[j])+pulses->GetT(orderT_id[k]);
      //      if (dR>100) hPairwiseSpeed->Fill(fabs((clight*dT*2*1e-9)/dR),1);
      //      hDR->Fill(dR,1);
      //      hDT->Fill(fabs(dT),1);
      
      //count candidate pairs: signal is always propagating faster than the speed of light
      //fill vector with modules to form triplet or more
      if (fabs(clight*dT*2*1e-9)/dR<1.1) {
	npair++;
	//try to find another fitting hit and form triplet
	if (!isInTriplet_chan[chanID_first]) {
	  triplet_id.push_back(orderT_id[j]);
	  isInTriplet_chan[chanID_first]=true;
	}
	if (!isInTriplet_chan[chanID_sec]) {
	  triplet_id.push_back(orderT_id[k]);
	  isInTriplet_chan[chanID_sec]=true;
	}
      }
      
    }
  }
  //  std::cout<<"size of triplet: "<<triplet_id.size()<<std::endl;
  return triplet_id;
  
}

std::vector<int> makeQualityTriplet(BGeomTel2015 geom, BExtractedImpulseTel* pulses, std::vector<int> triplet_id, int vtx_id)
{
  //check direction of hits from given vertex, selected from preliminary triplet vector
  //ensure that all pairs satisfy the timing cut

  //preliminary quality triplet
  std::vector<int> tripletPreQual_id;
  //final quality triplet
  std::vector<int> tripletQual_id;
  tripletQual_id.resize(0);
  //hNPairs->Fill(npair,1);
  
  if (triplet_id.size()>2) {
    //      vert++;
    //    std::cout<<"NEXT EVENT       size of triplet: "<<triplet_id.size()<<std::endl;
    
    tripletPreQual_id.push_back(triplet_id[vtx_id]);
    int chanID_vtx=pulses->GetNch(triplet_id[vtx_id]);
    std::vector <float> vtx;
    vtx=getHitLocation(geom,chanID_vtx);
    
    float dRprev=0.1;
    int nTaken=0;
    std::vector<float> dirPrev;
    
    for (int k=vtx_id; k<triplet_id.size(); k++){
      int chanID_test=pulses->GetNch(triplet_id[k]);
      std::vector <float> test;
      test=getHitLocation(geom,chanID_test);
      
      std::vector<float> direction;
      direction.push_back(test[0]-vtx[0]);
      direction.push_back(test[1]-vtx[1]);
      direction.push_back(test[2]-vtx[2]);
      float modDir=sqrt(pow(direction[0],2)+pow(direction[1],2)+pow(direction[2],2));
      //      std::cout<<"direction: "<<direction[0]<<"  "<<direction[1]<<"  "<<direction[2]<<std::endl;
      
      float cos=1;
      if (nTaken>0) {
	float modDirPrev=sqrt(pow(dirPrev[0],2)+pow(dirPrev[1],2)+pow(dirPrev[2],2));
	
	cos=(direction[0]*dirPrev[0]+direction[1]*dirPrev[1]+direction[2]*dirPrev[2])/(modDir*modDirPrev);
	//	std::cout<<"chanID_vtx: "<<chanID_vtx<<"   chanID_test: "<<chanID_test<<"   modDir: "<<modDir<<"   modDirPrev: "<<modDirPrev<<"   cos: "<<cos<<std::endl;
      }
      
      //      std::cout<<"chanID_vtx: "<<chanID_vtx<<"   chanID_test: "<<chanID_test<<"   modDir: "<<modDir<<"   dRprev: "<<dRprev<<"   cos: "<<cos<<std::endl;
      if (modDir>0&&modDir>=dRprev&&dRprev>0&&cos>=-0.1) {
	nTaken++;
	//	  std::cout<<"take into triplet"<<std::endl;
	tripletPreQual_id.push_back(triplet_id[k]);
	dirPrev=direction;
	dRprev=modDir;
      }
      
    }

    //check if all pairwise combinations within the triplet satisfy timing cut
    //vector of bad channels
    std::vector<bool> isBad;
    isBad.resize(tripletPreQual_id.size());
    for (int k=0; k<isBad.size(); k++){
      isBad[k]=true;
    }
    
    if (tripletPreQual_id.size()>2){
      for (int k=0; k<tripletPreQual_id.size(); k++){
	int chanID_A=pulses->GetNch(triplet_id[k]);
	std::vector <float> A=getHitLocation(geom, chanID_A);
	for (int p=k+1; p<tripletPreQual_id.size(); p++){
	  int chanID_B=pulses->GetNch(triplet_id[p]);
	  std::vector <float> B=getHitLocation(geom, chanID_B);
	  float dR=sqrt(pow(A[0]-B[0],2)+pow(A[1]-B[1],2)+pow(A[2]-B[2],2));
	  float dT=-pulses->GetT(tripletPreQual_id[p])+pulses->GetT(tripletPreQual_id[k]);
	  if (clight*dT*2*1e-9/dR<1.1) isBad[p]=false; 
	}
      }
    }
    
    //final quality triplet
    for (int k=0; k<tripletPreQual_id.size(); k++){
      if (!isBad[k]) tripletQual_id.push_back(tripletPreQual_id[k]);
    }
    
    if (tripletQual_id.size()<3) tripletQual_id.resize(0);
    
  }
  //  std::cout<<"bloblo  "<<tripletQual_id.size()<<std::endl;
  return tripletQual_id;
}

//maximization of LL
//likelyhood* P(x1, y1, z1, x2, y2, z2, T1, T2)*P(x2, y2, z2, x3, y3, z3, T2, T3)*P(...)...
//Ti is computed propagation time
//P is gaussian with ~resolution width or any other function peaked at measured time

//function to calculate likelyhood value 
//takes geometry, array of impulses, vector of hit id, zero point and direction of trajectory

float getLvalue(BGeomTel2015 geom, BExtractedImpulseTel* pulses, std::vector<int> tripletQual_id, std::vector<float> zero, std::vector<float> direction)
{
  //array of time of impulse arrival as calculated from given trajectory 
  std::vector<float> timeEstimate;
  //array of time of impulse arrival as measured
  std::vector<float> timeMeasured;  
  //calculate time of propagation to each channel
  for (int i=0; i<tripletQual_id.size(); i++){
    //get location of hit:
    int chanID=pulses->GetNch(tripletQual_id[i]);
    std::vector<float> M = getHitLocation(geom, chanID);
    //fill time vectors
    timeEstimate.push_back(getTimeEstimate(zero, direction, M));
    timeMeasured.push_back(pulses->GetT(tripletQual_id[i]));
  }
  
  //return value
  float L=1;

  //define "resolution function"
  for (int i=1; i<tripletQual_id.size(); i++) {
    float sigma=10;
    float mean=2*(timeMeasured[i]-timeMeasured[i-1]);
    float estimate=(timeEstimate[i]-timeEstimate[i-1]);
    //      std::cout<<"mean: "<<mean<<"   estimate: "<<estimate<<std::endl;   
    float P=pow(sqrt(2*TMath::Pi())*sigma,-1)*exp(-pow(estimate-mean,2)/(2*pow(sigma,2)));
    L=L*P;
  }
  
  return L;
}


int main(int argc, char *argv[])
{
  TApplication theApp("App", &argc, argv);
  
  int nThreads=1;
  
  //deltaT*clight/R for pairwise module combinations
  TH1F* hPairwiseSpeed=new TH1F("hPairwiseSpeed","hPairwiseSpeed",1000,-1,100);
  
  TH1F* hPS_maxQ=new TH1F("hPS_maxQ","hPS_maxQ",10000,-1,1000);

  TH1F* hDR=new TH1F("hDR","hDR",1000,0,1000);
  TH1F* hDT=new TH1F("hDT","hDT",1000,0,1000);

  TH1F* hNPairs=new TH1F("hNPairs","hNPairs",100,0,100);

  TH1F* hTripletSize=new TH1F("hTripletSize","hTripletSize",100,0,100);

  TH1F* hIsSameSection=new TH1F("hIsSameSection","hIsSameSection",2,0,2);
  TH1F* hIsSameSectionQual=new TH1F("hIsSameSectionQual","hIsSameSectionQual",2,0,2);
    
  //get geometry
  
  //  geom.SetGeometry("BGeomTel2015");
 
  BGeomTel2015 geom;
  for (int i=0; i<192; i++)
      std::cout<<" channel: "<<i<<"       "<<(geom.At(i))->GetX()<<"  "
	       <<(geom.At(i))->GetY()<<"  "<<(geom.At(i))->GetZ()<<std::endl;
  
  TFile* fIN = new TFile("/home/local1/work/baikal/bars/BARS/tags/0.2.0/MYAN/e0461_copy.root");
  //TFile* fIN = new TFile("/home/local1/work/baikal/bars/BARS/tags/0.2.0/macros/OUT/e0556_01-31_orig.root");

  //std::cout<<"yoprst"<<std::endl;

  TTree* event_tree = (TTree*)fIN->Get("Events");
  TBranch* pulses_branch=event_tree->GetBranch("BExtractedImpulseTel.");
  
  TBranch* master_branch=event_tree->GetBranch("BRawMasterHeader.");

  std::cout<<"master entries: "<<master_branch->GetEntries()<<std::endl;

  std::cout<<"pointer 0: "<<event_tree<<std::endl;
  
  std::cout<<"number of baskets: "<<pulses_branch->GetMaxBaskets()<<std::endl;

  //Set branch address
  BExtractedImpulseTel* pulses = new BExtractedImpulseTel();
  pulses_branch->SetAddress(&pulses);

  //set header master header branch
  BRawMasterHeader* masters = new BRawMasterHeader();
  master_branch->SetAddress(&masters);

  int event_n=pulses_branch->GetEntries();
  
  int event_n_m=master_branch->GetEntries();

  std::cout<<"number of events: "<<event_n<<"   "<<event_n_m<<std::endl;
  
  int vert=0;

  TFile* fOUT=new TFile("fOUT.root","RECREATE");
  
  TH1F* hChanID=new TH1F("hChanID","hChanID",200,0,200);  
  TH1F* hAmpl80=new TH1F("hAmpl80","hAmpl80",10000,0,10000);

  for (Int_t i = 0; i< event_n; i++){
    pulses_branch->GetEntry(i);
    if (i%100000==0) std::cout<<" event "<<i<<std::endl;
    Int_t impulse_n = pulses->GetNimpulse();
    int nImp=0;

    for (int j=0; j<impulse_n; j++){
      int iChan=pulses->GetNch(j);
      int iSec=(geom.At(iChan))->GetSecNum();
      int iSecRef=0;
      if (j==0) iSecRef=(geom.At(iChan))->GetSecNum();
      if (iSec==iSecRef) hIsSameSection->Fill(1.0,1);
      else hIsSameSection->Fill(0.0,1);
    }

    std::vector<int> orderT_id=filterHits(pulses);
    
    std::vector<int> triplet_id=makePrelimTriplet(geom, pulses, orderT_id);

    std::vector<int> tripletQual_id;

    if (triplet_id.size()>1){
      for (int k=0; k<triplet_id.size()-2; k++){
	 std::vector<int> tripletQualBuf_id=makeQualityTriplet(geom, pulses, triplet_id, k);
	if(tripletQualBuf_id.size()>0) 
	  {
	    //	    std::cout<<"event "<<i<<" got quality triplet"<<std::endl;
	    hTripletSize->Fill(tripletQualBuf_id.size(),1);
	    tripletQual_id=tripletQualBuf_id;
	  }
	else continue;
      }
    }
    else continue;

    if(tripletQual_id.size()<3) continue;

    for (int j=0; j<tripletQual_id.size(); j++){
      int iChan=pulses->GetNch(tripletQual_id[j]);
      int iSec=(geom.At(iChan))->GetSecNum();
      int iSecRef=0;
      if (j==0) iSecRef=(geom.At(iChan))->GetSecNum();
      if (iSec==iSecRef) hIsSameSectionQual->Fill(1.0,1);
      else hIsSameSectionQual->Fill(0.0,1);
    }
    

    //perform fit
    //define initial trajectory:
    //take r13 as initial vector
    //r1-r13 as initial 0 point

    //initial point vector A
    std::vector<float> A;
    //direction vector s
    std::vector<float> s;
    
    //    std::cout<<"IDs:  "<<tripletQual_id[0]<<"  "<<tripletQual_id[1]<<"  "<<tripletQual_id[2]<<std::endl;

    //M1
    int chanID_M1=pulses->GetNch(tripletQual_id[0]);
    std::vector<float> M1;
    M1=getHitLocation(geom, chanID_M1);
    //M2
    int chanID_M2=pulses->GetNch(tripletQual_id[1]);
    std::vector<float> M2;
    M2=getHitLocation(geom, chanID_M2);
    //M3
    int chanID_M3=pulses->GetNch(tripletQual_id[2]);
    std::vector<float> M3;
    M3=getHitLocation(geom, chanID_M3);
     
    //construct A
    A.push_back(M1[0]-(M3[0]-M1[0]));
    A.push_back(M1[1]-(M3[1]-M1[1]));
    A.push_back(M1[2]-(M3[2]-M1[2]));

    //construct direction
    s.push_back(M3[0]-M1[0]);
    s.push_back(M3[1]-M1[1]);
    s.push_back(M3[2]-M1[2]);

    //calculate T_i for a given trajectory
    //test: get time for M2:
    float T2 = getTimeEstimate(A,s,M2);
    
    float L=getLvalue(geom, pulses, tripletQual_id, A, s);
    //    std::cout<<"NEXT EVENT , propagation to M2: "<<T2<<"    likelyhood: "<<L<<std::endl;
    
  }
  
  std::cout<<vert<<std::endl;
  fOUT->cd();
  hChanID->Write();
  hAmpl80->Write();
  hPairwiseSpeed->Write();
  hPS_maxQ->Write();
  hDR->Write();
  hDT->Write();
  hNPairs->Write();
  hTripletSize->Write();
  hIsSameSection->Write();
  hIsSameSectionQual->Write();
  fOUT->Close();

  fIN->Close();

  return 0;
}
