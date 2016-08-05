#include "BMyRecoMC.h"
#include "BGeomTel.h"
#include "BEvent.h"
#include "BExtractedImpulse.h"
#include "BExtractedImpulseTel.h"
#include "BJoinExtractedImpulseTel.h"
#include "BEventMask.h"
#include "BChannelMask.h"
#include "BExtractedHeader.h"
#include "BMCEvent.h"

//#include "BMyRecParam.h"

#include "MParList.h"
#include "MLog.h"
#include "MLogManip.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TMath.h"

#include <math.h>
#include <vector>

ClassImp(BMyRecoMC);


//Add functions from BMyReco

BMyRecoMC::BMyRecoMC(string fname)
{
  iEvent=0;
   
  fOUT=new TFile(fname.c_str(),"RECREATE");
  hMuonN=new TH1F("hMuonN","hMuonN",200,0,200);
  hTotalMuonN=new TH1F("hTotalMuonN","hTotalMuonN",100,0,100);
  hMuCutN=new TH1F("hMuCutN","hMuCutN",100,0,100);
  hMuNotNull=new TH1F("hMuNotNull","hMuNotNull",100,0,100);
  
  hMuonE=new TH1F("hMuonE","hMuonE",1000,0,1000);

  hThetaPrimary=new TH1F("hThetaPrimary","hThetaPrimary",1000,0,1000);
  hNHits_ThetaPrimary=new TH2F("hNHits_ThetaPrimary","hNHits_ThetaPrimary",100,0,100,100,0,100);
  hNHits_EnergyPrimary=new TH2F("hNHits_EnergyPrimary","hNHits_EnergyPrimary",1000,0,100000,100,0,100);


  hChannelN=new TH1F("hChannelN","hChannelN",1000,0,1000);

  hChannelN_bevt=new TH1F("hChannelN_bevt","hChannelN_bevt",1000,0,1000);
  hPulseN_bevt=new TH1F("hPulseN_bevt","hPulseN_bevt",1000,0,1000);
  
  hMapOfHits1=new TH2F("hMapOfHits1","hMapOfHits1",21,-0.5,20.5,25,-0.5,24.5);
  hMapOfHits3=new TH2F("hMapOfHits3","hMapOfHits3",21,-0.5,20.5,25,-0.5,24.5);
  hGeom=new TH2F("hGeom","hGeom",2000,0.5,20.5,10000,0,1000);

  clight=3e8*pow(1.33,-1);
  cmu=3e8;
}

BMyRecoMC::~BMyRecoMC()
{
}

Int_t BMyRecoMC::PreProcess(MParList * pList)
{
  std::cout<<"WE ARE IN BMyRecoMC::PreProcess"<<std::endl;
  
  fMCEvent=(BMCEvent*)pList->FindObject("BMCEvent","BMCEvent");

  fEvent=(BEvent*)pList->FindObject("BEvent","BEvent");
  
  fGeomTel = (BGeomTel*)pList->FindObject(AddSerialNumber("BGeomTel"));

  std::cout<<"N OM: "<<fGeomTel->GetNumOMs()<<std::endl;
  std::cout<<"section number of 100th module: "<<fGeomTel->At(100)->GetSecNum()<<std::endl;
  std::cout<<"coordinates of 100th module: "<<fGeomTel->At(100)->GetX()<<"  "<<fGeomTel->At(100)->GetY()<<"  "<<fGeomTel->At(100)->GetZ()<<std::endl;

  return kTRUE;
}



Int_t BMyRecoMC::Process()
{
  iEvent++;
  // std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<WE ARE IN BMyReco::Process"<<std::endl;
  if (iEvent%10000==0) std::cout<<"eventNumber: "<<iEvent<<std::endl;

  int nMuon=fMCEvent->GetMuonsN();
  hMuonN->Fill(nMuon,1);
  int nTotalMuon=fMCEvent->GetResponseMuonsN();
  hTotalMuonN->Fill(nTotalMuon,1);
  int nMu=0;
  int nMuNotNull=0;
  float muE=0;
  for (int i=0; i<nTotalMuon; i++){
    nMuNotNull++;
    muE=(fMCEvent->GetTrack(i))->GetMuonEnergy();
    if (muE>200) nMu++;
    //hNHits_ThetaPrimary->Fill(180-fMCEvent->GetPrimaryParticlePolar(), fMCEvent->GetChannelN(), 1);
    //hNHits_EnergyMuon->Fill((fMCEvent->GetTrack(i))->GetMuonEnergy(), 1);
  }
  
  hThetaPrimary->Fill(fMCEvent->GetPrimaryParticlePolar(),1);

   hNHits_ThetaPrimary->Fill(180-fMCEvent->GetPrimaryParticlePolar(), fMCEvent->GetChannelN(), 1);
   //std::cout<<"energy primary: "<<fMCEvent->GetSumEnergyBundleSurf()<<std::endl;
   if (180-fMCEvent->GetPrimaryParticlePolar() < 20) hNHits_EnergyPrimary->Fill(fMCEvent->GetSumEnergyBundleSurf(), fMCEvent->GetChannelN(), 1);

  hMuonE->Fill(muE,1);
  hMuCutN->Fill(nMu,1);
  hMuNotNull->Fill(nMuNotNull,1);

  int nChannels=fMCEvent->GetChannelN();
  if (nTotalMuon==1) hChannelN->Fill(nChannels,1);

  int nPulses=fEvent->GetTotImpulses();

  std::vector<int> orderT_id=filterHits(1);

  

  

  //if(vOrderedHits.size()>0) std::cout<<vOrderedHits.size()<<std::endl;
  //  for (int i =0; i< nPulses; i++){
  //    fEvent->GetImpulse(i)->GetNimpulse();
  //  }
 
  hChannelN_bevt->Fill(fEvent->GetChannelN(),1);
  hPulseN_bevt->Fill(fEvent->GetTotImpulses(),1);
  //std::cout<<fEvent->GetChannelsN()<<"    "<<fEvent->GetTotImpulses()<<std::endl;
  
 

  // for (int i=0; i<nChannels; i++){
    
  //  }

  // std::cout<<"number of nucleons in primary partile: "<<fMCEvent->GetNucleonN()<<std::endl;

  return kTRUE;
}

Int_t BMyRecoMC::PostProcess()
{
  std::cout<<"WE ARE IN BMyReco::PostProcess"<<std::endl;
  fOUT->cd();
  hMuonN->Write();
  hTotalMuonN->Scale(pow(hTotalMuonN->GetEntries(),-1));
  hTotalMuonN->Write();
  hMuCutN->Write();
  hMuNotNull->Write();
  hMuonE->Write();

  hThetaPrimary->Write();
  hNHits_ThetaPrimary->Write();
  hNHits_EnergyPrimary->Write();

  hChannelN->Write();

  hChannelN_bevt->Write();
  hPulseN_bevt->Write();

  hMapOfHits1->Write();
  hMapOfHits3->Write();
  hGeom->Write();
  return kTRUE;
}

/*
std::vector<float> BMyReco::getHitLocation(int hitID)
{
  GetHitChannel(Int_t index)
  int chanID=
    fExtractedImpulse->GetNch(hitID);
  std::vector<float> point;
  point.push_back((fGeomTel->At(chanID))->GetX());
  point.push_back((fGeomTel->At(chanID))->GetY());
  point.push_back((fGeomTel->At(chanID))->GetZ());
  return point;
}
*/

//filter hits based on amplitude
//order them in time
std::vector<int> BMyRecoMC::filterHits(float minAmpl)
{
  int impulse_n=fEvent->GetTotImpulses();
  //test geometry
  /*
  std::cout<<"N OM: "<<fGeomTel->GetNumOMs()<<std::endl;
  std::cout<<"section number of 100th module: "<<fGeomTel->At(100)->GetSecNum()<<std::endl;
  std::cout<<"coordinates of 100th module: "<<fGeomTel->At(100)->GetX()<<"  "<<fGeomTel->At(100)->GetY()<<"  "<<fGeomTel->At(100)->GetZ()<<std::endl;
  */
  //if (impulse_n>1) std::cout<<"imp: "<<impulse_n<<std::endl;
  std::vector<int> selectByAmpl_id; 
  for (int j=0; j<impulse_n; j++){
    //118 and 119 are noisy. Take hits with integrated charge > 200
    //fill simple histograms
    //hCharge->Fill(fEvent->GetImpulse(j)->GetAmplitude(),1);
    int nch=fEvent->GetImpulse(j)->GetChannelID();
    //    std::cout<<"yobls: "<<nch<<"   "<<fGeomTel->At(100)->GetZ()<<std::endl;
    if (fEvent->GetImpulse(j)->GetAmplitude()>1) hMapOfHits1->Fill(floor((nch)/24),(nch)%24+1);
    if (fEvent->GetImpulse(j)->GetAmplitude()>3) hMapOfHits3->Fill(floor((nch)/24),(nch)%24+1);
    hGeom->Fill(floor((nch)/24),-(fGeomTel->At(nch)->GetZ()),1);
    //fGeomTel->At(nch)->GetZ(),1);
    //(fGeomTel->At(chanID))->GetZ()
    //std::cout<<"ampl:   "<<fEvent->GetImpulse(j)->GetAmplitude()<<std::endl;
    if (nch!=118&&nch!=119&&fEvent->GetImpulse(j)->GetAmplitude()>minAmpl) selectByAmpl_id.push_back(j); 
  }

  std::vector<int> orderT_id;
  orderT_id.resize(selectByAmpl_id.size());

  for (Int_t j = 0; j < selectByAmpl_id.size(); j++){
    //hSecNum->Fill(floor(fExtractedImpulse->GetNch(selectByAmpl_id[j])/24),1);
    int nBefore=0;
    for (int k=0; k<selectByAmpl_id.size(); k++) {
      if (fEvent->GetImpulse(selectByAmpl_id[j])->GetTime()>fEvent->GetImpulse(selectByAmpl_id[k])->GetTime()) nBefore++;
    }
    orderT_id[nBefore]=selectByAmpl_id[j];
  }
  return orderT_id;
}

std::vector<float> BMyRecoMC::getHitLocation(int hitID)
{
  int chanID=fEvent->GetImpulse(hitID)->GetChannelID();
  std::vector<float> point;
  point.push_back((fGeomTel->At(chanID))->GetX());
  point.push_back((fGeomTel->At(chanID))->GetY());
  point.push_back((fGeomTel->At(chanID))->GetZ());
  return point;
}

float BMyRecoMC::getDeltaT_ns(int hit1_id, int hit2_id)
{
  int tim1=fEvent->GetImpulse(hit1_id)->GetTime();
  int tim2=fEvent->GetImpulse(hit2_id)->GetTime();
  float deltaT;
  deltaT=5*(tim2-tim1); //one unit is 5ns
  return deltaT;
}

float BMyRecoMC::getDeltaR_m(int hit1_id, int hit2_id)
{
  std::vector<float> Hit1=getHitLocation(hit1_id);
  
  std::vector<float> Hit2=getHitLocation(hit2_id);
  
  float deltaR=sqrt(pow(Hit1[0]-Hit2[0],2)+pow(Hit1[1]-Hit2[1],2)+pow(Hit1[2]-Hit2[2],2));
  
  return deltaR;
}

std::vector<float> BMyRecoMC::getDirection(int hit1_id, int hit2_id)
{
  std::vector <float> sec=getHitLocation(hit2_id);
  
  //previous hit on cluster
  std::vector<float> first=getHitLocation(hit1_id);
    
  //direction from previous hit
  std::vector<float> dir;
  dir.push_back(sec[0]-first[0]);
  dir.push_back(sec[1]-first[1]);
  dir.push_back(sec[2]-first[2]);

  return dir;
}


