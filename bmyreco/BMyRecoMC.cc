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
#include "TF1.h"
#include "TProfile.h"
#include "TFile.h"
#include "TMath.h"
#include "TRotation.h"
#include "TVector3.h"

#include <math.h>
#include <vector>

ClassImp(BMyRecoMC);


//Add functions from BMyReco

BMyRecoMC::BMyRecoMC(string fname)
{
  fName  = "BMyRecoMC";
  fTitle = "RMyRecoMC";
  
  iEvent=0;

  fNCalibTracks=0;
  
  fOUT=new TFile(fname.c_str(),"RECREATE");

  hNmuResp=new TH1F("hNmuResp","hNmuResp",50,0,50);
  hNmuTotal=new TH1F("hNmuTotal","hNmuTotal",50,0,50);
  //  hNmuResp_N5_amp3=new TH1F("hNmuResp_N5_amp3","hNmuResp_N5_amp3",50,0,50);
  
  hPolar1muMC=new TH1F("hPolar1muMC","hPolar1muMC",1000,-500,500);
  hRho1muMC=new TH1F("hRho1muMC","hRho1muMC",1000,-500,500);
  hNhit1muMC=new TH1F("hNhit1muMC","hNhit1muMC",100,0,100);
  //hNhit1muSignalMC=new TH1F("hNhit1muSignalMC","hNhit1muSignalMC",100,0,100);
  //hNhit1muSignal3peMC=new TH1F("hNhit1muSignal3peMC","hNhit1muSignal3peMC",100,0,100);
  hNdirect1muHit=new TH1F("hNdirect1muHit","hNdirect1muHit",100,0,100);
  hNshower1muHit=new TH1F("hNshower1muHit","hNshower1muHit",100,0,100);
  hNnoiseHit=new TH1F("hNnoiseHit","hNnoiseHit",100,0,100);
  hNhitPerChannel=new TH1F("hNhitPerChannel","hNhitPerChannel",193,0,193);
  hDistToTrackMC_direct=new TH1F("hDistToTrackMC_direct","hDistToTrackMC_direct",1000,0,1000);
  hDistToTrackMC_shower=new TH1F("hDistToTrackMC_shower","hDistToTrackMC_shower",1000,0,1000);
  hTimeDifference_norm=new TH1F("hTimeDifference_norm","hTimeDifference_norm",2000,-2000,2000);
  hDiffReal_DiffEst=new TH2F("hDiffReal_DiffEst","hDiffReal_DiffEst",2000,-1000,1000,2000,-1000,1000);
  hPulsesPerChannel=new TH1F("hPulsesPerChannel","hPulsesPerChannel",10,0,10);
  hYieldVsAmpl_ref=new TH1F("hYieldVsAmpl_ref","hYieldVsAmpl_ref",50,1,11);

  char tmp[100];
  for (int i=0; i<192; i++){
    sprintf(tmp,"hChanOffset_%d",i+1);
    hChanOffset[i]=new TH1F(tmp,tmp,2000,-1000,1000);
  }

  hClusterOffsets=new TH1F("hClusterOffsets","hClusterOffsets",192,0.5,192.5);
  
  hDebugTimeShowerHits=new TH2F("hDebugTimeShowerHits","hDebugTimeShowerHits",5000,0,5000,20,0,20);
  
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
  hPulseTime=new TH1F("hPulseTime","hPulseTime",10000,-10,9990);

  
  cWater=3e8*pow(1.33,-1);
  cVacuum=3e8;
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

  RunMCAnalysis();
  
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

  for (int i=0; i<nChannels; i++){
    for (int j=0; j<fMCEvent->GetHitChannel(i)->GetPulseN(); j++){
      if (fMCEvent->GetHitChannel(i)->GetPulse(j)->GetMagic()!=1) hPulseTime->Fill(fMCEvent->GetHitChannel(i)->GetPulse(j)->GetTime(),1);
    }
  }
  
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
  
  hNmuResp->Write();
  hNmuTotal->Write();
  //  hNmuResp_amp3->Write();
  hPolar1muMC->Write();
  hRho1muMC->Write();
  hNhit1muMC->Write();
  hNdirect1muHit->Write();
  hNshower1muHit->Write();
  hNnoiseHit->Write();
  hNhitPerChannel->Write();
  hDistToTrackMC_direct->Write();
  hDistToTrackMC_shower->Write();
  hTimeDifference_norm->Write();
  hDiffReal_DiffEst->Write();
  hPulsesPerChannel->Write();
  hYieldVsAmpl_ref->Write();

  fOUT->mkdir("chanOffset");
  fOUT->cd("chanOffset");

  for (int i=0; i<192; i++) {
    //fit the peak and find offset
    float meanEst=hChanOffset[i]->GetMean();
    float rmsEst=hChanOffset[i]->GetRMS();

    TF1* ffOffset=new TF1("ffOffset","gaus");
    ffOffset->SetParameter(0,hChanOffset[i]->GetEntries());
    ffOffset->SetParameter(1,hChanOffset[i]->GetMean());
    ffOffset->SetParameter(2,hChanOffset[i]->GetRMS());

    float fitMargin=1.5;
    float fitMean=0;
    float fitError=0;
    if (hChanOffset[i]->Fit(ffOffset,"Q","",meanEst-fitMargin*rmsEst,meanEst+fitMargin*rmsEst)){
      fitMean=ffOffset->GetParameter(1);
      fitError=ffOffset->GetParError(1);
    }

    std::cout<<"time offset for channel #"<<i+1<<" :    "
	     <<fitMean<<" +- "
	     <<fitError<<std::endl;
    hChanOffset[i]->Write();

    hClusterOffsets->SetBinContent(i+1,fitMean);
    hClusterOffsets->SetBinError(i+1,fitError);
  }
  
  fOUT->cd();

  
  hClusterOffsets->Write();
  hDebugTimeShowerHits->Write();
  
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
  hPulseTime->Write();

  std::cout<<"TOTAL OF "<<fNCalibTracks<<" EVENTS USED FOR CAIBRATION"<<std::endl;
  
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

float BMyRecoMC::getTrackDistanceToOM(TVector3 initialPoint, TVector3 direction, TVector3 xyzOM)
{
  TVector3 diff(xyzOM.X()-initialPoint.X(),xyzOM.Y()-initialPoint.Y(),xyzOM.Z()-initialPoint.Z());

  float absDir=direction.Mag();

  float scalarProd=diff.Dot(direction);

  float cosAlpha=scalarProd/(diff.Mag()*absDir);

  cosAlpha=round(cosAlpha*10000000)*pow(10000000,-1); 

  float dist=diff.Mag()*sqrt(1-cosAlpha*cosAlpha);

  return dist;
}

//return time of propagation from arbitrary point A on trajectory with direction s to coordinates M
float BMyRecoMC::getTimeEstimate_ns(TVector3 A, TVector3 s, TVector3 M)
{
  TVector3 AM(M.X()-A.X(),M.Y()-A.Y(),M.Z()-A.Z());
  float modAM=AM.Mag();
  float modS=s.Mag();
  //  std::cout<<"AM: "<<modAM<<"  s: "<<modS<<std::endl;

  float cosAlpha=AM.Dot(s)/(modAM*modS);
  cosAlpha=round(cosAlpha*1000000)*pow(1000000,-1); //rounding to avoid nan in the next line
  float sinAlpha=sqrt(1-pow(cosAlpha,2)); 

  //water refraction index
  float refWat=1.33;
  float cos_c=1/1.33;
  float sin_c=sqrt(1-cos_c*cos_c);
  float tan_c=sin_c/cos_c;

  
  //  std::cout<<"cos, sin: "<<cosAlpha<<"   "<<sinAlpha<<"  "<<cos_c<<"  "<<sin_c<<"  "<<tan_c<<std::endl;
  /*
  float dist=modAM*(cosAlpha - 2.47*sinAlpha);
  
  //  std::cout<<"dist:  "<<dist<<std::endl;
  float time=1e9*dist/(clight);
  */

  //separate distance in ligh and muon parts

  float dMuon=modAM*cosAlpha - modAM*sinAlpha/tan_c;
  float tMuon=1e9*dMuon/cVacuum;

  float dLight=modAM*sinAlpha/sin_c;
  float tLight=1e9*dLight/cWater;
  
  float time=tMuon+tLight;

  //  std::cout<<"ddddddddddddd"<<std::endl;

  return time;
}


int BMyRecoMC::RunMCAnalysis()
{
  
  float primThetaRad=M_PI*fMCEvent->GetPrimaryParticlePolar()/180;
  float primPhiRad=M_PI*fMCEvent->GetPrimaryParticleAzimuth()/180;

  hNmuResp->Fill(fMCEvent->GetResponseMuonsN(),1);
  hNmuTotal->Fill(fMCEvent->GetTotalMuonsN(),1);
  
  hPolar1muMC->Fill(fMCEvent->GetPrimaryParticlePolar(),1);
  
  float angle=-1;

  TVector3 genVec(sin(primThetaRad)*cos(primPhiRad),
		  sin(primThetaRad)*sin(primPhiRad),
		  cos(primThetaRad));

  TVector3 inPoTrack(fMCEvent->GetTrack(0)->GetX()-1000*genVec.X(), fMCEvent->GetTrack(0)->GetY()-1000*genVec.Y(), fMCEvent->GetTrack(0)->GetZ()-1000*genVec.Z());
  
  
  if (fMCEvent->GetResponseMuonsN()!=1) return 0;

  bool useEventForCalib=false;

  //count hits in the event:
  
  for (int iAmpl=0; iAmpl<50; iAmpl++){
    float AMPL=1+0.2*iAmpl;
    int nHits=0;
    int nDirectHits=0;
    for (int iCh=0; iCh<fMCEvent->GetChannelN(); iCh++){
      bool countChan=false;
      bool countChanDirect=false;
      for (int iPu=0; iPu<fMCEvent->GetHitChannel(iCh)->GetPulseN(); iPu++){
	if (fMCEvent->GetHitChannel(iCh)->GetPulse(iPu)->GetMagic()!=1){
	  if(fMCEvent->GetHitChannel(iCh)->GetPulse(iPu)->GetAmplitude()>AMPL) countChan=true;
	}
	
	if (fMCEvent->GetHitChannel(iCh)->GetPulse(iPu)->GetMagic()==-999&&iAmpl==5){
	  countChanDirect=true;
	}
      }
      if (countChan) nHits++;
      if (countChanDirect) nDirectHits++;
    }
    if (nHits>=5) hYieldVsAmpl_ref->Fill(AMPL+0.1,1);
    if (nDirectHits>=4) useEventForCalib=true;
  }
  
  TVector3 zero(0,0,0);
  float rho=getTrackDistanceToOM(inPoTrack, genVec, zero);
  hRho1muMC->Fill(rho,1);
  
  
  for (int i=0; i<fMCEvent->GetResponseMuonsN(); i++){
    int nStrongInt=0;
    for (int j=0; j<fMCEvent->GetTrack(0)->GetInteractionN(); j++){
      if (fMCEvent->GetTrack(0)->GetInteraction(j)->GetEnergy()>0.1) nStrongInt++;
    }
  }
  
  int nChans=fMCEvent->GetChannelN();
  int nPulses=0;
  int nDirectPulses=0;
  int nShowerPulses=0;
  int nNoisePulses=0;

  for (int i=0; i<nChans; i++){

    int idch=fMCEvent->GetHitChannel(i)->GetChannelID()-1;
    idch=24*floor(idch/24)+(24-idch%24);
    idch=idch-1;
    
    TVector3 chanXYZ(fGeomTel->At(idch)->GetX(),fGeomTel->At(idch)->GetY(),fGeomTel->At(idch)->GetZ());

    float time0=0;
    
    bool useChanDirect=false;
    bool useChanShower=false;
    int nPulsesPerChannel=0;
    for (int iP=0; iP<fMCEvent->GetHitChannel(i)->GetPulseN(); iP++){
      if (fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetMagic()+1000==1&&fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetAmplitude()>2) {
	useChanDirect=true;
	nDirectPulses++;
	time0+=fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetTime();
	nPulsesPerChannel++;
	nPulses++;
      }
      
      if (fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetMagic()>1&&fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetAmplitude()>2) {
	useChanShower=true;
	time0+=fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetTime();
	nPulsesPerChannel++;
	nShowerPulses++;
	if (iEvent==9507) {
	  hDebugTimeShowerHits->Fill(fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetTime(), (fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetMagic()-1)/1000, 1);
	  std::cout<<fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetMagic()-1000<<std::endl;
	}
	nPulses++;
      }
      
      if (fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetMagic()==1&&fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetAmplitude()>2) {
	nNoisePulses++;
	nPulses++;
      }
      hPulsesPerChannel->Fill(nPulsesPerChannel,1);
    }

    time0=time0/nPulsesPerChannel;
    
    hNhit1muMC->Fill(nPulses,1);
    hNdirect1muHit->Fill(nDirectPulses,1);
    hNshower1muHit->Fill(nShowerPulses,1);
    hNnoiseHit->Fill(nNoisePulses,1);

    if (!useChanDirect&&!useChanShower) continue;

    //if (!useChanDirect) continue;
    
    hNhitPerChannel->Fill(idch,1);

    //TVector3 chanXYZ(fGeomTel->At(idch)->GetX(),fGeomTel->At(idch)->GetY(),fGeomTel->At(idch)->GetZ());
    //hRealHitsALL->Fill(floor((idch)/24),(idch)%24+1);

    float dist=getTrackDistanceToOM(inPoTrack, genVec, chanXYZ);
    if (useChanDirect) hDistToTrackMC_direct->Fill(dist,1);
    if (useChanShower) hDistToTrackMC_shower->Fill(dist,1);

    //calculate time difference between two adjacent channels 
    //compare to actual time difference
    if (!useEventForCalib) continue;
    fNCalibTracks++;
    //std::cout<<"evt: "<<iEvent<<"    number of shower hits in channel #"<<idch<<":  "<<nPulsesPerChannel<<std::endl;
    //   if (fMCEvent->GetHitChannel(i+1)->GetPulseN()!=1) continue;
    //   if (fMCEvent->GetHitChannel(i+1)->GetPulse(0)->GetAmplitude()<3) continue;
    
    //    if (fMCEvent->GetHitChannel(i)->GetPulse(0)->GetMagic()==-999){

      //find tzero of muon not using this channel
    float tZERO=0;
    int nPulseTZ=0;
    for (int iCh=0; iCh<fMCEvent->GetChannelN(); iCh++){
      if (iCh==i) continue;
      int idch=fMCEvent->GetHitChannel(iCh)->GetChannelID()-1;
      idch=24*floor(idch/24)+(24-idch%24);
      idch=idch-1;
      
      TVector3 chanXYZ(fGeomTel->At(idch)->GetX(),fGeomTel->At(idch)->GetY(),fGeomTel->At(idch)->GetZ());
      
      for (int iPu=0; iPu<fMCEvent->GetHitChannel(iCh)->GetPulseN(); iPu++){
	if (fMCEvent->GetHitChannel(iCh)->GetPulse(iPu)->GetMagic()!=1&&fMCEvent->GetHitChannel(iCh)->GetPulse(iPu)->GetAmplitude()>2) {
	  tZERO+=fMCEvent->GetHitChannel(iCh)->GetPulse(iPu)->GetTime()-getTimeEstimate_ns(inPoTrack,genVec,chanXYZ);
	  nPulseTZ++;
	  //	    break;
	  
	}
      }
    }
    tZERO=tZERO/nPulseTZ;
    /////
    
    //      time0=fMCEvent->GetHitChannel(i)->GetPulse(0)->GetTime();
    
    float time1=tZERO+getTimeEstimate_ns(inPoTrack,genVec, chanXYZ);
    
    hTimeDifference_norm->Fill(time0-time1,1);
    hChanOffset[idch]->Fill(time0-time1,1);
    hDiffReal_DiffEst->Fill(time0,time1,1);

    if (time0-time1<-100) std::cout<<time0<<"   "<<time1<<std::endl;
    
    /*
      for (int iNext=i+1; iNext<nChans; iNext++){
	
	if (fMCEvent->GetHitChannel(iNext)->GetPulseN()==1&&
	    fMCEvent->GetHitChannel(iNext)->GetPulse(0)->GetMagic()!=1&&
	    fMCEvent->GetHitChannel(iNext)->GetPulse(0)->GetAmplitude()>=3){
	
	  time1=fMCEvent->GetHitChannel(i+1)->GetPulse(0)->GetTime();
	  
	  float dTim=time1-time0;
	  
	  int idchNext=fMCEvent->GetHitChannel(i+1)->GetChannelID()-1;
	  idchNext=24*floor(idchNext/24)+(24-idchNext%24);
	  idchNext=idchNext-1;
	  TVector3 chanXYZnext(fGeomTel->At(idchNext)->GetX(),fGeomTel->At(idchNext)->GetY(),fGeomTel->At(idchNext)->GetZ());
	  
	  //calculate time of propagation from zero point

	  //first define time in zero point averaging estimatied time of propagation to all other hits except for considered hit
	  	  
	  float timeZero=getTimeEstimate_ns(inPoTrack, genVec, chanXYZ);
	  float timeOne=getTimeEstimate_ns(inPoTrack, genVec, chanXYZnext);
	  float dTimEstimate=timeOne-timeZero;
	  if (time0!=0&&time1!=0&&timeZero!=0&&timeOne!=0) {
	    hTimeDifference_norm->Fill(dTim-dTimEstimate,1);
	    hChanOffset[idch]->Fill(dTim-dTimEstimate,1);
	    hDiffReal_DiffEst->Fill(dTim,dTimEstimate,1);
	  }
	  break;
	}
	
      }
      */
      
  }
  //     std::cout<<dTim<<"   "<<dTimEstimate<<"        debug: "<<time0<<" "<<time1<<"    est: "<<timeZero<<"  "<<timeOne<<std::endl;
  
  
  return 0;
}
