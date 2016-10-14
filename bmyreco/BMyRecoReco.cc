#include "BMyRecoReco.h"
#include "BGeomTel.h"
#include "BEvent.h"
#include "BMCEvent.h"
#include "BExtractedImpulse.h"
#include "BExtractedImpulseTel.h"
#include "BJoinExtractedImpulseTel.h"
#include "BEventMask.h"
#include "BChannelMask.h"
#include "BExtractedHeader.h"

#include "BRecParameters.h"

#include "MParList.h"
#include "MLog.h"
#include "MLogManip.h"
 
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TMath.h"
#include "TRotation.h"
#include "TVector3.h"

#include <math.h>
#include <vector>

ClassImp(BMyRecoReco);

BMyRecoReco::BMyRecoReco(string fname, int cChan, bool useMCEvent)
{
  iEvent=0;
 
  fOUT=new TFile(fname.c_str(),"RECREATE");

  cWater=3e8*pow(1.33,-1);
  cVacuum=3e8;
  
  calibChanID=cChan;

  //plots in order
  hMuMult=new TH1F("hMuMult","hMuMult",50,0,50);
  hPolar1muMC=new TH1F("hPolar1muMC","hPolar1muMC",1000,-500,500);
  hPolar1muRec=new TH1F("hPolar1muRec","hPolar1muRec",1000,-500,500);
  hRho1muMC=new TH1F("hRho1muMC","hRho1muMC",1000,-500,500);
  hNhitMC=new TH1F("hNhitMC","hNhitMC",100,0,100);
  hNhit=new TH1F("hNhit","hNhit",100,0,100);
  hNdirectHit=new TH1F("hNdirectHit","hNdirectHit",100,0,100);
  hNshowerHit=new TH1F("hNshowerHit","hNshowerHit",100,0,100);
  hNnoiseHit=new TH1F("hNnoiseHit","hNnoiseHit",100,0,100);
  hNhitPerChannel=new TH1F("hNhitPerChannel","hNhitPerChannel",193,0,193);
  
  
  hNhit1mu=new TH1F("hNhit1mu","hNhit1mu",100,0,100);
  hChi2=new TH1F("hChi2","hChi2",10000,0,1000);
  hChi2_zoom1=new TH1F("hChi2_zoom1","hChi2_zoom1",1000,0,1);
  hTimeDiff=new TH1F("hTimeDiff","hTimeDiff",10000,-5000,5000);
  hTimeDiff_usedHits=new TH1F("hTimeDiff_usedHits","hTimeDiff_usedHits",10000,-5000,5000);
  hDistToTrack=new TH1F("hDistToTrack","hDistToTrack",1000,0,1000);
  hDistToTrack_usedHits=new TH1F("hDistToTrack_usedHits","hDistToTrack_usedHits",1000,0,1000);
  hDistToTrack1mu_usedHits=new TH1F("hDistToTrack1mu_usedHits","hDistToTrack1mu_usedHits",1000,0,1000);

  hDistToTrackMC=new TH1F("hDistToTrackMC","hDistToTrackMC",1000,0,1000);

  hDistToTrack_debug=new TH1F("hDistToTrack_debug","hDistToTrack_debug",10000,0,100);
  
  hChi2_rhoAvg=new TH2F("hChi2_rhoAvg","hChi2_rhoAvg",1000,0,100,1000,0,300);
  hChi2_rhoMax=new TH2F("hChi2_rhoMax","hChi2_rhoMax",1000,0,100,1000,0,300);

  hAngle1muGenRec=new TH1F("hAngle1muGenRec","hAngle1muGenRec",1000,-500,500);

  hMagNum_usedHits=new TH1F("hMagNum_usedHits","hMagNum_usedHits",3,0,3);

  hHitsDR_MC=new TH1F("hHitsDR_MC","hHitsDR_MC",1000,0,1000);

  hExtrapolatedHits=new TH2F("hExtrapolatedHits","hExtrapolatedHits",21,-0.5,20.5,25,-0.5,24.5);
  hRealHits=new TH2F("hRealHits","hRealHits",21,-0.5,20.5,25,-0.5,24.5);

  
  hExtrapolatedHitsALL=new TH2F("hExtrapolatedHitsALL","hExtrapolatedHitsALL",21,-0.5,20.5,25,-0.5,24.5);
  hRealHitsALL=new TH2F("hRealHitsALL","hRealHitsALL",21,-0.5,20.5,25,-0.5,24.5);

  hNumMC_NumReco=new TH2F("hNumMC_NumReco","hNumMC_NumReco",50,0,50,50,0,50);
  hDMC_dReco=new TH2F("hDMC_dReco","hDMC_dReco",312,-11,300, 312,-11,312);
}




BMyRecoReco::~BMyRecoReco()
{
  fOUT->cd();
  hMuMult->Write();
  hPolar1muMC->Write();
  hPolar1muRec->Write();
  hRho1muMC->Write();
  hNhitMC->Write();
  hNhit->Write();
  hNdirectHit->Write();
  hNshowerHit->Write();
  hNnoiseHit->Write();
  hNhitPerChannel->Write();
  
  
  hNhit1mu->Write();
  hTimeDiff->Write();
  hTimeDiff_usedHits->Write();
  hDistToTrack->Write();
  hDistToTrack_usedHits->Write();
  hDistToTrack1mu_usedHits->Write();
  hChi2->Write();
  hChi2_zoom1->Write();
  hChi2_rhoAvg->Write();
  hChi2_rhoMax->Write();
  hAngle1muGenRec->Write();
  hMagNum_usedHits->Write();

  hDistToTrackMC->Write();

  hDistToTrack_debug->Write();

  hHitsDR_MC->Write();

  hExtrapolatedHits->Write();
  hRealHits->Write();

  hExtrapolatedHitsALL->Write();
  hRealHitsALL->Write();

  hNumMC_NumReco->Write();
  hDMC_dReco->Write();
  
  fOUT->Close();
}

Int_t BMyRecoReco::PreProcess(MParList * pList)
{
  std::cout<<"WE ARE IN BMyReco Reco::PreProcess"<<std::endl;
  
  fGeomTel = (BGeomTel*)pList->FindObject(AddSerialNumber("BGeomTel"));
  if (!fGeomTel){
    * fLog << err << AddSerialNumber("BGeomTel") << " not found... aborting." << endl;
    return kFALSE;
  }
  
  for (int i=0; i<192; i++){
    std::cout<<"iChan: "<<i<<"   floor(nch/24): "<<floor(i/24)<<"  nch%24+1: "<<i%24+1<<"   X: "<<fGeomTel->At(i)->GetX()<<"   Y: "<<fGeomTel->At(i)->GetY()<<"   Z: "<<fGeomTel->At(i)->GetZ()<<std::endl;
  }
  

  fEvent=(BEvent*)pList->FindObject("BEvent","BEvent");
  if (!fEvent){
    * fLog << err << AddSerialNumber("BEvent") << " not found... aborting." << endl;
    return kFALSE;
  }
  
  fRecParam = (BRecParameters*)pList->FindObject("BRecParameters");
  //if(!fRecParam){ * fLog << err << AddSerialNumber("BRecParameters") << " not found... aborting." << endl;
  //    return kFALSE;
  //  }

  fMCEvent=(BMCEvent*)pList->FindObject("BMCEvent","BMCEvent");
  if (!fEvent){
    * fLog << err << AddSerialNumber("BMCEvent") << " not found... aborting." << endl;
    return kFALSE;
  }

  fEventMask= (BEventMask*)pList->FindObject("MuonCriterionFilterMask");
  
  return kTRUE;
}


Int_t BMyRecoReco::Process()
{
  //std::cout<<"start"<<std::endl;
  //debug distance to track
  TVector3 vtxDeb(10,0,1000);
  TVector3 speed(0,0,-1);
  
  for (int id=168; id<192; id++){
    TVector3 cooOM(fGeomTel->At(id)->GetX(),fGeomTel->At(id)->GetY(),fGeomTel->At(id)->GetZ());
    hDistToTrack_debug->Fill(getTrackDistanceToOM(vtxDeb, speed, cooOM),1);
  }
  /////////////////
  iEvent++;
  if (iEvent%10000==0) std::cout<<"eventNumber: "<<iEvent<<std::endl;

  hNhit->Fill(fRecParam->GetNhit(),1);
  hChi2->Fill(fRecParam->GetFuncValue(),1);
  hChi2_zoom1->Fill(fRecParam->GetFuncValue(),1);

  float thetaRad=fRecParam->GetThetaRec();
  float phiRad=fRecParam->GetPhiRec();

  float primThetaRad=M_PI*fMCEvent->GetPrimaryParticlePolar()/180;
  float primPhiRad=M_PI*fMCEvent->GetPrimaryParticleAzimuth()/180;
  
  float angle=-1;

  TVector3 genVec(sin(primThetaRad)*cos(primPhiRad),
		  sin(primThetaRad)*sin(primPhiRad),
		  cos(primThetaRad));
  
  
  TVector3 recVec(sin(thetaRad)*cos(phiRad),
		  sin(thetaRad)*sin(phiRad),
		  cos(thetaRad));


  RunMCAnalysis();

  
  //  if (fMCEvent->GetResponseMuonsN()==1&&fMCEvent->GetMuonsN()==1&&fMCEvent->GetTotalMuonsN()==1){
  if (fMCEvent->GetResponseMuonsN()==1){
    
    //filter muon based on energy and total interaction energy
    
    hPolar1muRec->Fill(180*fRecParam->GetThetaRec()/M_PI,1);
    
    angle = 180*genVec.Angle(recVec)/M_PI;
    hAngle1muGenRec->Fill(angle,1);
    if (angle<10&&angle>0) {
      hNhit1mu->Fill(fRecParam->GetNhit(),1);
      // std::cout<<"fill"<<std::endl;
    }  
  }    
  //return kTRUE;

  //implement track selection wrt particular OM passed in the constructor
  
  //get track parameters:
  //1). Some initial point - problem, ask Fedor
  //2). direction
  //3). Time - problem, ask Fedor  
  
  //suppose we have distance: rTrackOM

  //transform from X0, Y0 to X,Y,Z, copy-paste from BMuonX0Y0::DefineXYZ()
  /*
  TRotation m;
  m.RotateY(thetaRad - M_PI);
  m.RotateZ(phiRad);
  
  TVector3 R(fRecParam->GetX0Rec(), fRecParam->GetY0Rec(), 0);
  */	
  TVector3 R = xyToXYZ(recVec,fRecParam->GetX0Rec(), fRecParam->GetY0Rec());  
  //m * R;
  
  /////////////////////////////////

  //  std::cout<<(pow(fRecParam->GetX0Rec(),2)+pow(fRecParam->GetY0Rec(),2))<<"      "<<pow(initialX,2)+pow(initialY,2)+pow(initialZ,2)<<std::endl;

  //std::cout<<"scalar prod after back rotation: "<<
  
  TVector3 initialPoint(R.X(),R.Y(),R.Z());
 
  TVector3 trackDirection(sin(thetaRad)*cos(phiRad), sin(thetaRad)*sin(phiRad), cos(thetaRad));
  std::cout<<std::setprecision(2);
  //<<fRecParam->GetThetaRec()<<"   "<<thetaRad<<"   "<<cos(thetaRad)<<"   "<<fRecParam->GetPhiRec()<<"   "<<phiRad<<std::endl;

  //std::cout<<"scalar prod after back rotation: "<<initialX*trackDirection[0]+initialY*trackDirection[1]+initialZ*trackDirection[2]<<std::endl;
  
  //std::cout<<"track direction:  "<<trackDirection.X()<<"   "<<trackDirection.Y()<<"   "<<trackDirection.Z()<<std::endl;
    
  TVector3 tZeroPoint(initialPoint.X()-trackDirection.X()*1000, initialPoint.Y()-trackDirection.Y()*1000, initialPoint.Z()-trackDirection.Z()*1000);

  //std::cout<<"tZeroPoint:  "<<tZeroPoint.X()<<"   "<<tZeroPoint.Y()<<"   "<<tZeroPoint.Z()<<std::endl;
  //DEBUG:
  TVector3 center(0,0,0);

//  std::cout<<"distance to center:  "
// 	   <<sqrt(pow(initialX,2)+pow(initialY,2)+pow(initialZ,2))
//  	   <<"     analytic from tZeroPoint: "
    // <<getTrackDistanceToOM(tZeroPoint, trackDirection, center)
  //  	   <<std::endl;
  
  TVector3 xyzOM(fGeomTel->At(calibChanID)->GetX(), fGeomTel->At(calibChanID)->GetY(), fGeomTel->At(calibChanID)->GetZ());


  //  std::cout<<"reco: "<<getTrackDistanceToOM(tZeroPoint, trackDirection, center)<<"   "<<sqrt(pow(fRecParam->GetX0Rec(),2)+pow(fRecParam->GetY0Rec(),2))<<std::endl;
  //std::cout<<"xyzOM: "<<xyzOM.X()<<"  "<<xyzOM.Y()<<"  "<<xyzOM.Z()<<std::endl;
  //float rTrackOM=getTrackDistanceToOM(tZeroPoint, trackDirection, xyzOM); //meters

  //hDistToTrack->Fill(rTrackOM,1);
  
  //  if (rTrackOM>10)
  // {
      //  std::cout<<"MISS channel 10"<<std::endl;
      //return kTRUE;
  //  }
  // std::cout<<"hit channel 10: "<<rTrackOM<<std::endl;
  
  //find absolute (inside the event) time of hit in the OM of interest
  float time0=0; //time0=GetTimeRec() BUT THERE IS NO SUCH FUNCTION!!! ASK FEDOR

  //get time from residuals
  
  float T0;
  float rhoAvg=0;
  float rhoMax=-1; 
  for (int i=0; i<fRecParam->GetNhit(); i++){
    float Tres=fRecParam->GetTres(i);
    int nch=fRecParam->GetNchGeom(i);
    TVector3 xyzHit(fGeomTel->At(nch)->GetX(), fGeomTel->At(nch)->GetY(), fGeomTel->At(nch)->GetZ());
    
    //propagationTime to hit from tZeroPoint:
    float Tprop=getTimeEstimate_ns(tZeroPoint,trackDirection,xyzHit);
    
    //DEBUG: distance from track to hit:
    float distUsed=getTrackDistanceToOM(tZeroPoint, trackDirection, xyzHit);
    if (fMCEvent->GetResponseMuonsN()==1){
      //std::cout<<"hits: "<<distUsed<<" "<<angle<<std::endl;
      //      if (angle<10&&angle>0){
      hDistToTrack1mu_usedHits->Fill(distUsed,1);
				       //fRecParam->GetRho(i),1);
	//	std::cout<<"hits fill"<<std::endl;
	//      }
    }
    hDistToTrack_usedHits->Fill(distUsed,1);
				//fRecParam->GetRho(i),1);
    
    //study chi2 vs Rho
    rhoAvg+=fRecParam->GetRho(i);
    if (fRecParam->GetRho(i)>rhoMax) rhoMax=fRecParam->GetRho(i);
    
    //experimental time:
    int impID=fRecParam->GetImpulseNumber(i);
    float Texp=fEvent->GetImpulseTime(impID);

    int magNum=-1;
    
    for (int iMCch=0; iMCch<fMCEvent->GetChannelN(); iMCch++){
      for (int iPulse=0; iPulse<fMCEvent->GetHitChannel(iMCch)->GetPulseN(); iPulse++){
	if (Texp==fMCEvent->GetHitChannel(iMCch)->GetPulse(iPulse)->GetTime()) magNum=fMCEvent->GetHitChannel(iMCch)->GetPulse(iPulse)->GetMagic();
      }
    }

    if (magNum==1) hMagNum_usedHits->Fill(0.5,1);
    if (magNum!=1&&magNum!=-1) hMagNum_usedHits->Fill(1.5,1);
    if (magNum==-1) hMagNum_usedHits->Fill(2.5,1);
    
    T0=Texp-Tres-Tprop;
    //std::cout<<"T0 from channel "<<nch<<"  and impulse "<<impID<<" :   "<<T0<<std::endl;
  }

  rhoAvg=rhoAvg/fRecParam->GetNhit();
  hChi2_rhoAvg->Fill(fRecParam->GetFuncValue(),rhoAvg,1);
  hChi2_rhoMax->Fill(fRecParam->GetFuncValue(),rhoMax,1);
  
  float timeOfOMHit_estimate=T0+getTimeEstimate_ns(tZeroPoint,trackDirection, xyzOM);

  //here one should look at BMCEvent or BEvent
  //loop over hits
  float timeOfOMHit_actual=-1;
  int nImpulseInCalibChan=0;
  
  for (int i=0; i<fEvent->GetTotImpulses(); i++) {
    if (calibChanID==fEvent->GetImpulse(i)->GetChannelID()){
      timeOfOMHit_actual=fEvent->GetImpulseTime(i);
      nImpulseInCalibChan++;
    }
  }
  
  if (nImpulseInCalibChan==1) hTimeDiff->Fill(timeOfOMHit_estimate-timeOfOMHit_actual,1);  

  //DEBUG
  for (int i=0; i<fRecParam->GetNhit(); i++){
    int timeOfUsedHit=fEvent->GetImpulseTime(fRecParam->GetImpulseNumber(i));
    int usedChanID=fRecParam->GetNchGeom(i);
    TVector3 xyzUsedHit(fGeomTel->At(usedChanID)->GetX(),fGeomTel->At(usedChanID)->GetY(),fGeomTel->At(usedChanID)->GetZ());
    float timeOfUsedHit_theor=T0+getTimeEstimate_ns(tZeroPoint, trackDirection, xyzUsedHit);
    hTimeDiff_usedHits->Fill(timeOfUsedHit_theor-timeOfUsedHit,1);
  }
  
  return kTRUE;
  
}

//Float_t getDistance(std::vector<float> 


Int_t BMyRecoReco::PostProcess()
{
  return kTRUE;
}

float BMyRecoReco::getTrackDistanceToOM(TVector3 initialPoint, TVector3 direction, TVector3 xyzOM)
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
float BMyRecoReco::getTimeEstimate_ns(TVector3 A, TVector3 s, TVector3 M)
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

TVector3 BMyRecoReco::xyToXYZ(TVector3 s, float X, float Y)
{
  //std::cout<<"mag: "<<s.Mag()<<std::endl;
  float thetaRad=acos(s.Z());
  float phiRad=acos(s.X()/sqrt(s.X()*s.X()+s.Y()*s.Y()));
  
  //transform from X0, Y0 to X,Y,Z, copy-paste from BMuonX0Y0::DefineXYZ()
  TRotation m;
  m.RotateY(thetaRad - M_PI);
  m.RotateZ(phiRad);
  
  TVector3 R(X, Y, 0);
 	
  R = m * R;

  return R;
}
  //  float initialX = R.X();
  //  float initialY = R.Y();
  //  float initialZ = R.Z();
  /////////////////////////////////


int BMyRecoReco::RunMCAnalysis()
{
  float thetaRad=fRecParam->GetThetaRec();
  float phiRad=fRecParam->GetPhiRec();

  float primThetaRad=M_PI*fMCEvent->GetPrimaryParticlePolar()/180;
  float primPhiRad=M_PI*fMCEvent->GetPrimaryParticleAzimuth()/180;

  hMuMult->Fill(fMCEvent->GetResponseMuonsN(),1);
  hPolar1muMC->Fill(fMCEvent->GetPrimaryParticlePolar(),1);
  
  float angle=-1;

  TVector3 genVec(sin(primThetaRad)*cos(primPhiRad),
		  sin(primThetaRad)*sin(primPhiRad),
		  cos(primThetaRad));
  
  if (fMCEvent->GetResponseMuonsN()!=1) return 0;

  TVector3 inPoTrack(fMCEvent->GetTrack(0)->GetX()-1000*genVec.X(), fMCEvent->GetTrack(0)->GetY()-1000*genVec.Y(), fMCEvent->GetTrack(0)->GetZ()-1000*genVec.Z());

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
    bool useChan=false;
    for (int iP=0; iP<fMCEvent->GetHitChannel(i)->GetPulseN(); iP++){
      if (fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetMagic()+1000==1) {
	useChan=true;
	nDirectPulses++;
	nPulses++;
      }
      if (fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetMagic()>1) {
	//useChan=true;
	nShowerPulses++;
	nPulses++;
      }
      if (fMCEvent->GetHitChannel(i)->GetPulse(iP)->GetMagic()==1) {
	nNoisePulses++;
	nPulses++;
      }
    }
    
    hNhitMC->Fill(nPulses,1);
    hNdirectHit->Fill(nDirectPulses,1);
    hNshowerHit->Fill(nShowerPulses,1);
    hNnoiseHit->Fill(nNoisePulses,1);

    if (!useChan) continue;

    int idch=fMCEvent->GetHitChannel(i)->GetChannelID()-1;
    idch=24*floor(idch/24)+(24-idch%24);
    idch=idch-1;

    hNhitPerChannel->Fill(idch,1);

    TVector3 chanXYZ(fGeomTel->At(idch)->GetX(),fGeomTel->At(idch)->GetY(),fGeomTel->At(idch)->GetZ());
    //hRealHitsALL->Fill(floor((idch)/24),(idch)%24+1);

    float dist=getTrackDistanceToOM(inPoTrack, genVec, chanXYZ);
    hDistToTrackMC->Fill(dist,1);    
  }

  

  return 0;
}
