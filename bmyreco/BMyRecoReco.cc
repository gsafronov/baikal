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
  fName  = "BMyRecoReco";
  fTitle = "RMyRecoReco";
  
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
  hDistToTrackMC=new TH1F("hDistToTrackMC","hDistToTrackMC",1000,0,1000);
  hYieldVsAmpl_reco=new TH1F("hYieldVsAmpl_reco","hYieldVsAmpl_reco",50,1,11);
  
  hDebugRCenter=new TH1F("hDebugRCenter","hDebugRCenter",1000,0,10);
  hDebugRhoFedor=new TH1F("hDebugRhoFedor","hDebugRhoFedor",1000,0,1000);

  hAngleGenRec=new TH1F("hAngleGenRec","hAngleGenRec",600,-300,300);
  hAngle1muGenRec=new TH1F("hAngle1muGenRec","hAngle1muGenRec",600,-300,300);

  
  hNhit1mu=new TH1F("hNhit1mu","hNhit1mu",100,0,100);
  hChi2=new TH1F("hChi2","hChi2",10000,0,1000);
  hChi2_zoom1=new TH1F("hChi2_zoom1","hChi2_zoom1",1000,0,1);
  hTimeDiff=new TH1F("hTimeDiff","hTimeDiff",10000,-5000,5000);
  hTimeDiff_usedHits=new TH1F("hTimeDiff_usedHits","hTimeDiff_usedHits",10000,-5000,5000);
  hDistToTrack=new TH1F("hDistToTrack","hDistToTrack",1000,0,1000);
  hDistToTrack_usedHits=new TH1F("hDistToTrack_usedHits","hDistToTrack_usedHits",1000,0,1000);
  hDistToTrack1mu_usedHits=new TH1F("hDistToTrack1mu_usedHits","hDistToTrack1mu_usedHits",1000,0,1000);

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
  /*
  hRho1muMC->Write();
  hNhitMC->Write();
  hNhit->Write();
  hNdirectHit->Write();
  hNshowerHit->Write();
  hNnoiseHit->Write();
  hNhitPerChannel->Write();
  hYieldVsAmpl_reco->Write();
  */

  
  //hDebugRCenter->Write();
  //hDebugRhoFedor->Write();

  hAngleGenRec->Write();

  //hAngle1muGenRec->Write();

  /*
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
  */

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

  fEventMask= (BEventMask*)pList->FindObject("TimeClusterFilterMask");
  
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

  float primThetaRad=M_PI*(180-fMCEvent->GetPrimaryParticlePolar())/180;
  float primPhiRad=M_PI*(fMCEvent->GetPrimaryParticleAzimuth())/180;
  
  float angle=-1;

  TVector3 genVec(sin(primThetaRad)*cos(primPhiRad),
		  sin(primThetaRad)*sin(primPhiRad),
		  cos(primThetaRad));
  
  
  TVector3 recVec(sin(thetaRad)*cos(phiRad),
		  sin(thetaRad)*sin(phiRad),
		  cos(thetaRad));

  hPolar1muMC->Fill(180-fMCEvent->GetPrimaryParticlePolar());
  hPolar1muRec->Fill(180*fRecParam->GetThetaRec()/M_PI,1);
  hAngleGenRec->Fill(180*genVec.Angle(recVec)/M_PI,1);


  //Analyse fEventMask
  //1) find how many tracks are matched to the cluster
  
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
 
