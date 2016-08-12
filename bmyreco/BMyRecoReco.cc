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

  hTimeDiff=new TH1F("hTimeDiff","hTimeDiff",10000,-5000,5000);

  hDistToTrack=new TH1F("hDistToTrack","hDistToTrack",1000,0,1000);  
}


BMyRecoReco::~BMyRecoReco()
{
  fOUT->cd();
  hTimeDiff->Write();
  hDistToTrack->Write();
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
    std::cout<<"iChan: "<<i<<"   X: "<<fGeomTel->At(i)->GetX()<<"   Y: "<<fGeomTel->At(i)->GetY()<<"   Z: "<<fGeomTel->At(i)->GetZ()<<std::endl;
  }

  fEvent=(BEvent*)pList->FindObject("BEvent","BEvent");
  if (!fEvent){
    * fLog << err << AddSerialNumber("BMCEvent") << " not found... aborting." << endl;
    return kFALSE;
  }
  
  fRecParam = (BRecParameters*)pList->FindObject("BRecParameters");
  //if(!fRecParam){ * fLog << err << AddSerialNumber("BRecParameters") << " not found... aborting." << endl;
  //    return kFALSE;
  //  }
  return kTRUE;
}


Int_t BMyRecoReco::Process()
{
  iEvent++;
  if (iEvent%10000==0) std::cout<<"eventNumber: "<<iEvent<<std::endl;

  //  std::cout<<fRecParam->GetNumPar()<<std::endl;
  //return kTRUE;

  //implement track selection wrt particular OM passed in the constructor
  
  //get track parameters:
  //1). Some initial point - problem, ask Fedor
  //2). direction
  //3). Time - problem, ask Fedor  
  
  //suppose we have distance: rTrackOM

  float thetaRad=M_PI*fRecParam->GetThetaRec()/180;
  float phiRad=M_PI*fRecParam->GetPhiRec()/180;
    
  //transform from X0, Y0 to X,Y,Z, copy-paste from BMuonX0Y0::DefineXYZ()
  TRotation m;
  m.RotateY(thetaRad - M_PI);
  m.RotateZ(phiRad);
  
  TVector3 R(fRecParam->GetX0Rec(), fRecParam->GetY0Rec(), 0);
 	
  R = m * R;
 	
  float initialX = R.X();
  float initialY = R.Y();
  float initialZ = R.Z();
  /////////////////////////////////

  std::vector<float> initialPoint;
  initialPoint.push_back(initialX);  //calculate according to 1).
  initialPoint.push_back(initialY);  //calculate according to 1).
  initialPoint.push_back(initialZ);  //calculate according to 1).

  std::vector<float> trackDirection;
  std::cout<<std::setprecision(2);
  //<<fRecParam->GetThetaRec()<<"   "<<thetaRad<<"   "<<cos(thetaRad)<<"   "<<fRecParam->GetPhiRec()<<"   "<<phiRad<<std::endl;
  trackDirection.push_back(sin(thetaRad)*cos(phiRad));
  trackDirection.push_back(sin(thetaRad)*sin(phiRad));
  trackDirection.push_back(cos(thetaRad));

  //  std::cout<<"track direction:  "<<trackDirection[0]<<"   "<<trackDirection[1]<<"   "<<trackDirection[2]<<std::endl;
    
  std::vector<float> tZeroPoint;
  tZeroPoint.push_back(initialX-trackDirection[0]*1000);
  tZeroPoint.push_back(initialY-trackDirection[1]*1000);
  tZeroPoint.push_back(initialZ-trackDirection[2]*1000);

  //  std::cout<<"tZeroPoint:  "<<tZeroPoint[0]<<"   "<<tZeroPoint[1]<<"   "<<tZeroPoint[2]<<std::endl;

  std::vector<float> xyzOM;
  xyzOM.push_back(fGeomTel->At(calibChanID)->GetX());
  xyzOM.push_back(fGeomTel->At(calibChanID)->GetY());
  xyzOM.push_back(fGeomTel->At(calibChanID)->GetZ());

  float rTrackOM=getTrackDistanceToOM(tZeroPoint, trackDirection, xyzOM); //meters

  hDistToTrack->Fill(rTrackOM,1);
  
  if (rTrackOM>10)
    {
      //  std::cout<<"MISS channel 10"<<std::endl;
      return kTRUE;
    }
  // std::cout<<"hit channel 10: "<<rTrackOM<<std::endl;
  
  //find absolute (inside the event) time of hit in the OM of interest
  float time0=0; //time0=GetTimeRec() BUT THERE IS NO SUCH FUNCTION!!! ASK FEDOR

  //get time from residuals
  
  float T0;
  for (int i=0; i<fRecParam->GetNhit(); i++){
    float Tres=fRecParam->GetTres(i);
    int nch=fRecParam->GetNchGeomMC(i);
    std::vector<float> xyzHit;
    xyzHit.push_back(fGeomTel->At(nch)->GetX());
    xyzHit.push_back(fGeomTel->At(nch)->GetY());
    xyzHit.push_back(fGeomTel->At(nch)->GetZ());

    //propagationTime to hit from tZeroPoint:
    float Tprop=getTimeEstimate_ns(tZeroPoint,trackDirection,xyzHit);

    //experimental time:
    int impID=fRecParam->GetImpulseNumber(i);
    float Texp=fEvent->GetImpulseTime(impID);

    T0=Texp-Tres-Tprop;
    //std::cout<<"T0 from channel "<<nch<<"  and impulse "<<impID<<" :   "<<T0<<std::endl;
  }
  
  float timeOfOMHit_estimate=T0+getTimeEstimate_ns(tZeroPoint,trackDirection, xyzOM);

  //here one should look at BMCEvent or BEvent
  //loop over hits
  float timeOfOMHit_actual=-1;
  int nImpulseInCalibChan=0;
  
  for (int i=0; i<fEvent->GetTotImpulses(); i++) {
    if (calibChanID==fEvent->GetImpulse(i)->GetChannelID()){
      timeOfOMHit_actual=fEvent->GetImpulse(i)->GetTime();
      nImpulseInCalibChan++;
    }
  }
  
  hTimeDiff->Fill(timeOfOMHit_estimate-timeOfOMHit_actual,1);  
  
  return kTRUE;
  
}

//Float_t getDistance(std::vector<float> 


Int_t BMyRecoReco::PostProcess()
{
  return kTRUE;
}

//distance from track to OM
float BMyRecoReco::getTrackDistanceToOM(std::vector<float> A, std::vector<float> a, std::vector<float> xyzOM)
{
  float modAOM=sqrt(pow(A[0]-xyzOM[0],2)+pow(A[1]-xyzOM[1],2)+pow(A[2]-xyzOM[2],2));
  float moda=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  float scalarMult=((xyzOM[0]-A[0])*a[0]+(xyzOM[1]-A[1])*a[1]+(xyzOM[2]-A[2])*a[2]);
  float cosAlpha=scalarMult/(modAOM*moda);
  cosAlpha=round(cosAlpha*10000000)*pow(10000000,-1); //rounding to avoid nan in the next line
  float dist=modAOM*sqrt(1-cosAlpha*cosAlpha);

  //if (dist<10) std::cout<<modAOM<<"   "<<moda<<"   "<<scalarMult<<"   "<<cosAlpha<<"   "<<sqrt(1-cosAlpha*cosAlpha)<<"   "<<dist<<std::endl;
  
  return dist;
}


//return time of propagation from arbitrary point A on trajectory with direction s to coordinates M
float BMyRecoReco::getTimeEstimate_ns(std::vector<float> A, std::vector<float> s, std::vector<float> M)
{
  std::vector<float> AM;
  AM.push_back(M[0]-A[0]);
  AM.push_back(M[1]-A[1]);
  AM.push_back(M[2]-A[2]);
  float modAM=sqrt(pow(AM[0],2)+pow(AM[1],2)+pow(AM[2],2));
  float modS=sqrt(pow(s[0],2)+pow(s[1],2)+pow(s[2],2));
  //  std::cout<<"AM: "<<modAM<<"  s: "<<modS<<std::endl;

  float cosAlpha=(AM[0]*s[0]+AM[1]*s[1]+AM[2]*s[2])/(modAM*modS);
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
