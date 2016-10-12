#include "BMyRecoChanSetter.h"
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

#include <math.h>
#include <vector>

ClassImp(BMyRecoChanSetter);

BMyRecoChanSetter::BMyRecoChanSetter(const char *name, const char *title)
{
  //  chanIDs=chans;
  // fMaskNoise=maskNoise;
  fInputMaskName ="";
  fOutputMaskName="";

  fName  = name ? name : "BMyRecoChanSetter";
  fTitle = title ? title : "MCNoiseFilter";
  
  //  fOutputMaskName = "AmplitudeFilterMask";
}

BMyRecoChanSetter::~BMyRecoChanSetter()
{
}

Int_t BMyRecoChanSetter::PreProcess(MParList * pList)
{
  //  std::cout<<"WE ARE IN BMyRecoBEvt::PreProcess"<<std::endl;
  
  fMCEvent=(BMCEvent*)pList->FindObject("BMCEvent","BMCEvent");
  if (!fEvent){
    * fLog << err << AddSerialNumber("BMCEvent") << " not found... aborting." << endl;
    return kFALSE;
  }

  fEvent=(BEvent*)pList->FindObject("BEvent","BEvent");
  if (!fEvent){
    * fLog << err << AddSerialNumber("BEvent") << " not found... aborting." << endl;
    return kFALSE;
  }
  
  fChannelMask = (BChannelMask*)pList->FindObject("BChannelMask");
  if (!fChannelMask)
    {
      * fLog << err << AddSerialNumber("BChannelMask") << " not found... aborting." << endl;
      return kFALSE;
    }
  
  fInputEventMask = (BEventMask*)pList->FindObject(fInputMaskName);
  if (!fInputEventMask)
    {
      * fLog << err << AddSerialNumber(fInputMaskName) << " not found... aborting." << endl;
      return kFALSE;
    }
  
  /*
  std::cout<<"BChanSetter: exclude calibrated channel from reconstruction"<<std::endl;
  for (int i=0; i<chanIDs.size(); i++){
    if (fChannelMask->GetFlag(chanIDs[i])==1) fChannelMask->SetFlag(chanIDs[i], BChannelMask::EMode::kOff);
  }
  */
  if (!BFilter::PreProcess(pList)) {
    return kFALSE;
  }

  return kTRUE;
}


Bool_t BMyRecoChanSetter::Filter()
{
  //  std::cout<<"in filter"<<std::endl;
  int n_impulse=fMCEvent->GetChannelN();
  
  //  std::cout<<"MASK BEFORE"<<std::endl;
  for (int i=0; i<fEvent->GetTotImpulses(); i++){
    bool isNoise=false;
    //std::cout<<"impulse: "<<i<<"   channel: "<<fEvent->GetImpulse(i)->GetChannelID()<<"   mask: "<<fInputEventMask->GetFlag(i);
    
    for (int j=0; j<n_impulse; j++){
      BMCHitChannel* fHitChan=fMCEvent->GetHitChannel(j);
      for (int k=0; k<fHitChan->GetPulseN(); k++){
	if (fEvent->GetImpulse(i)->GetTime()==fHitChan->GetPulse(k)->GetTime()){
	  if (fHitChan->GetPulse(k)->GetMagic()==1) isNoise=true;
	  fOutputEventMask->SetFlag(i,0);
	  // std::cout<<"   magic number: "<<fHitChan->GetPulse(k)->GetMagic()<<std::endl;
	}
      }
    }
    if (isNoise) fOutputEventMask->SetFlag(i,0);
    else fOutputEventMask->SetFlag(i,fInputEventMask->GetFlag(i));
  }
    
    
  int noisePulseID=-100;
  for (int i=0; i<n_impulse; i++) {
    BMCHitChannel* fHitChan=fMCEvent->GetHitChannel(i);
    for (int j=0; j<fHitChan->GetPulseN(); j++){
      //std::cout<<fHitChan->GetPulse(j)->GetMagic()<<std::endl;
      if (fHitChan->GetPulse(j)->GetMagic()==1) {
	//	int noisePulseID=0;
	float time=fHitChan->GetPulse(j)->GetTime();
	for (int k=0; k<fEvent->GetTotImpulses(); k++){
	  if (fEvent->GetImpulse(k)->GetTime()==time) {
	    // std::cout<<fEvent->GetImpulse(k)->GetTime()<<"   "<<time<<std::endl;
	      noisePulseID=k;
	  }
	}
	//	fInputEventMask->SetFlag(noisePulseID,-2);
	//std::cout<<"impulse: "<<noisePulseID<<"   channel: "<<fEvent->GetImpulse(noisePulseID)->GetChannelID()<<"   mask: "<<fInputEventMask->GetFlag(i)<<"   magic: "<<fHitChan->GetPulse(j)->GetMagic()<<std::endl;
      }
      //std::cout<<"impulse: "<<noisePulseID<<"   channel: "<<fEvent->GetImpulse(noiePulseID)->GetChannelID()<<<<std::endl;
    }
  }

  /*
  std::cout<<"MASK AFTER"<<std::endl;
  for (int i=0; i<fEvent->GetTotImpulses(); i++){
    std::cout<<"impulse: "<<i<<"    channel: "<<fEvent->GetImpulse(i)->GetChannelID()<<"   mask: "<<fOutputEventMask->GetFlag(i)<<"   "<<noisePulseID<<std::endl;
  }
  */
  
  /*
  std::cout<<"BChanSetter: exclude calibrated channel from reconstruction"<<std::endl;
  for (int i=0; i<chanIDs.size(); i++){
    if (fChannelMask->GetFlag(chanIDs[i])==1) fChannelMask->SetFlag(chanIDs[i], BChannelMask::EMode::kOff);
  }
  */
  
  
  
  return kTRUE;
}

/*
Int_t BMyRecoChanSetter::PostProcess()
{
  return kTRUE;
}
*/
