#include "BTimeClusters.h"
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

ClassImp(BTimeClusters);

BTimeClusters::BTimeClusters(const char *name, const char *title)
{
  //  chanIDs=chans;
  // fMaskNoise=maskNoise;
  fInputMaskName ="";
  fOutputMaskName="";

  fName  = name ? name : "BMyRecoChanSetter";
  fTitle = title ? title : "MCNoiseFilter";
  
  //  fOutputMaskName = "AmplitudeFilterMask";
}

BTimeClusters::~BTimeClusters()
{
}

Int_t BTimeClusters::PreProcess(MParList * pList)
{
  fMCEvent=(BMCEvent*)pList->FindObject("BMCEvent","BMCEvent");
  
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


Bool_t BTimeClusters::Filter()
{
  int n_impulse=fMCEvent->GetChannelN();
  
  for (int i=0; i<fEvent->GetTotImpulses(); i++){
    bool isNoise=false;

    for (int j=0; j<n_impulse; j++){
      BMCHitChannel* fHitChan=fMCEvent->GetHitChannel(j);
      for (int k=0; k<fHitChan->GetPulseN(); k++){
	if (fEvent->GetImpulse(i)->GetTime()==fHitChan->GetPulse(k)->GetTime()){
	  if (fHitChan->GetPulse(k)->GetMagic()==1) isNoise=true;
	  fOutputEventMask->SetFlag(i,0);
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
      if (fHitChan->GetPulse(j)->GetMagic()==1) {
	float time=fHitChan->GetPulse(j)->GetTime();
	for (int k=0; k<fEvent->GetTotImpulses(); k++){
	  if (fEvent->GetImpulse(k)->GetTime()==time) {
	      noisePulseID=k;
	  }
	}
      }
    }
  }

  //detector-level impulses
  
  std::vector<int> string_impulses[8];

  int max_ampl_id[8] = {0};
  float max_ampl[8] = {0};
  float max_adjacent_amlp[8] = {0};
  
  int impulse_n=fEvent->GetTotImpulses();
  for (int j=0; j<impulse_n; j++){
    int iChannel=fEvent->GetImpulse(j)->GetChannelID();
    int iString=floor(iChannel/24)+1;
    string_impulses[iString].push_back(j);    
    if (fEvent->GetImpulse(j)->GetAmplitude()>max_ampl[iString]){
      max_ampl[iString]=fEvent->GetImpulse(j)->GetAmplitude();
      //max_adjacent_amlp[iString]=
    }
  }
  
  return kTRUE;
}

/*
Int_t BMyRecoChanSetter::PostProcess()
{
  return kTRUE;
}
*/
