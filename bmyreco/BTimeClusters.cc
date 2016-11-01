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

  cVacuum=3e8;
  cWater=3e8/1.33;
  
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
  
  fGeomTel = (BGeomTel*)pList->FindObject(AddSerialNumber("BGeomTel"));

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
  /*
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
  */
  //detector-level impulses
  
  std::vector<int> string_impulses[8];

  int impulse_n=fEvent->GetTotImpulses();
  for (int iPulse=0; iPulse<impulse_n; iPulse++){
    int iChannel=fEvent->GetImpulse(iPulse)->GetChannelID();
    string_impulses[int(floor(iChannel/24))].push_back(iPulse);    
  }

  for (int iString=0; iString<8; iString++){
    std::vector<int> stringCluster=BuildStringCluster(iString, string_impulses[iString]);
  }
  
  
  return kTRUE;
}

std::vector<int> BTimeClusters::BuildStringCluster(int iString, std::vector<int> string_impulses)
{
  std::vector<int> stringCluster;

  //find seed channel - the one with maximum signal
  float max_ampl=-1;
  int max_ampl_id=-1;
  
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
    if (max_ampl<fEvent->GetImpulse(string_impulses[iPulse])->GetAmplitude()){
      max_ampl=fEvent->GetImpulse(string_impulses[iPulse])->GetAmplitude();
      max_ampl_id=iPulse;
    }
  }

  //check signal in seed neighbours
  float adjacent_ampl=0;
  int adjacent_id=-1;
  int seed_channel_id=fEvent->GetImpulse(string_impulses[max_ampl_id])->GetChannelID();
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
    int channel_id=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
    if (!(channel_id==seed_channel_id+1||channel_id==seed_channel_id-1)) continue;
    if (adjacent_ampl<fEvent->GetImpulse(string_impulses[iPulse])->GetAmplitude()){
      adjacent_ampl=fEvent->GetImpulse(string_impulses[iPulse])->GetAmplitude();
      adjacent_id=iPulse;
    }
  }

  //fail if no adjacent channel with ampl > 2 pe. But need to check other seeds then. Possibly build pairs of adjacent channels with A > 2pe and treat the one with higher signal as seed.
  
  if (adjacent_ampl<2) return stringCluster;

  //hot spot is found, run cluster buiding procedure
  float window = 50;
  stringCluster=AddClusterImpulses(max_ampl_id, adjacent_id, string_impulses, window);

  return stringCluster;
}

std::vector<int> BTimeClusters::AddClusterImpulses(int max_ampl_id, int adjacent_id, std::vector<int> string_impulses, float window)
{
  std::vector<int> stringCluster;
  stringCluster.push_back(string_impulses[max_ampl_id]);
  stringCluster.push_back(string_impulses[adjacent_id]);

  float seed_time=fEvent->GetImpulse(string_impulses[max_ampl_id])->GetTime();
  int seed_channel_id=fEvent->GetImpulse(string_impulses[max_ampl_id])->GetChannelID();
  float adjacent_time=fEvent->GetImpulse(string_impulses[adjacent_id])->GetTime();
  int adjacent_channel_id=fEvent->GetImpulse(string_impulses[adjacent_id])->GetChannelID();

  float deltaZ=fGeomTel->At(adjacent_channel_id)->GetZ()-fGeomTel->At(seed_channel_id)->GetZ();
  float deltaT=adjacent_time-seed_time;
  int id_increment=adjacent_channel_id-seed_channel_id;

  //add +-1 channels
  float time_early=min(seed_time-(deltaT/fabs(deltaT))*fabs(deltaZ)/cWater-window, seed_time-deltaT-window); //experimental
  float time_late=max(seed_time-(deltaT/fabs(deltaT))*fabs(deltaZ)/cWater+window, seed_time-deltaT+window);
  
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
    if (iPulse==max_ampl_id||iPulse==adjacent_id) continue;
    int chanID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
    if (chanID==seed_channel_id-id_increment){
      float chan_time=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
      if (chan_time>time_early&&chan_time<time_late) stringCluster.push_back(string_impulses[iPulse]);
    }
  }

  return stringCluster;  
}

/*
Int_t BMyRecoChanSetter::PostProcess()
{
  return kTRUE;
}
*/
