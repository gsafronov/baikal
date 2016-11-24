#include "BStringCluster.h"
#include "BEvent.h"
#include "BMCEvent.h"

ClassImp(BStringCluster);

BStringCluster::BStringCluster():
  fSignalShower(0), fSignalDirect(0), fNoise(0), fTrackID_hit(0), fTrackID_cluster(0), fElements(0), fNTracksPerHit(0)  
{
  fMCEvent=NULL;
  fGeomTel=NULL;
  
  fSize=0;
  fString=0;

  fCenterZ=0;
  fCenterTime=0;
  
  fHotSpotAmpl=0;
  fSumAmpl=0;

  fUseMCevent=false;
}

BStringCluster::~BStringCluster()
{
}


BStringCluster::BStringCluster(int iString, std::vector<int> hotspot, BEvent* event, BGeomTel* geomtel, bool useMCevent, BMCEvent* mcevent):
  fSignalShower(0), fSignalDirect(0), fNoise(0), fTrackID_hit(0), fTrackID_cluster(0), fNTracksPerHit(0)
{
  fEvent=event;
  fMCEvent=mcevent;
  fGeomTel=geomtel;
  
  fElements.push_back(hotspot[0]);
  fElements.push_back(hotspot[1]);

  fHotSpot.push_back(hotspot[0]);
  fHotSpot.push_back(hotspot[1]);
  
  fSize=fElements.size();
  fString=iString;

  std::pair<float, float> center=getClusterCenter(fElements);
  fCenterZ=center.first;
  fCenterTime=center.second;
  
  fHotSpotAmpl=calculateHotSpotAmpl();
  fSumAmpl=calculateSumAmpl();

  fUseMCevent=useMCevent;

  if (fUseMCevent){
    doMCMatching();
  }
  
}

int BStringCluster::AddImpulse(int id)
{
  fElements.push_back(id);
  fSize=fElements.size();
  std::pair<float, float> center=getClusterCenter(fElements);
  fCenterZ=center.first;
  fCenterTime=center.second;
  fSumAmpl=calculateSumAmpl();
  if (fUseMCevent){
    doMCMatching();
  }
  return 1;
}


float BStringCluster::calculateSumAmpl()
{
  float amplSum=0;
  for (int i=0; i<fElements.size(); i++){
    amplSum+=fEvent->GetImpulse(fElements[i])->GetAmplitude();
  }
  return amplSum;
}

float BStringCluster::calculateHotSpotAmpl()
{
  float amplSum=0;
  for (int i=0; i<fHotSpot.size(); i++){
    amplSum+=fEvent->GetImpulse(fHotSpot[i])->GetAmplitude();
  }
  return amplSum;
}

std::pair<float, float> BStringCluster::getClusterCenter(std::vector<int> stringCluster)
{
  float weightedSum=0;
  float amplSum=0;
  for (int i=0; i<stringCluster.size(); i++){
    float zPulse=fGeomTel->At(fEvent->GetImpulse(stringCluster[i])->GetChannelID())->GetZ();
    float amplPulse=fEvent->GetImpulse(stringCluster[i])->GetAmplitude();
    weightedSum+=zPulse*amplPulse;
    amplSum+=amplPulse;
  }
  
  if (weightedSum==0&&amplSum==0) return std::make_pair(0,0);
  float zCenter=weightedSum/amplSum;
  
  //initialise time with preliminary estimate from central pulse, then find closest to the center;
  float tCenter=fEvent->GetImpulse(stringCluster[int(floor(stringCluster.size()/2))])->GetTime();
  float zMin=1000;
  for (int i=0; i<stringCluster.size(); i++){
    float zPulse=fGeomTel->At(fEvent->GetImpulse(stringCluster[i])->GetChannelID())->GetZ();
    if (fabs(zPulse-zCenter)<zMin){
      zMin=fabs(zPulse-zCenter);
      tCenter=fEvent->GetImpulse(stringCluster[i])->GetTime();
    }
  }
  return std::make_pair(zCenter, tCenter);
}

BImpulse*  BStringCluster::GetConstituent(int id)
{
  return fEvent->GetImpulse(fElements[id]);
}
  
int BStringCluster::GetHotSpotMax()
{
  return max(fEvent->GetImpulse(fHotSpot[0])->GetAmplitude(),
	     fEvent->GetImpulse(fHotSpot[1])->GetAmplitude());
}

int BStringCluster::GetHotSpotMin()
{
  return min(fEvent->GetImpulse(fHotSpot[0])->GetAmplitude(),
	     fEvent->GetImpulse(fHotSpot[1])->GetAmplitude());
}


int BStringCluster::doMCMatching()
{
  int mc_n_channels=fMCEvent->GetChannelN();
  fNTracksPerCluster=0;
  int globalTrackCounter[100]={0};
  fNTracksPerHit.erase(fNTracksPerHit.begin(), fNTracksPerHit.end());
  //  fNTracksPerHit.clear();
  fTrackID_hit.clear();
  fTrackID_cluster.clear();
   
  for (int iPulse=0; iPulse<fSize; iPulse++){
    //   std::cout<<"NEXT PULSE"<<std::endl;
    int chID=fEvent->GetImpulse(fElements[iPulse])->GetChannelID();
    float time=fEvent->GetImpulse(fElements[iPulse])->GetTime();
    float noise=0;
    float direct=0;
    float shower=0;
    int trackIDs[100]={0};
    
    for (int iMCChan=0; iMCChan<mc_n_channels; iMCChan++){
      int mcchID=fMCEvent->GetHitChannel(iMCChan)->GetChannelID()-1;
      mcchID=24*floor(mcchID/24)+(24-mcchID%24);
      mcchID=mcchID-1;
      if (mcchID!=chID) continue;
      BMCHitChannel* mcChan=fMCEvent->GetHitChannel(iMCChan);
      int mc_n_pulses=mcChan->GetPulseN();
      
      for (int iMCPulse=0; iMCPulse<mc_n_pulses; iMCPulse++){
     	if (fabs(mcChan->GetPulse(iMCPulse)->GetTime()-time)<0.0001) {
	  int magic=mcChan->GetPulse(iMCPulse)->GetMagic();
	  //std::cout<<"yo!   "<<time<<"   "<<mcChan->GetPulse(iMCPulse)->GetTime()<<"  "<<magic<<std::endl;
	  if (magic==1) noise+=mcChan->GetPulse(iMCPulse)->GetAmplitude();
	  if (magic<0) {
	    direct+=mcChan->GetPulse(iMCPulse)->GetAmplitude();
	    trackIDs[magic+1000-1]++;
	    globalTrackCounter[magic+1000-1]++;
	  }
	  if (magic>=1000) {
	    shower+=mcChan->GetPulse(iMCPulse)->GetAmplitude();
	    trackIDs[magic%1000-1]++;
	    globalTrackCounter[magic%1000-1]++;
	  }
	}
      }
    }
    
    fSignalShower.push_back(shower);
    fSignalDirect.push_back(direct);
    fNoise.push_back(noise);
    std::vector<int> tracks;
    fNTracksPerHit.push_back(0);
    for (int i=0; i<100; i++){
      if (trackIDs[i]>=1){
	tracks.push_back(i);
	fNTracksPerHit[iPulse]++;
      }
    }
    //    if (tracks.size()>0) std::cout<<"bug: "<<chID<<"  "<<fMCEvent->GetResponseMuonsN()<<std::endl;
    fTrackID_hit.push_back(tracks);
  }

  for (int i=0; i<100; i++){
    if (globalTrackCounter[i]>=1) {
      fNTracksPerCluster++;
      fTrackID_cluster.push_back(i);
    }
  }
}

