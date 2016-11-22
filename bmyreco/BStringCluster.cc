#include "BStringCluster.h"
#include "BEvent.h"
#include "BMCEvent.h"
#include "BGeomTel.h"

ClassImp(BStringCluster);

BStringCluster::BStringCluster():
  fHasShowerHit(0), fHasDirectHit(0), fHasNoiseHit(0), fTrackID(0), fElements(0)  
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
  fHasShowerHit(0), fHasDirectHit(0), fHasNoiseHit(0), fTrackID(0) 
{
  fEvent=event;
  fMCEvent=mcevent;
  fGeomTel=geomtel;
  
  fElements.push_back(hotspot[0]);
  fElements.push_back(hotspot[1]);
  
  fSize=fElements.size();
  fString=iString;

  std::pair<float, float> center=getClusterCenter(fElements);
  fCenterZ=center.first;
  fCenterTime=center.second;
  
  fHotSpotAmpl=calculateHotSpotAmpl();
  fSumAmpl=calculateSumAmpl();

  fUseMCevent=useMCevent;

  /*
  if (fUseMCevent){
    for (int i=0; i<hotspot.size(); i++){
      fHasShowerHit=matchShowers(hotspot[i]);
      fHasDirectHit=matchDirect(hotspot[i]);
      hHasNoiseHit=matchNoise(hotspot[i]);
      
      fTrackID.push_back(getMatchedTrackIDs(hotspot[i])); 
    }
  }
  */
}

int BStringCluster::AddImpulse(int id)
{
  fElements.push_back(id);
  std::pair<float, float> center=getClusterCenter(fElements);
  fCenterZ=center.first;
  fCenterTime=center.second;
  fSumAmpl=calculateSumAmpl();
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
  for (int i=0; i<fElements.size(); i++){
    amplSum+=fEvent->GetImpulse(fElements[i])->GetAmplitude();
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


