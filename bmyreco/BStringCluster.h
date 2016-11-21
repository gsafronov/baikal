#ifndef BARS_BStringCluster
#define BARS_BStringCluster

#include <math.h>
#include <vector>
#include <cstddef>
#include "TObject.h"

class BEvent;
class BMCEvent;
class BGeomTel;
class BImpulse;

class BStringCluster : public TObject
{
 private:
  BEvent* fEvent;
  BMCEvent* fMCEvent;
  BGeomTel* fGeomTel;
  
  int fSize;
  int fString;
  float fCenterZ;
  float fCenterTime;
  float fHotSpotAmpl;
  float fSumAmpl;
  std::vector<int> fHotSpot;
  std::vector<int> fElements;
  
  bool fUseMCevent;
  
  //  std::vector<int> fMagicNumbers;
  
  std::vector<bool> fHasShowerHit;
  std::vector<bool> fHasDirectHit;
  std::vector<bool> fHasNoiseHit;
  
  std::vector<std::vector<int> > fTrackID; 
  
  std::pair<float,float> getClusterCenter(std::vector<int> stringCluster);
  
  float calculateSumAmpl();
  float calculateHotSpotAmpl();
  
 public:
  
  BStringCluster();
  BStringCluster(int iString, std::vector<int> hotspot, BEvent* event, BGeomTel* geomtel,  bool useMCevent = false, BMCEvent *mcevent = NULL);
  ~BStringCluster();
  
  int AddImpulse(int id);
  //  int RemoveImpulse(int id);
  int GetSize() {return fSize;}
  int GetStringID() {return fString;}
  float GetCenterZ() {return fCenterZ;}
  float GetCenterTime() {return fCenterTime;}
  float GetHotSpotAmpl() {return fHotSpotAmpl;}
  float GetSumAmpl() {return fSumAmpl;}
  
  BImpulse* GetConstituent(int id);
  int GetImpulseID(int id) {return fElements[id];}
  
  int GetHotSpotMax();
  int GetHotSpotMin();
  
  //  int GetNtracks(int id) {return fTrackID[id].size();}
  
  //  bool HasTrack(int id);
  
  ClassDef(BStringCluster,0);
    
};

#endif
