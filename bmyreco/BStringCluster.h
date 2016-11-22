#ifndef BARS_BStringCluster
#define BARS_BStringCluster

#include <math.h>
#include <vector>
#include <cstddef>
#include "TObject.h"
#include "BGeomTel.h"

class BEvent;
class BMCEvent;
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
  
  std::vector<float> fSignalShower;
  std::vector<float> fSignalDirect;
  std::vector<float> fNoise;
  
  std::vector<std::vector<int> > fTrackID_hit;
  std::vector<int> fTrackID_cluster;
  int fNTracksPerCluster;
  
  std::pair<float,float> getClusterCenter(std::vector<int> stringCluster);
  
  float calculateSumAmpl();
  float calculateHotSpotAmpl();

  int doMCMatching();
  
 public:
  
  BStringCluster();
  BStringCluster(int iString, std::vector<int> hotspot, BEvent* event, BGeomTel* geomtel,  bool useMCevent = false, BMCEvent *mcevent = NULL);
  ~BStringCluster();
  
  int AddImpulse(int id);
  //  int RemoveImpulse(int id);
  int GetSize() {return fSize;}
  int GetStringID() {return fString;}
  float GetX() {return fGeomTel->At(fString*24)->GetX();}
  float GetY() {return fGeomTel->At(fString*24)->GetY();} 
  float GetCenterZ() {return fCenterZ;}
  float GetCenterTime() {return fCenterTime;}
  float GetHotSpotAmpl() {return fHotSpotAmpl;}
  float GetSumAmpl() {return fSumAmpl;}
  
  BImpulse* GetConstituent(int id);
  int GetImpulseID(int id) {return fElements[id];}
  
  int GetHotSpotMax();
  int GetHotSpotMin();

  int GetNSignalTracksHit(int id) {return fTrackID_hit[id].size();}
  std::vector<int> GetTracksHit(int id) {return fTrackID_hit[id];}
  
  float GetNoise(int id) {return fNoise[id];}
  float GetSignalDirect(int id) {return fSignalDirect[id];}
  float GetSignalShower(int id) {return fSignalShower[id];}

  int GetNSignalTracks() {return fNTracksPerCluster;}
  std::vector<int> GetTracks() {return fTrackID_cluster;}
  
  ClassDef(BStringCluster,0);
    
};

#endif
