#ifndef BARS_BClusterProducer
#define BARS_BClusterProducer

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BTimeClusters
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

//#include "MTask.h"
#include "BFilter.h"
#include "BStringCluster.h"

class TH1F;
class TH2F;
class TProfile2D;
class TFile;
class TVector3;
class TProfile;
class BEvent;
class BMCEvent;
class BEventMask;
class BGeomTel;
class BChannelMask;
class BExtractedImpulseTel;
class BJoinExtractedImpulseTel;
class BExtractedHeader;
class BJoinExtractedHeader;
class BRecParameters;
//class BStringCluster;


class BClusterProducer : public BFilter
{
 private:
  BChannelMask *      fChannelMask;   //!
  BEventMask *        fInputEventMask;
  BMCEvent *          fMCEvent;
  BEvent *            fEvent;
  BGeomTel *          fGeomTel; 
  
  // MTask
  Int_t   PreProcess(MParList * pList);
  //Int_t         Process();
  virtual Int_t   PostProcess();
  Bool_t          Filter();

  TString     fInputMaskName;
  
  std::vector<BStringCluster> buildStringClusters(int iString, std::vector<int> string_impulses);
  int addImpulses(BStringCluster* hotspot, std::vector<int> string_impulses);

  float getTrackDistanceToOM(TVector3 initialPoint, TVector3 direction, TVector3 xyzOM);

  int runWPscan(std::vector<int>* string_impulses, std::vector<bool> noiseOMs_gen, std::vector<bool> signalOMs_gen, std::vector<bool> signalStrings_gen);
  
  int sgn(float x);

  int findTriggerString(float high_threshold, float low_threshold);
  
  int countMCStringGaps(std::vector<int> string_impulses_gen);
  
  std::vector<BStringCluster> findHotSpots(int iString, std::vector<int> string_impulses);

  std::pair<float,float> getPolarEstimate_string(BStringCluster stringCluster);

  std::vector<BStringCluster> buildGlobalCluster(std::vector<BStringCluster>* stringClusters);
  
  std::vector<int> chanIDs;
  bool fMaskNoise;
  
  float cVacuum;
  float cWater;

  float fSignalCut_gen;
  float fSeedSignalCut_gen;
  float fSignalCut_hotspot1;
  float fSignalCut_hotspot2;
  float fGen_rhoCut;
  float fSafetyWindow;
  float fGen_minAngle;
  float fGen_maxAngle;

  int fEventCounter;
  int fNoiseOMs;
  int fSignalOMs;
  int fSignalStrings;
  
  TFile* fOUT;

  TH1F* h_ntracks;
  TH1F* h_ntracks_6hits3strings;
  TH1F* h_debug;
  TH2F* h_gloClu_vs_firedStrings;
  TH1F* h_strClu_ntracks;
  TH1F* h_gloClu_ntracks;
  TH1F* h_pulse_ntracks;
  TH1F* h_strClu_noiseFrac;

  
 public:
  BClusterProducer(const char *name = NULL, const char *title = NULL);
  ~BClusterProducer();

  void        SetInputMaskName(TString name) {fInputMaskName = name; }
  void        SetOutputMaskName(TString name) {fOutputMaskName = name; }
  
  ClassDef(BClusterProducer, 0);
};
    
#endif
