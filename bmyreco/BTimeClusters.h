#ifndef BARS_BTimeClusters
#define BARS_BTimeClusters

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BTimeClusters
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

//#include "MTask.h"
#include "BFilter.h"

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


class BTimeClusters : public BFilter
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
  
  std::vector<std::vector<int> > buildStringClusters(int iString, std::vector<int> string_impulses);
  std::vector<int> addImpulses(std::vector<int> hotspot, std::vector<int> string_impulses);

  float getTrackDistanceToOM(TVector3 initialPoint, TVector3 direction, TVector3 xyzOM);

  int runWPscan(std::vector<int>* string_impulses, std::vector<bool> noiseOMs_gen, std::vector<bool> signalOMs_gen, std::vector<bool> signalStrings_gen);
  
  int sgn(float x);

  int findTriggerString(float high_threshold, float low_threshold);
  
  int countMCStringGaps(std::vector<int> string_impulses_gen);
  
  std::vector<std::vector<int>> findHotSpots(std::vector<int> string_impulses);

  std::pair<float,float> getPolarEstimate_string(std::vector<int> stringCluster);

  std::pair<float,float> getClusterCenter(std::vector<int> stringCluster);

  std::vector<int> buildGlobalCluster(std::vector<std::vector<int> >* stringClusters);
  
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
  TH1F* h_hits_per_string_2pe;
  TH1F* h_hits_2pe;
  TH1F* h_fired_strings_2pe;
  
  TH1F* h_1mu_hits_per_string_2pe;
  TH1F* h_1mu_hits_2pe;
  TH1F* h_1mu_fired_strings_2pe;

  TH1F* hreco_hits_per_string;
  TH1F* hreco_fired_strings;
  TH1F* hreco_hits;

  TH2F* h_strings_reco_vs_gen;
  TH2F* h_strings_bevt_vs_gen;

  TH2F* h_hitsPerString_reco_vs_gen;
  TH2F* h_hitsPerString_bevt_vs_gen;

  TH2F* h_smuons_polar_vs_rho;
  TH2F* h_rcand_polar_vs_rho;
  TH2F* h_rcandFrac_polar_vs_rho;
  TProfile2D* h_strings_polar_vs_rho;
  
  TH1F* h_muon_energy;
  TH1F* h_muon_energy_rcand_gen;
  TH1F* h_muon_energy_rcand;

  TH2F* h_hitSignal_reco_vs_gen;

  TH2F* h_noiseFrac_clustered;
  TH2F* h_signalFrac_clustered;
  TH2F* h_stringSignalFrac_clustered;
  
  TH2F* h_clusteredFrac_noise;
  TH2F* h_clusteredFrac_signal;

  TH2F* h_effpur_mult;
  
  TH2F* h_clustered;
  TH2F* h_clustered_highMult;

  TH1F* hWhatIsOM;

  TH2F* hitMap_gen;
  TH2F* hitMap_det;
  
  //TH2F* h_spots_per_string;
  
 public:
  BTimeClusters(const char *name = NULL, const char *title = NULL);
  ~BTimeClusters();

  void        SetInputMaskName(TString name) {fInputMaskName = name; }
  void        SetOutputMaskName(TString name) {fOutputMaskName = name; }
  
  ClassDef(BTimeClusters, 0);
};
    
#endif
