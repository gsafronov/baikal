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
  
  std::vector<int> BuildStringCluster(int iString, std::vector<int> string_impulses);
  std::vector<int> AddClusterImpulses(int max_ampl_id, int adjacent_id, std::vector<int> string_impulses, float window);

  float getTrackDistanceToOM(TVector3 initialPoint, TVector3 direction, TVector3 xyzOM);

  int runWPscan(std::vector<int>* string_impulses, std::vector<int> noiseOMs_gen, std::vector<int> signalOMs_gen);
  
  std::vector<int> chanIDs;
  bool fMaskNoise;
  
  float cVacuum;
  float cWater;

  float fSignalCut_gen;
  float fSeedSignalCut_gen;
  float fSignalCut_hotspot;
  float fGen_rhoCut;
  float fTimeMargin;
  float fGen_minAngle;
  float fGen_maxAngle;

  int fEventCounter;
  int fNoiseOMs;
  int fSignalOMs;
  
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

  TH2F* h_clusteredFrac_noise;
  TH2F* h_clusteredFrac_signal;

  TH2F* h_effpur_mult;
  
  TH2F* h_clustered;

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
