#ifndef BARS_BMyRecoMC
#define BARS_BMyRecoMC

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BMyRecoMC
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef MARS_MTask
#include "MTask.h"
#endif

class TFile;
class TH1F;
class TH2F;
class TProfile;
class BEvent;
class BEventMask;
class BGeomTel;
class BChannelMask;
class BExtractedImpulseTel;
class BJoinExtractedImpulseTel;
class BExtractedHeader;
class BJoinExtractedHeader;
class BMCEvent;

class BMyRecoMC : public MTask
{
 protected:
  BEvent *     fEvent;             //! Generic event
  BMCEvent*    fMCEvent;           
  BEventMask * fEventMask;         //! EventMask
  BGeomTel *   fGeomTel;           //! Geometry information
  BExtractedImpulseTel* fExtractedImpulse;
  BJoinExtractedImpulseTel* fJoinExtractedImpulse;
  BExtractedHeader* fExtractedHeader;
  BJoinExtractedHeader* fJoinExtractedHeader;
  //  BMyRecParam *    fMyRecParam;              //!
  BChannelMask *      fChannelMask;   //!
  //  BSecInteraction *   fZeroSecInteraction;
  
  // MTask
  virtual Int_t   PreProcess(MParList * pList);
  virtual Int_t   Process();
  virtual Int_t   PostProcess();

  //output
  TFile* fOUT;
  TH1F* hMuonN;
  TH1F* hTotalMuonN;
  TH1F* hMuCutN;
  TH1F* hMuNotNull;
  TH1F* hMuonE;

  TH1F* hThetaPrimary;
  TH2F* hNHits_ThetaPrimary;
  TH2F* hNHits_EnergyPrimary;

  TH1F* hChannelN;

  TH1F* hChannelN_bevt;
  TH1F* hPulseN_bevt;
  TH2F* hMapOfHits1;
  TH2F* hMapOfHits3;  
  TH2F* hGeom;

 public:
  BMyRecoMC(string fname);
  ~BMyRecoMC();
  
  std::vector<int> filterHits(float minAmpl);
  std::vector<float> getHitLocation(int hitID);
  float getDeltaT_ns(int hit1_id, int hit2_id);
  float getDeltaR_m(int hit1_id, int hit2_id);
  std::vector<float> getDirection(int hit1_id, int hit2_id);
  


 private:
  //counter
  float clight;
  float cmu;

  Int_t iEvent;

  ClassDef(BMyRecoMC, 0);
};
    
#endif
