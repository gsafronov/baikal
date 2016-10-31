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
class TVector3;
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
  TH1F* hNmuResp;
  TH1F* hNmuTotal;
  TH1F* hPolar1muMC;
  TH1F* hRho1muMC;
  TH1F* hNhit1muMC;
  TH1F* hNdirect1muHit;
  TH1F* hNshower1muHit;
  TH1F* hNnoiseHit;
  TH1F* hDistToTrackMC_direct;
  TH1F* hDistToTrackMC_shower;
  TH1F* hNhitPerChannel;
  TH1F* hTimeDifference_norm;
  TH2F* hDiffReal_DiffEst;
  TH1F* hPulsesPerChannel;
  TH1F* hYieldVsAmpl_ref;

  TH1F* hChanOffset[192];
  TH1F* hClusterOffsets;
  TH1F* hChannelSignalLength;

  
  TH2F* hDebugTimeShowerHits;  
  TH2F* hSeedChanID_dr_vs_response;
  TH1F* hSeedChanID_dr_minus_response;
  TH2F* hPropagation_z_vs_time;
  TH2F* hDelaySeed_z_vs_delay;
  TH2F* hDelaySeedAnalytic_z_vs_delay;
  TH2F* hStringDelaySeed_z_vs_delay[8];
  TH2F* hStringDelaySeedAn_z_vs_delay[8];
  
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
  TH1F* hPulseTime;

  
 public:
  BMyRecoMC(string fname);
  ~BMyRecoMC();
  
 private:
  //  TVector3 fGenVec;
  //  TVector3 fInitialPoint;
  //  float fPrimThetaRad;
  //  float fPrimPhiRad;
  
  std::vector<int> filterHits(float minAmpl);
  std::vector<float> getHitLocation(int hitID);
  float getDeltaT_ns(int hit1_id, int hit2_id);
  float getDeltaR_m(int hit1_id, int hit2_id);
  std::vector<float> getDirection(int hit1_id, int hit2_id);
  
  int RunMCAnalysis();
  int RunClusteringAnalysis();
  float getTrackDistanceToOM(TVector3 initialPoint, TVector3 trajectoryDirection, TVector3 xyzOM);
  TVector3 getTrackClosestApproach(TVector3 initialPoint, TVector3 trajectoryDirection, TVector3 xyzOM);
  float getTimeEstimate_ns(TVector3 initialPoint, TVector3 trajectoryDirection, TVector3 pointInSpace);

  //counter
  float cWater;
  float cVacuum;

  Int_t iEvent;
  Int_t fNCalibTracks;

  ClassDef(BMyRecoMC, 0);
};
    
#endif
