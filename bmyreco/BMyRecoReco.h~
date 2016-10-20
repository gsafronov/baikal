#ifndef BARS_BMyRecoReco
#define BARS_BMyRecoReco

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BMyRecoReco
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
class BMCEvent;
class BEventMask;
class BGeomTel;
class BChannelMask;
class BExtractedImpulseTel;
class BJoinExtractedImpulseTel;
class BExtractedHeader;
class BJoinExtractedHeader;
class BRecParameters;



class BMyRecoReco : public MTask
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
  BRecParameters *    fRecParam;              //!
  BChannelMask *      fChannelMask;   //!
  //  BSecInteraction *   fZeroSecInteraction;
  
  // MTask
  virtual Int_t   PreProcess(MParList * pList);
  virtual Int_t   Process();
  virtual Int_t   PostProcess();
  
 public:
  BMyRecoReco(string fname, int cChan, bool useMC);
  ~BMyRecoReco();
  
 private:
  //counter
  int iEvent;
  float cWater;
  float cVacuum;

  TFile* fOUT;
  TH1F* hMuMult;
  TH1F* hPolar1muMC;
  TH1F* hPolar1muRec;
  TH1F* hRho1muMC;
  TH1F* hNhitMC;
  TH1F* hNhit;
  TH1F* hNdirectHit;
  TH1F* hNshowerHit;
  TH1F* hNnoiseHit;
  TH1F* hDistToTrackMC;
  TH1F* hNhitPerChannel;
  
  TH1F* hNhit1mu;
  TH1F* hChi2;
  TH1F* hChi2_zoom1;
  TH2F* hChi2_rhoAvg;
  TH2F* hChi2_rhoMax;
  TH1F* hTimeDiff;
  TH1F* hTimeDiff_usedHits;
  TH1F* hDistToTrack;
  TH1F* hDistToTrack_usedHits;
  TH1F* hDistToTrack1mu_usedHits;

  TH1F* hDistToTrack_debug;
  
  TH1F* hAngle1muGenRec;
  TH1F* hMagNum_usedHits;

  TH1F* hHitsDR_MC;

  TH2F* hExtrapolatedHits;
  TH2F* hRealHits;

  TH2F* hExtrapolatedHitsALL;
  TH2F* hRealHitsALL;

  TH2F* hNumMC_NumReco;
  TH2F* hDMC_dReco;
  
  //for MC debug

  
  //  float getTrackDistanceToOM(std::vector<float> initialPoint, std::vector<float> trajectoryDirection, std::vector<float> xyzOM);
  float getTrackDistanceToOM(TVector3 initialPoint, TVector3 trajectoryDirection, TVector3 xyzOM);
  float getTimeEstimate_ns(TVector3 initialPoint, TVector3 trajectoryDirection, TVector3 pointInSpace);
  int calibChanID;

  int RunMCAnalysis();
  
  TVector3 xyToXYZ(TVector3 s, float X, float Y);
  
  ClassDef(BMyRecoReco, 0);
};
    
#endif
