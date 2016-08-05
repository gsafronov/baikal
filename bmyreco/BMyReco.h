#ifndef BARS_BMyReco
#define BARS_BMyReco

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BMyReco
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


class BMyReco : public MTask
{
 protected:
  BEvent *     fEvent;             //! Generic event
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
  TH1F* hCharge;
  TH1F* hChargeNoLED;
  TH2F* hGeom;
  TH2F* hMapOfHits100;
  TH2F* hMapOfHits600;
  TH1F* hPairwiseSpeed;
  TH1F* hDR;
  TH1F* hDT;
  TH1F* hNPairs;
  TH1F* hSecNum;
  TH1F* hSecNumQual;
  TH1F* hClusterSize;
  TH1F* hClusterSizeSS;
  TH1F* hDirection;
  TH1F* hDeltaTvert;
  TH1F* hDeltaTvertNorm;
  TH1F* hDeltaT_DEBUG;
  TH1F* hDeltaTNorm_DEBUG;
  TH1F* hLLvert;
  TH1F* hLOverEst;
  TH1F* hLOverEstNoFac;
  TH1F* hUsedHitFrac;
  TH1F* hNTracksVsCut;
  TProfile* pLLvertVsCut;
  TH2F* h_LL3_nHits;
  TH2F* h_LL4_nHits;
  TH2F* h_LL5_nHits;
  
 public:
  BMyReco(string fname);
  ~BMyReco();
  
 private:
  //counter
  float clight;

  Int_t iEvent;
  std::vector<float> getHitLocation(int chanID);
  float getTimeEstimate_ns(std::vector<float> initialPoint, std::vector<float> trajectoryDirection, std::vector<float> pointInSpace);
  std::vector<int> filterHits(float minAmpl);
  std::vector<int> makeQualityCluster(std::vector<int> clusterIDs, int initialHitID);
  std::vector<int> makeQualitySSCluster(std::vector<int> clusterIDs, int initialHitID);
  float getLvalue(std::vector<int> clusterQual_id, std::vector<float> initialPoint, std::vector<float> trajectoryDirection);

  float getDeltaT_ns(int hit1_id, int hit2_id);
  float getDeltaR_m(int hit1_id, int hit2_id);
  std::vector<float> getDirection(int hit1_id, int hit2_id);


  ClassDef(BMyReco, 0);
};
    
#endif
