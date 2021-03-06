#ifndef BARS_BMyRecoBEvt
#define BARS_BMyRecoBEvt

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BMyRecoBEvt
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
class BMCEvent;
class BEventMask;
class BGeomTel;
class BChannelMask;
class BExtractedImpulseTel;
class BJoinExtractedImpulseTel;
class BExtractedHeader;
class BJoinExtractedHeader;


class BMyRecoBEvt : public MTask
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

  //BEvent plots
  TFile* fOUT;
  TH1F* hCharge;
  TH1F* hChargeNoLED;
  TH1F* hChargeOverAmpl;
  TH1F* hTime;
  TH1F* hTimeOverCharge;
  TH2F* hGeom;
  TH2F* hMapOfHits1pe;
  TH2F* hMapOfHitsNpe;
  TH2F* hMapOfClusteredHits;
  TH1F* hPairwiseSpeed;
  TH1F* hPairwiseSpeedMu;
  TH1F* hNClusters;
  TH1F* hSecNum;
  TH1F* hSecNumQual;
  TH1F* hClusterSize;
  TH1F* hClusterSizeSS;
  TH1F* hClustersPerHit;
  TH1F* hDirection;
  TH1F* hDeltaT_verticalMu;
  TH1F* hDeltaT_verticalMu_Norm;
  TH1F* hDeltaT_DEBUG;
  TH1F* hDeltaTNorm_DEBUG;
  TH1F* hLL_verticalMu;
  TH1F* hTime_overestimation;

  //plots using BMCEvent 
  TH1F* hThetaPrimary;
  TH1F* hMuonDelay;
  TH2F* hNHits_EnergyPrimary;
  TH2F* hNHits_NRespMuon;
  TH1F* hNRespMuons;
  TH1F* hNRespMuons_vertical;
  TH1F* hNMuonsTot;
  TH1F* hNMuonsTot_vertical;
  TH1F* hNRespMuonsBelow;
  TH2F* hClustersRespMatrix;
  TH2F* hNClusters_EnergyPrimary;
  TH2F* hNClusters_NRespMuon;
  TH1F* hSumPolarAngle;
  TH1F* hAngleRes;
  TH1F* hDirDebug;
  TH1F* hPolarVertical;
  TH1F* hPolarAll;
  TH1F* hDeltaT_verticalMuGen;
  TH1F* hDeltaT_verticalMuDelta;
  TH1F* hDeltaZ_vert;
  TH1F* hNPulsePerChannel;

  TH2F* hModMH_NRespMuons;
  
  //array of vector of histograms for calibration
  //  std::vector<TH1F> hDeltaT_pairChan[23];

 public:
  BMyRecoBEvt(string fname, bool useMC);
  ~BMyRecoBEvt();
  
 private:
  //counter
  float cWater;
  float cVacuum;

  bool fUseMCEvent;

  std::vector<bool> fHitIsUsed;

  Int_t iEvent;
  std::vector<float> getHitLocation(int chanID);
  float getTimeEstimate_ns(std::vector<float> initialPoint, std::vector<float> trajectoryDirection, std::vector<float> pointInSpace);
  std::vector<int> filterHits(float minAmpl);
  std::vector<int> makeQualityCluster(std::vector<int> clusterIDs, int initialHitID);
  std::vector<int> makeQualitySSCluster(std::vector<int> clusterIDs, int initialHitID);
  float getLvalue(std::vector<int> clusterQual_id, std::vector<float> initialPoint, std::vector<float> trajectoryDirection);

  int PlotClusters(std::vector<int> ff);

  float getDeltaT_ns(int hit1_id, int hit2_id);
  float getDeltaR_m(int hit1_id, int hit2_id);
  std::vector<float> getDirection(int hit1_id, int hit2_id);
  std::vector<float> GetPrimaryDirection();


  void SetUseMCEvent();

  int iLastBins;

  ClassDef(BMyRecoBEvt, 0);
};
    
#endif
