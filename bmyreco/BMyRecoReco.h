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
  BMyRecoReco(string fname, bool useMC);
  ~BMyRecoReco();
  
 private:
  //counter
  int iEvent;
  float cWater;
  float cVacuum;

  TFile* fOUT;

  ClassDef(BMyRecoReco, 0);
};
    
#endif
