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
  //virtual Int_t   PostProcess();
  Bool_t          Filter();

  TString     fInputMaskName;
  
  std::vector<int> BuildStringCluster(int iString, std::vector<int> string_impulses);
  std::vector<int> AddClusterImpulses(int max_ampl_id, int adjacent_id, std::vector<int> string_impulses, float window);
  
  std::vector<int> chanIDs;
  bool fMaskNoise;
  
  float cVacuum;
  float cWater;

  
 public:
  BTimeClusters(const char *name = NULL, const char *title = NULL);
  ~BTimeClusters();

  void        SetInputMaskName(TString name) {fInputMaskName = name; }
  void        SetOutputMaskName(TString name) {fOutputMaskName = name; }
  

  ClassDef(BTimeClusters, 0);
};
    
#endif
