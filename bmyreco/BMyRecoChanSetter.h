#ifndef BARS_BMyRecoChanSetter
#define BARS_BMyRecoChanSetter

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BMyRecoChanSetter
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


class BMyRecoChanSetter : public BFilter
{
 private:
  BChannelMask *      fChannelMask;   //!
  BEventMask *        fInputEventMask;
  BMCEvent *          fMCEvent;
  BEvent *            fEvent;
  
  // MTask
  Int_t   PreProcess(MParList * pList);
  //Int_t         Process();
  //virtual Int_t   PostProcess();
  Bool_t          Filter();

  TString     fInputMaskName;
  
  
  std::vector<int> chanIDs;
  bool fMaskNoise;
  
 public:
  BMyRecoChanSetter(const char *name = NULL, const char *title = NULL);
  ~BMyRecoChanSetter();

  void        SetInputMaskName(TString name) {fInputMaskName = name; }
  void        SetOutputMaskName(TString name) {fOutputMaskName = name; }
  

  ClassDef(BMyRecoChanSetter, 0);
};
    
#endif
