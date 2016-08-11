#ifndef BARS_BMyRecoChanSetter
#define BARS_BMyRecoChanSetter

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BMyRecoChanSetter
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef MARS_MTask
#include "MTask.h"
#endif

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



class BMyRecoChanSetter : public MTask
{
 protected:
  BChannelMask *      fChannelMask;   //!
  
  // MTask
  virtual Int_t   PreProcess(MParList * pList);
  virtual Int_t   Process();
  virtual Int_t   PostProcess();
  
 public:
  BMyRecoChanSetter(int ch);
  ~BMyRecoChanSetter();
  
 private:
  int chanID;

  ClassDef(BMyRecoChanSetter, 0);
};
    
#endif
