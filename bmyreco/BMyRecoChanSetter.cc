#include "BMyRecoChanSetter.h"
#include "BGeomTel.h"
#include "BEvent.h"
#include "BMCEvent.h"
#include "BExtractedImpulse.h"
#include "BExtractedImpulseTel.h"
#include "BJoinExtractedImpulseTel.h"
#include "BEventMask.h"
#include "BChannelMask.h"
#include "BExtractedHeader.h"

#include "BRecParameters.h"

#include "MParList.h"
#include "MLog.h"
#include "MLogManip.h"
 
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TMath.h"

#include <math.h>
#include <vector>

ClassImp(BMyRecoChanSetter);

BMyRecoChanSetter::BMyRecoChanSetter(int ch)
{
  chanID=ch;
}

BMyRecoChanSetter::~BMyRecoChanSetter()
{
}

Int_t BMyRecoChanSetter::PreProcess(MParList * pList)
{
  //  std::cout<<"WE ARE IN BMyRecoBEvt::PreProcess"<<std::endl;
   
  fChannelMask = (BChannelMask*)pList->FindObject("BChannelMask");
  if (!fChannelMask)
    {
      * fLog << err << AddSerialNumber("BChannelMask") << " not found... aborting." << endl;
      return kFALSE;
    }

  std::cout<<"BChanSetter: exclude calibrated channel from reconstruction"<<std::endl;
  if (fChannelMask->GetFlag(chanID)==1) fChannelMask->SetFlag(chanID, BChannelMask::EMode::kOff);
  
  return kTRUE;
}


Int_t BMyRecoChanSetter::Process()
{
  return kTRUE;
}

Int_t BMyRecoChanSetter::PostProcess()
{
  return kTRUE;
}
