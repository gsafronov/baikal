#include "BMyRecoReco.h"
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

ClassImp(BMyRecoReco);

BMyRecoReco::BMyRecoReco(string fname, bool useMCEvent)
{
  iEvent=0;
 
  fOUT=new TFile(fname.c_str(),"RECREATE");

  cWater=3e8*pow(1.33,-1);
  cVacuum=3e8;
}

BMyRecoReco::~BMyRecoReco()
{
}

Int_t BMyRecoReco::PreProcess(MParList * pList)
{
  std::cout<<"WE ARE IN BMyRecoBEvt::PreProcess"<<std::endl;
  
  fGeomTel = (BGeomTel*)pList->FindObject(AddSerialNumber("BGeomTel"));
  if (!fGeomTel){
    * fLog << err << AddSerialNumber("BGeomTel") << " not found... aborting." << endl;
    return kFALSE;
  }

  for (int i=0; i<192; i++){
    std::cout<<"iChan: "<<i<<"   X: "<<fGeomTel->At(i)->GetX()<<"   Y: "<<fGeomTel->At(i)->GetY()<<"   Z: "<<fGeomTel->At(i)->GetZ()<<std::endl;
  }
  
  fRecParam = (BRecParameters*)pList->FindObject("BRecParameters");
  //if(!fRecParam){ * fLog << err << AddSerialNumber("BRecParameters") << " not found... aborting." << endl;
  //    return kFALSE;
  //  }
  return kTRUE;
}


Int_t BMyRecoReco::Process()
{
  iEvent++;
  if (iEvent%10000==0) std::cout<<"eventNumber: "<<iEvent<<std::endl;

  std::cout<<fRecParam->GetNumPar()<<std::endl;
  return kTRUE;
}


Int_t BMyRecoReco::PostProcess()
{
  return kTRUE;
}
