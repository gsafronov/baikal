//////////////////////////////////////////////////////////////////////////////
//
//  BReconstructLaser
//
//
//////////////////////////////////////////////////////////////////////////////
#include "BReconstructLaser.h"


#include "MLog.h"
#include "MLogManip.h"

#include "MParList.h"
#include "TMath.h"
#include "TTimeStamp.h"

#include <fstream>
#include <iostream>

#include "BGeomTel.h"
#include "BCalibTel.h"
#include "BEvent.h"
#include "BEventMask.h"
#include "BRecParameters.h"
#include "BSecInteraction.h"

ClassImp(BReconstructLaser);

using namespace std;

//
//
//
BReconstructLaser::BReconstructLaser(const char *name, const char *title)
{
    fName  = name  ? name  : "BReconstructLaser";
    fTitle = title ? title : "ReconstructLaser";    
}

//
//
//
Int_t BReconstructLaser::PreProcess(MParList *pList)
{	
	if(!BReconstruct::PreProcess(pList)) {
		return kFALSE;
	}
    
    fSecInteraction = (BShower*)pList->FindCreateObj(AddSerialNumber("BShower"));
    if(!fSecInteraction) {
		*fLog << err << "Cannot create " << AddSerialNumber("BShower") << endl;
		return kFALSE;
    }
    
    // Set Reconstruction Parameters
/*    int n = 3;
    fSecInteraction->SetNumPar(n);
    
    for(int i = 0; i < n; i++) {
		fSecInteraction->SetIsFixed(i, 0);
		fSecInteraction->SetInitial(i, 0);
		fSecInteraction->SetMinLimit(i, 0);
		fSecInteraction->SetMaxLimit(i, 0);
		fSecInteraction->SetStep(i, 0.01);
	}
*/            
    return kTRUE;
}
