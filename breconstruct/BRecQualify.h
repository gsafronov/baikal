#ifndef BARS_BRecQualify
#define BARS_BRecQualify

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BRecQualify
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef MARS_MTask
#include "MTask.h"
#endif

class BSourceEAS;
class BEvent;
class BEventMask;
class BGeomTel;
class BCalibTel;
class BSecInteraction;
class BChannelMask;
class BRecParameters;

class BRecQualify : public MTask {

private:
	BSourceEAS*      fEventSource;
	BEvent*          fEvent;                 //! Generic event
	BEventMask*      fEventMask;             //! EventMask
	BGeomTel*        fGeomTel;               //! Geometry information
	BSecInteraction* fSecInteraction; //!	 
	BChannelMask*    fChannelMask;    //!
	 
	BRecParameters*  fRecParameters;
	BEventMask*      fMCEventMask;         // filter event mask	
	
	TString fFilterEventMaskName;
	
    Int_t   PreProcess(MParList *pList);
	Int_t   Process();
	Int_t   PostProcess();
	
	Bool_t fVerbose;
	Int_t  fNchIsRes;

public:
	BRecQualify(const char *name = NULL, const char *title = NULL);
	
	void SetFilterEventMaskName(TString name) { fFilterEventMaskName = name; }
	void MCResiduals();
	void MCResidualsX0Y0();
	void SetVerbose(Bool_t v) { fVerbose = v; }
	void SetNchIsRes(Int_t nch) { fNchIsRes = nch; }    
	
	ClassDef(BRecQualify, 0) // Task to 
};

#endif

