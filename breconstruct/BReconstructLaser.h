#ifndef BARS_BReconstructLaser
#define BARS_BReconstructLaser

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BReconstructLaser                                                       //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef MARS_MTask
#include "MTask.h"
#endif

#include "BReconstruct.h"

class BReconstructLaser : public BReconstruct
{
	private:
	
		// MTask
		Int_t   PreProcess(MParList *pList);
	
	public:
		BReconstructLaser(const char *name = NULL, const char *title = NULL);
		
		ClassDef(BReconstructLaser, 0) // BReconstructLaser
};

#endif


