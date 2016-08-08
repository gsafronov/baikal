#ifndef BARS_BRecQuality
#define BARS_BRecQuality

#ifndef MARS_MParContainer
#include "MParContainer.h"
#endif

class BRecQuality : public MParContainer {

protected:
    Float_t fSourceAzimuthAngle;   // source azimuth angle in rad
    Float_t fSourcePolarAngle;     // source polar angle in rad

    Float_t fRecAzimuthAngle;      // rec azimuth angle in rad
    Float_t fRecPolarAngle;        // rec polar angle in rad
    
    Float_t fPsi;                  // angle between rec and source directions in rad

	void CalcPsi();

public:
	BRecQuality(const char *name = NULL, const char *title = NULL);
	
	void SetSourceAzimuthAngle(Float_t ang) { fSourceAzimuthAngle = ang; }
	void SetSourcePolarAngle(Float_t ang) { fSourcePolarAngle = ang; }
	
	void SetRecAzimuthAngle(Float_t ang);
	void SetRecPolarAngle(Float_t ang);

	virtual void   SetReadyToSave(Bool_t flag = kTRUE);

	ClassDef(BRecQuality, 1) // Base
};

#endif
