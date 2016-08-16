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
#include "Math/Interpolator.h"

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

  double fX;
  double fY;
  double fZ;
  double fTheta;
  double fPhi;
  double fTelLength;

  ROOT::Math::Interpolator *fInterpolator;
  int fMuonNumber;
  const char* fQ10FileName;
  bool fHitCriterionFlag;
  double fRhoUpLim;
  int HitCriterion();
  double GetProbability(double rho, double energy, double angle, bool hit);

  double    fXshift;
  double    fYshift;
  double    fZshift;         // coordinate orgin shift along Z axis

  void ZDist();
  double GetZDist(double theta, double phi);
  Int_t GetLargestMuonNum() const;
  double GetP(double energy);

 public:
	BRecQualify(const char *name = NULL, const char *title = NULL);

  void SetQ10Path(const char* path) {fQ10FileName = path;};
  void SetHitCriterionFlag(bool flag=true) {fHitCriterionFlag = flag;}
	void SetFilterEventMaskName(TString name) { fFilterEventMaskName = name; }
	void MCResiduals();
	void MCResidualsX0Y0();
	void SetVerbose(Bool_t v) { fVerbose = v; }
	void SetNchIsRes(Int_t nch) { fNchIsRes = nch; }

	ClassDef(BRecQualify, 0) // Task to
};

#endif
