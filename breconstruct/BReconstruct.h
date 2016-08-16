#ifndef BARS_BReconstruct
#define BARS_BReconstruct

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BReconstruct
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef MARS_MTask
#include "MTask.h"
#endif

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Interpolator.h"

#include <vector>

#include <limits>

class BEvent;
class BEventMask;
class BGeomTel;
class BCalibTel;
class BSecInteraction;
class BRecParameters;
class BChannelMask;
/* class BSourceEAS; */
class BMuonCriterionParameters;

class BReconstruct : public MTask
{
 protected:
  BEvent *     fEvent;             //! Generic event
  BEventMask * fEventMask;         //! EventMask
  BGeomTel *   fGeomTel;           //! Geometry information
  BSecInteraction *   fSecInteraction;             //!
  BRecParameters *    fRecParameters;              //!
  BChannelMask *      fChannelMask;   //!
  BSecInteraction *   fZeroSecInteraction;
  /* BSourceEAS * fEventSource; */
  BMuonCriterionParameters * fMuonCriterionParams;
  /* BEventMask*      fMCEventMask; */

  TString fEventMaskName;

  Bool_t      fVerbose;
  //Bool_t      fIsWOAzimuth;  // use reconstruction without azimuth angle (usually for one string rec)
  Bool_t      fFixAzimuth;     // fix azimuth angle
  Double_t    fXshift;
  Double_t    fYshift;
  Double_t    fZshift;         // coordinate orgin shift along Z axis
  Int_t       fNchIsRes;       // channel to exclude from reconstruction
  int         fInitialsType;

  ROOT::Math::Minimizer * fMinimizer;
  ROOT::Math::Functor *   fFunctor;

  // MTask
  virtual Int_t   PreProcess(MParList  *pList);
  virtual Int_t   Process();
  virtual Int_t   PostProcess();

  virtual Int_t Minimize();
  virtual Double_t Chi2(const Double_t  *x);
  virtual Double_t M_estimator(const Double_t  *x);
  virtual void Residuals();

  BSecInteraction * ZeroRecAnalytic() const;
  void Tauinv(Double_t txf, Double_t tyf, Double_t tzf, Double_t &theta, Double_t &phi) const;
  void To_cos(Double_t the, Double_t phi, Double_t &t1, Double_t &t2, Double_t &f1, Double_t &f2) const;
  void InitialValues();

  struct pair {
    int first, second;
  };


 public:
  BReconstruct(const char * name = NULL, const char * title = NULL);
  ~BReconstruct();

  void SetVerbose(Bool_t verbose) { fVerbose = verbose; }
  //void SetWOAzimuth(Bool_t a) { fIsWOAzimuth = a; }
  void FixAzimuth(Bool_t a) { fFixAzimuth = a; }
  void SetNchIsRes(Int_t nch) { fNchIsRes = nch; }
  void SetInitialsType(Int_t type) {fInitialsType = type; }

  using MTask::Print;
  void Print();

  void SetEventMaskName(TString name) { fEventMaskName = name; }

  ClassDef(BReconstruct, 0) // Reconstruct
    };

#endif
