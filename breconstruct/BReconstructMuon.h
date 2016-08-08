#ifndef BARS_BReconstructMuon
#define BARS_BReconstructMuon

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// BReconstructMuon                                                        //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef MARS_MTask
#include "MTask.h"
#endif

#include "BReconstruct.h"
#include "Math/Interpolator.h"

class TH2F;
class TFile;

class BReconstructMuon : public BReconstruct
{
 private:
  const char* fQ10Filename;
  bool fHitProbabilityCriterionFlag;
  // Hit probability criterion
  double fRho_up_lim;
  ROOT::Math::Interpolator *fInterpolator;
  // MTask
  virtual Int_t   PreProcess(MParList *pList);
  virtual Int_t   Process();
  virtual Int_t   PostProcess();
  int HitProbabilityCriterion();
  TFile* fou;
  TH2F* hMapOfUsedOM;
  double GetProbability(double rho, double angle, bool hit);

 public:
  void SetQ10Path(const char* path) {fQ10Filename = path;};
  BReconstructMuon(const char *name = NULL, const char *title = NULL);
  void SetHitProbabilityCriterionFlag(bool flag=true) {
    fHitProbabilityCriterionFlag = flag; }

  ClassDef(BReconstructMuon, 0) // BReconstructMuon
};

#endif
