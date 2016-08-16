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

class BReconstructMuon : public BReconstruct
{
 private:
  // MTask
  virtual Int_t   PreProcess(MParList *pList);
  virtual Int_t   Process();
  void GetX0Y0(Double_t &x0, Double_t &y0, Double_t theta, Double_t phi) const;
  BSecInteraction * PlaneWaveSingleMuonDirectionReconstruction();
  void Combinations(int pool, int need, unsigned long chosen, int at,
                    vector<pair> & set_of_pairs, double * parameters);
  int InitialValues();
  int Chi2Fit();

 public:
  BReconstructMuon(const char *name = NULL, const char *title = NULL);

  ClassDef(BReconstructMuon, 0) // BReconstructMuon
};

#endif
