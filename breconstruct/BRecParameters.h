#ifndef BRECPARAMETERS
#define BRECPARAMETERS

#include <MParContainer.h>

//_____________________________________________________________________________________
class BRecParameters : public MParContainer {
 protected:
  Double_t     fFuncValue;          // Значение функции минимизации
  Int_t        fQual;               // 0 - поиск минимума прошел удачно, остальные значения - с ошибкой
  Int_t        fNumPar;
  Int_t*       fIsFixed;            //[fNumPar]
  Double_t*    fInitial;            //[fNumPar] Начальные значения варьируемых параметров
  Double_t*    fMinLimit;           //[fNumPar] Ограничение снизу на значения варьируемых параметров
  Double_t*    fMaxLimit;           //[fNumPar] Ограничение сверху на значения варьируемых параметров
  Double_t*    fStep;               //[fNumPar] Начальный шаг по варьируемым параметрам
  Double_t*    fErrors;             //[fNumPar] Ошибки по варьируемым параметрам

  Int_t          fNhit;             // Число значений (каналов или импульсов), по которым была произведена минимизация
  Double_t      *fTres;             //[fNhit] Резидуалы в нс
  Int_t         *fNchGeom;          //[fNhit] порядковые номер каналов

  Int_t          fNhitMC;           // Число каналов, для которых были посчитаны МК резидуалы от мюона
  Double_t      *fTresMC;           //[fNhitMC] МК резидуалы в нс
  Int_t         *fNchGeomMC;        //[fNhitMC] МК порядковые номер каналов

  Double_t       fThetaRec;         // reconstructed zenith angle in degrees in (0-180)
  Double_t       fThetaMC;          // MC zenith angle in degrees in (0-180)

  Double_t       fPhiRec;           // reconstructed azimuth angle in degrees in (0-180)
  Double_t       fPhiMC;            // MC azimuth angle in degrees in (0-180)

  Double_t       fX0Rec;         // reconstructed x0 in meters
  Double_t       fX0MC;          // MC x0 in meters

  Double_t       fY0Rec;         // reconstructed y0 in meters
  Double_t       fY0MC;          // MC y0 in meters

  Double_t       fTimeRec;          // Rec time in ns
  Double_t       fTimeMC;           // MC time in ns

  double      fHit_p;
  double      fNo_hit_p;

  Double_t       fFuncValueMC;
 public:
  BRecParameters(const char* name = 0, const char* title = 0);
  virtual ~BRecParameters();

  Int_t GetNumPar() const { return fNumPar; }
  Int_t GetIsFixed(Int_t index) const { return fIsFixed[index]; }
  Double_t GetInitial(Int_t index) const { return fInitial[index]; }
  Double_t GetMaxLimit(Int_t index) const { return fMaxLimit[index]; }
  Double_t GetMinLimit(Int_t index) const { return fMinLimit[index]; }
  Double_t GetStep(Int_t index) const { return fStep[index]; }
  Double_t GetQual() const { return fQual; }
  Double_t GetFuncValue() const { return fFuncValue; }
  Double_t GetThetaRec() const { return fThetaRec; }
  Double_t GetThetaMC() const { return fThetaMC; }
  Double_t GetPhiMC() const { return fPhiMC; }
  Double_t GetPhiRec() const { return fPhiRec; }
  Double_t GetX0Rec() const { return fX0Rec; }
  Double_t GetY0Rec() const { return fY0Rec; }

  Double_t GetPsi() const;

  void SetNumPar(Int_t value);
  void SetIsFixed(Int_t index, Int_t value) { fIsFixed[index] = value; }
  void SetInitial(Int_t index, Double_t value) { fInitial[index] = value; }
  void SetMinLimit(Int_t index, Double_t value) { fMinLimit[index] = value; }
  void SetMaxLimit(Int_t index, Double_t value) { fMaxLimit[index] = value; }
  void SetStep(Int_t index, Double_t value) { fStep[index] = value; }

  void SetErrors(const Double_t *value);
  void SetFuncValue(Double_t value) { fFuncValue = value; }
  void SetQual(Int_t value) { fQual = value; }

  Int_t GetNhit() const { return fNhit; }
  Double_t GetTres(Int_t i) { return fTres[i]; }
  Double_t GetNchGeom(Int_t i) { return fNchGeom[i]; }
  Double_t GetTresNchGeom(Bool_t ismc, Int_t i) const;

  void SetNhit(Int_t nhit);
  void SetTres(Int_t index, Double_t value) { fTres[index] = value; }
  void SetNchGeom(Int_t index, Int_t value) { fNchGeom[index] = value; }
  void SetTresMC(Int_t nch, Double_t *tresmc, Int_t *nchgeom);

  void SetThetaRec(Double_t th) { fThetaRec = th; }
  void SetThetaMC(Double_t th) { fThetaMC = th; }

  void SetPhiRec(Double_t phi) { fPhiRec = phi; }
  void SetPhiMC(Double_t phi) { fPhiMC = phi; }

  void SetX0Rec(Double_t x) { fX0Rec = x; }
  void SetX0MC(Double_t x) { fX0MC = x; }

  void SetY0Rec(Double_t y) { fY0Rec = y; }
  void SetY0MC(Double_t y) { fY0MC = y; }

  void SetHitP(double p) { fHit_p = p; }
  void SetNoHitP(double np) { fNo_hit_p = np; }

  void     Clear(Option_t* o= "");

  void SetTimeRec(Double_t t) { fTimeRec = t; }
  void SetTimeMC(Double_t t) { fTimeMC = t; }

  void SetFuncValueMC(Double_t f) { fFuncValueMC = f; }

  ClassDef(BRecParameters, 2);   // BRecParameters
};


#endif /*BRECPARAMETERS*/
