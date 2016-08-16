//////////////////////////////////////////////////////////////////////////////
//
// BReconstruct
//
// The base class for reconstruction. It gives various instruments needed for event reconstruction.
// Practically, we realize all common methods here that could be used for any type of event
// Specific reconstruction procedure must be realized in derived class (i.e. choice of zero approximation,
// type of minimizer, etc)
//
//////////////////////////////////////////////////////////////////////////////
#include "BReconstruct.h"

#include "MLog.h"
#include "MLogManip.h"

#include "MParList.h"
#include "TMath.h"
#include "TTimeStamp.h"

#include <fstream>
#include <string.h>

#include "BGeomTel.h"
#include "BEvent.h"
#include "BEventMask.h"
#include "BSecInteraction.h"
#include "BRecParameters.h"
#include "BChannelMask.h"

#include <math.h>
#include <vector>
#include "TRotation.h"

// #include "BSource.h"
#include "BMuonCriterionParameters.h"
#include <map>
#include <limits>
#include <assert.h>

ClassImp(BReconstruct);

using namespace std;

//
//
//
BReconstruct::BReconstruct(const char * name, const char * title)
  : fEvent(NULL), fEventMask(NULL), fGeomTel(NULL), fSecInteraction(NULL), fRecParameters(NULL), fChannelMask(NULL), fZeroSecInteraction(NULL), fMuonCriterionParams(NULL) {
  //, fEventSource(NULL), fMCEventMask(NULL) {
  fName  = name  ? name  : "BReconstruct";
  fTitle = title ? title : "Reconstruct";

  fMinimizer = 0;
  fFunctor   = 0;
  fVerbose   = kFALSE;
  //fIsWOAzimuth = kFALSE;
  fFixAzimuth = kFALSE;
  fNchIsRes = -1;
  fInitialsType = 0;

  fEventMaskName = "";
}

//
//
//
Int_t BReconstruct::PreProcess(MParList * pList) {
  if(fEventMaskName == "") {
    * fLog << err << "BReconstruct::PreProcess: Mask name is not defined. Use BReconstruct::SetMaskName()" << endl;
    return kFALSE;
  }

  fEvent = (BEvent*)pList->FindObject("BEvent", "BEvent");
  if (!fEvent)
    {
      * fLog << err << AddSerialNumber("BEvent") << " not found... aborting." << endl;
      return kFALSE;
    }

  fEventMask = (BEventMask*)pList->FindObject(fEventMaskName, "BEventMask");
  if (!fEventMask)
    {
      * fLog << err << AddSerialNumber("BEventMask") << " not found... aborting." << endl;
      return kFALSE;
    }

  fGeomTel = (BGeomTel*)pList->FindObject(AddSerialNumber("BGeomTel"));
  if (!fGeomTel)
    {
      * fLog << err << AddSerialNumber("BGeomTel") << " not found... aborting." << endl;
      return kFALSE;
    }

  fChannelMask = (BChannelMask*)pList->FindObject("BChannelMask");
  if (!fChannelMask)
    {
      * fLog << err << AddSerialNumber("BChannelMask") << " not found... aborting." << endl;
      return kFALSE;
    }

  fZeroSecInteraction = (BMuon*)pList->FindCreateObj("BSecInteraction", "BZeroSecInteraction");
  if(!fZeroSecInteraction) {
    * fLog << err << "Cannot create " << AddSerialNumber("BSecInteraction") << endl;
    return kFALSE;
  }

  // fEventSource = (BSourceEAS*)pList->FindObject("MCEventSource");
  // if (!fEventSource) {
  //   *fLog << err << AddSerialNumber("BEventSource") << " not found... aborting." << endl;
  //   return kFALSE;
  // }

  // fMCEventMask = (BEventMask*)pList->FindObject("MCEventMask");
  // if (!fMCEventMask) {
  // *fLog << err << AddSerialNumber("MCEventMask") << " not found... aborting." << endl;
  // return kFALSE;
  // }

  // Minimizer
  //fMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  fMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "");
  fMinimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  fMinimizer->SetTolerance(0.001);
  if(fVerbose) {
    fMinimizer->SetPrintLevel(2);
  }
  else {
    fMinimizer->SetPrintLevel(-2);
  }


  fFunctor = new ROOT::Math::Functor(this, & BReconstruct::Chi2, 6);
  fMinimizer->SetFunction(* fFunctor);

  return kTRUE;
}


//
//
//
Int_t BReconstruct::Process()
{
  // shift of coordinate origin along Z
  fXshift = 0;
  fYshift = 0;
  fZshift = 0;

  for(int i = 0; i < (int)fGeomTel->GetNumOMs(); i++) {
    fXshift += fGeomTel->At(i)->GetX();
    fYshift += fGeomTel->At(i)->GetY();
    fZshift += fGeomTel->At(i)->GetZ();
  }
  fXshift = fXshift / fGeomTel->GetNumOMs();
  fYshift = fYshift / fGeomTel->GetNumOMs();
  fZshift = fZshift / fGeomTel->GetNumOMs();

  // static int nn = 0;
  // nn++;
  // if(nn%10 == 0) {
  //    cout << nn << endl;
  // }
  Int_t ii = Minimize();
  return ii;
  // return kTRUE;
}

//
//
//
Int_t BReconstruct::PostProcess() {
  return kTRUE;
}

//
// Вычисляет chi2 для любой частицы по времени первого сигнала на ОМ
// T0 вычисляется из требования нулевой производной, а не задается в виде свободного параметра
// Функция имеет строго заданную сигнатуру для того, чтобы ее можно было
// подключать к минимайзеру
//
Double_t BReconstruct::Chi2(const Double_t * x){
  const double sigmat = 2;

  if(fVerbose) {
    cout.precision(10);
    cout << fixed;
    for(int i = 0; i < 6; i++) {
      cout << i << " " << setw(5) << x[i] << endl;
    }
  }

  fSecInteraction->SetParameters(x);

  Int_t      totch = fEvent->GetChannelN();
  Int_t*   nchgeom = new Int_t[totch];
  Double_t* Ttheor = new Double_t[totch];
  Double_t* Texp   = new Double_t[totch];
  //  Double_t* A   = new Double_t[totch];

  // заполнение вектора Ttheor и вычисление T0.
  double T0 = 0;   // начальное время частицы
  int nchcur = -1;
  int ii = 0;
  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    BImpulse * imp = fEvent->GetImpulse(i);
    Int_t nch = imp->GetChannelID();
    BChannelMask::EMode channelflag = fChannelMask->GetFlag(nch);

    if(nch != nchcur && fEventMask->GetOrigin(i)->GetFlag() == 1 && (channelflag == BChannelMask::kOn || channelflag == BChannelMask::kBadChargeCalib)) {
      nchgeom[ii] = nch;
      Double_t xom = fGeomTel->At(nch)->GetX() - fXshift;
      Double_t yom = fGeomTel->At(nch)->GetY() - fYshift;
      Double_t zom = fGeomTel->At(nch)->GetZ() - fZshift;
      //			if(fIsWOAzimuth && fFixAzimuth) {
      //				cout << "Error" << endl;
      //				exit(1);
      //			}
      //			else if(fIsWOAzimuth) {
      //				Ttheor[ii] = fSecInteraction->GetStraightLightTimeWOAzimuth(nchgeom[ii], fGeomTel);
      //			}
      //			else if(fFixAzimuth) {
      //				Ttheor[ii] = fSecInteraction->GetStraightLightTimeFixAzimuth(nchgeom[ii], fGeomTel);
      //			}
      //			else {
      if(fFixAzimuth) {
        Ttheor[ii] = fSecInteraction->GetStraightLightTime(fGeomTel, 0, 0, zom);
      }
      else {
        Ttheor[ii] = fSecInteraction->GetStraightLightTime(fGeomTel, xom, yom, zom);
      }
      //			}

      //	  A[ii]   = imp->GetAmplitude();
      Texp[ii]   = imp->GetTime();
      T0 += Texp[ii] - Ttheor[ii];
      ii++;
      nchcur = imp->GetChannelID();
    }
  }

  if(ii > 0) {
    T0 /= ii;
  }

  double chi2 = 0;
  for(int i = 0; i < ii; i++) {
    chi2 += (Ttheor[i] + T0 - Texp[i]) * (Ttheor[i] + T0 - Texp[i]);
    if(fVerbose) {
      cout << fixed;
      //cout << i << " nchgeom = " << nchgeom[i] << " Ttheor + T0 = " << Ttheor[i] + T0 <<
      //  " Texp = " << Texp[i] << " A = " << A[i] << " dT = " << Ttheor[i] + T0 - Texp[i] << endl;
      cout << i << " nchgeom = " << nchgeom[i] << " Ttheor + T0 = " << Ttheor[i] + T0 <<
        " Texp = " << Texp[i] << " dT = " << Ttheor[i] + T0 - Texp[i] << endl;
    }
  }
  chi2 /= (sigmat * sigmat * (ii - 3));

  if(fVerbose) {
    cout << scientific;
    cout << "chi2 = " << chi2 << endl;
    string tmp;
    getline(cin, tmp);
  }

  if(ii < 4) {
    for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
      BImpulse * imp = fEvent->GetImpulse(i);
      int nchgeom = imp->GetChannelID();
      double Texp   = imp->GetTime();
      cout << fixed;
      cout << i << " nchgeom = " << nchgeom << " Texp = " << Texp
           << " filter = " << fEventMask->GetOrigin(i)->GetFlag() << endl;
    }
    cout << "ii < 4" << endl;
    int kkk;
    cin >> kkk;
  }

  delete [] Ttheor;
  delete [] Texp;
  delete [] nchgeom;
  //delete [] A;
  return chi2;
}

//______________________________________________________________________
Double_t BReconstruct::M_estimator(const Double_t * x) {
  // Функция имеет строго заданную сигнатуру для того, чтобы ее можно было
  // подключать к минимайзеру.
  // Calculates M-estimator function sqrt(1+(t_expected + t_experiment)^2).
  // Its minimization might be more stable than Chi2

  const double sigmat = 2;

  if(fVerbose) {
    cout.precision(10);
    cout << fixed;
    for(int i = 0; i < 6; i++) {
      cout << i << " " << setw(5) << x[i] << endl;
    }
  }

  fSecInteraction->SetParameters(x);

  Int_t      totch = fEvent->GetChannelN();
  Int_t*   nchgeom = new Int_t[totch];
  Double_t* Ttheor = new Double_t[totch];
  Double_t* Texp   = new Double_t[totch];

  // заполнение вектора Ttheor и вычисление T0.
  double T0 = 0;   // начальное время частицы
  int nchcur = -1;
  int ii = 0;
  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    BImpulse * imp = fEvent->GetImpulse(i);
    Int_t nch = imp->GetChannelID();
    BChannelMask::EMode channelflag = fChannelMask->GetFlag(nch);

    if(nch != nchcur && fEventMask->GetOrigin(i)->GetFlag() == 1 && (channelflag == BChannelMask::kOn || channelflag == BChannelMask::kBadChargeCalib)) {
      nchgeom[ii] = nch;
      Double_t xom = fGeomTel->At(nch)->GetX() - fXshift;
      Double_t yom = fGeomTel->At(nch)->GetY() - fYshift;
      Double_t zom = fGeomTel->At(nch)->GetZ() - fZshift;
      //			if(fIsWOAzimuth && fFixAzimuth) {
      //				cout << "Error" << endl;
      //				exit(1);
      //			}
      //			else if(fIsWOAzimuth) {
      //				Ttheor[ii] = fSecInteraction->GetStraightLightTimeWOAzimuth(nchgeom[ii], fGeomTel);
      //			}
      //			else if(fFixAzimuth) {
      //				Ttheor[ii] = fSecInteraction->GetStraightLightTimeFixAzimuth(nchgeom[ii], fGeomTel);
      //			}
      //			else {
      if(fFixAzimuth) {
        Ttheor[ii] = fSecInteraction->GetStraightLightTime(fGeomTel, 0, 0, zom);
      }
      else {
        Ttheor[ii] = fSecInteraction->GetStraightLightTime(fGeomTel, xom, yom, zom);
      }
      //			}

      Texp[ii]   = imp->GetTime();
      T0 += Texp[ii] - Ttheor[ii];
      ii++;
      nchcur = imp->GetChannelID();
    }
  }

  if(ii > 0) {
    T0 /= ii;
  }

  double chi2 = 0;
  for(int i = 0; i < ii; i++) {
    chi2 += (Ttheor[i] + T0 - Texp[i]) * (Ttheor[i] + T0 - Texp[i]);
    if(fVerbose) {
      cout << fixed;
      cout << i << " nchgeom = " << nchgeom[i] << " Ttheor + T0 = " << Ttheor[i] + T0 <<
        " Texp = " << Texp[i] << " dT = " << Ttheor[i] + T0 - Texp[i] << endl;
    }
  }
  chi2 /= (sigmat * sigmat * (ii - 3));

  if(fVerbose) {
    cout << scientific;
    cout << "chi2 = " << chi2 << endl;
    string tmp;
    getline(cin, tmp);
  }

  if(ii < 4) {
    for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
      BImpulse * imp = fEvent->GetImpulse(i);
      int nchgeom = imp->GetChannelID();
      double Texp   = imp->GetTime();
      cout << fixed;
      cout << i << " nchgeom = " << nchgeom << " Texp = " << Texp
           << " filter = " << fEventMask->GetOrigin(i)->GetFlag() << endl;
    }
    cout << "ii < 4" << endl;
    int kkk;
    cin >> kkk;
  }

  delete [] Ttheor;
  delete [] Texp;
  delete [] nchgeom;
  return chi2;
}

//
//
//
Int_t BReconstruct::Minimize() {
  // проверка на минимальное количество сработавших каналов, требуемых для восстановления
  /*  int ii = 0;
      int nchcur = -1;
      for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
      BImpulse * imp = fEvent->GetImpulse(i);
      if(imp->GetChannelID() != nchcur && fEventMask->GetOrigin(i)->GetFlag() == 1) {
      nchcur = imp->GetChannelID();
      ii++;
      }
      }

      if(ii < 5) {
      return kCONTINUE;
      }
  */
  if(fVerbose) {
    cout << "MINMIZER starts. parnum = " << fRecParameters->GetNumPar() << endl;
  }

  fMinimizer->Clear();

  // Set the free variables to be minimized!
  BRecParameters * recpar = fRecParameters;
  char parname[80];
  for(int i = 0; i < recpar->GetNumPar(); i++) {
    sprintf(parname, "x%d", i);
    if(recpar->GetIsFixed(i) == 1) {
      fMinimizer->SetFixedVariable(i, parname, recpar->GetInitial(i));
    }
    else {
      fMinimizer->SetVariable(i, parname, recpar->GetInitial(i), recpar->GetStep(i));
    }
    fMinimizer->SetVariableLimits(i, recpar->GetMinLimit(i), recpar->GetMaxLimit(i));
    if(fVerbose) {
      cout << parname << " isfix = " << recpar->GetIsFixed(i) << " init = " <<
        recpar->GetInitial(i) << " step = " << recpar->GetStep(i) << endl;
    }
  }

  // do the minimization
  fMinimizer->Minimize();

  const double * xs = fMinimizer->X();
  if(fVerbose) {
    cout << fixed;
    std::cout << "Minimum: f(" << xs[0] << ", " << xs[1] << ", " << xs[2] <<
      //xs[3] << ", " << 180* TMath::ASin(TMath::Sin(xs[4])) / 3.1416 << ", " << 180 * TMath::ASin(TMath::Sin(xs[5])) / 3.1416 << "): "
      xs[3] << ", " << xs[4] << ", " << xs[5] << "): "
              << fMinimizer->MinValue() << std::endl;
    cout << "ThetaRec= " << 180 - 180 * xs[5] / 3.1416 << endl;
  }

  fSecInteraction->SetParameters(xs);

  const double * xerr = fMinimizer->Errors();
  fRecParameters->SetErrors(xerr);
  fRecParameters->SetFuncValue(fMinimizer->MinValue());
  fRecParameters->SetQual(fMinimizer->Status());

  fRecParameters->SetThetaRec(180/M_PI*TMath::ACos(TMath::Cos(fSecInteraction->GetPolarAngle())));
  //	fRecParameters->SetThetaRec(180 * TMath::ACos(TMath::Cos(fSecInteraction->GetPolarAngle())) / TMath::Pi());
  double fnpi = fSecInteraction->GetAzimuthAngle() / (2 * TMath::Pi());
  double npi;
  std::modf(fnpi, & npi);
  double phi = fSecInteraction->GetAzimuthAngle() - 2 * TMath::Pi() * npi;
  if(phi < 0) {
    phi = phi + 2 * TMath::Pi();
  }
  //fRecParameters->SetPhiRec(180 * (fSecInteraction->GetAzimuthAngle() - 2 * TMath::Pi() * npi) / TMath::Pi());
  //fRecParameters->SetThetaRec(180 * fSecInteraction->GetPolarAngle() / TMath::Pi());
  fRecParameters->SetPhiRec(phi*180/M_PI);
  fRecParameters->SetX0Rec(((BMuonX0Y0*)fSecInteraction)->GetX0());
  fRecParameters->SetY0Rec(((BMuonX0Y0*)fSecInteraction)->GetY0());

  if (fVerbose) {
    cout << "X0Rec: " << ((BMuonX0Y0*)fSecInteraction)->GetX0() << '\t'
         << "Y0Rec: " << ((BMuonX0Y0*)fSecInteraction)->GetY0() << '\n';
  }
  Residuals();

  fRecParameters->SetReadyToSave();

  if(fVerbose) {
    string tmp;
    getline(cin, tmp);
  }
  return kTRUE;
}


//
//
//
BReconstruct::~BReconstruct() {
  if(fFunctor) {
    delete fFunctor;
  }
  if(fMinimizer) {
    delete fMinimizer;
  }
}

//______________________________________________________________________
void BReconstruct::Print() {

}

//______________________________________________________________________
void BReconstruct::Residuals() {
  const double sigmat = 2;

  //const double xs[6] = {((BMuonX0Y0*)fSecInteraction)->GetX0(), -((BMuonX0Y0*)fSecInteraction)->GetY0(), 0, 0, 0, fSecInteraction->GetPolarAngle()};
  //fSecInteraction->SetParameters(xs);

  Int_t      totch = fEvent->GetChannelN();
  Int_t*   nchgeom = new Int_t[totch];
  Double_t* Ttheor = new Double_t[totch];
  Double_t* Texp   = new Double_t[totch];

  // заполнение вектора Ttheor и вычисление T0.
  double T0 = 0;   // начальное время частицы
  int nchcur = -1;
  int ii = 0;
  std::vector<int> impulse_numbers;
  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    BImpulse * imp = fEvent->GetImpulse(i);
    Int_t nch = imp->GetChannelID();
    BChannelMask::EMode channelflag = fChannelMask->GetFlag(nch);

    //if(nch != nchcur && ((fEventMask->GetOrigin(i)->GetFlag() == 1 && channelflag) || (nch == fNchIsRes))) {
    //cout << nch << " " << fEventMask->GetOrigin(i)->GetFlag() << " " << channelflag << endl;
    if(nch != nchcur && fEventMask->GetOrigin(i)->GetFlag() == 1 && (channelflag == BChannelMask::kOn || channelflag == BChannelMask::kBadChargeCalib || nch == fNchIsRes)) {
      nchgeom[ii] = nch;
      Double_t xom = fGeomTel->At(nch)->GetX() - fXshift;
      Double_t yom = fGeomTel->At(nch)->GetY() - fYshift;
      Double_t zom = fGeomTel->At(nch)->GetZ() - fZshift;

      if(fFixAzimuth) {
        Ttheor[ii] = fSecInteraction->GetStraightLightTime(fGeomTel, 0, 0, zom);
      }
      else {
        Ttheor[ii] = fSecInteraction->GetStraightLightTime(fGeomTel, xom, yom, zom);
      }

      Texp[ii]   = imp->GetTime();
      T0 += Texp[ii] - Ttheor[ii];
      ii++;
      nchcur = imp->GetChannelID();
      impulse_numbers.push_back(i);
    }
  }
  //int kkk;
  //cin >> kkk;
  if(ii > 0) {
    T0 /= ii;
  }

  double chi2 = 0;
  for(int i = 0; i < ii; i++) {
    chi2 += (Ttheor[i] + T0 - Texp[i]) * (Ttheor[i] + T0 - Texp[i]);
  }
  chi2 /= (sigmat * sigmat * (ii - 3));
  //if(chi2 != fRecParameters->GetFuncValue()) {
  //cout << "chi2 = " << chi2 << " fFuncValue = " << fRecParameters->GetFuncValue() << endl;
  //int kkk;
  //cin >> kkk;
  //}

  fRecParameters->SetNhit(ii);
  fRecParameters->SetTimeRec(T0);

  int *pulses = new int[ii];
  for(int i = 0; i < ii; i++) {
    fRecParameters->SetTres(i, Ttheor[i] + T0 - Texp[i]);
    fRecParameters->SetNchGeom(i, nchgeom[i]);
    pulses[i] = impulse_numbers[i];
  }

  fRecParameters->SetImpulseNumbers(pulses);
  delete [] Ttheor;
  delete [] Texp;
  delete [] nchgeom;
  delete [] pulses;
}

//______________________________________________________________________
BSecInteraction* BReconstruct::ZeroRecAnalytic() const {
  // algoritm moved 'as is' from fortran subroutine zero_rec (file: fitting_gvd_z0.f)

  //const Double_t xpow = 1;

  Double_t sa  = 0;
  Double_t st  = 0;
  Double_t sx  = 0;
  Double_t sy  = 0;
  Double_t sz  = 0;
  Double_t stt = 0;
  Double_t sxt = 0;
  Double_t syt = 0;
  Double_t szt = 0;
  Double_t vx = 0;
  Double_t vy = 0;
  Double_t vz = 0;

  int nchcuri = -1;
  int nchcurj = -1;

  for(Int_t i = 0; i < fEvent->GetTotImpulses(); i++) {
    BImpulse * impi = fEvent->GetImpulse(i);
    Int_t nchi = impi->GetChannelID();
    BChannelMask::EMode chflagi = fChannelMask->GetFlag(nchi);
    Int_t impflagi = fEventMask->GetOrigin(i)->GetFlag();

    if(nchi != nchcuri &&  impflagi == 1 && (chflagi == BChannelMask::kOn || chflagi == BChannelMask::kBadChargeCalib)) {
      //Double_t wi   = 1.** xpow;
      Double_t wi   = 1;

      Double_t xi = fGeomTel->At(nchi)->GetX();
      Double_t yi = fGeomTel->At(nchi)->GetY();
      Double_t zi = fGeomTel->At(nchi)->GetZ();
      Double_t ti = impi->GetTime();

      sx += wi * xi;
      sy += wi * yi;
      sz += wi * zi;
      sa += wi;
      st += wi * ti;

      for(Int_t j = 0; j < fEvent->GetTotImpulses(); j++) {
        if (j == i) {
          continue;
        }

        BImpulse * impj = fEvent->GetImpulse(j);
        Int_t nchj = impj->GetChannelID();
        BChannelMask::EMode chflagj = fChannelMask->GetFlag(nchj);
        Int_t impflagj = fEventMask->GetOrigin(j)->GetFlag();

        if(nchj != nchcurj &&  impflagj == 1 && (chflagj == BChannelMask::kOn || chflagj == BChannelMask::kBadChargeCalib)) {
          //Double_t wj   = 1.** xpow;
          Double_t wj   = 1;
          Double_t xj   = fGeomTel->At(nchj)->GetX();
          Double_t yj   = fGeomTel->At(nchj)->GetY();
          Double_t zj   = fGeomTel->At(nchj)->GetZ();
          Double_t tj   = impj->GetTime();
          Double_t dt   = ti - tj;
          //Double_t wtij = wi * wj * dt;

          stt += wi * wj * dt * dt;
          sxt += wi * wj * (xi * (tj * tj - ti * tj) + xj * (ti * ti - ti * tj));
          syt += wi * wj * (yi * (tj * tj - ti * tj) + yj * (ti * ti - ti * tj));
          szt += wi * wj * (zi * (tj * tj - ti * tj) + zj * (ti * ti - ti * tj));
          vx  += wi * wj * (xi - xj) * dt;
          vy  += wi * wj * (yi - yj) * dt;
          vz  += wi * wj * (zi - zj) * dt;

          nchcurj = impj->GetChannelID();
        }
      }

      nchcuri = impi->GetChannelID();
    }
  }

  vx /= stt;
  vy /= stt;
  vz /= stt;

  Double_t v_zero = TMath::Sqrt(vx * vx + vy * vy + vz * vz);

  if(fVerbose) {
    cout << "BReconstruct::ZeroRec():" << endl;
    cout << "v_zero = " << v_zero << endl;
  }

  Double_t x0g = sxt / stt;
  Double_t y0g = syt / stt;
  Double_t z0g = szt / stt;
  Double_t t0  = st  / sa;
  Double_t dcx = vx / v_zero;
  Double_t dcy = vy / v_zero;
  Double_t dcz = vz / v_zero;

  //transform from point (t=0) to (t=t0)

  x0g += vx * t0;
  y0g += vy * t0;
  z0g += vz * t0;

  if(fVerbose) {
    cout << "x0g = " << x0g << endl;
    cout << "y0g = " << y0g << endl;
    cout << "z0g = " << z0g << endl;
  }

  Double_t t_zero = t0;

  Double_t theta, phi;
  Tauinv(dcx, dcy, dcz, theta, phi);

  Double_t theta_zero = 180. - theta;
  Double_t phi_zero = phi - 180.;

  Double_t t1, t2, f1, f2;
  To_cos(theta_zero, phi_zero, t1, t2, f1, f2);

  // for the comparision with multiple zero approximation

  Double_t t11 = t1 * f1;
  Double_t t21 = t1 * f2;
  Double_t t31 = -t2;
  Double_t t12 = -f2;
  Double_t t22 = f1;
  Double_t t32 = 0;
  //Double_t t13 = t2 * f1;
  //Double_t t23 = t2 * f2;
  //Double_t t33 = t1;

  Double_t x_zero = x0g * t11 + y0g * t21 + z0g * t31;
  Double_t y_zero = x0g * t12 + y0g * t22 + z0g * t32;

  if(fVerbose) {
    cout << "guess:    x,y,z = " << x_zero << " " << y_zero;
    cout << "      the,phi,t = " << theta_zero << " " << phi_zero << " " << t_zero << endl;
  }

  BSecInteraction * sec = new BSecInteraction();
  sec->SetX(x_zero);
  sec->SetY(y_zero);
  sec->SetPolarAngle(theta_zero * TMath::Pi()/180);
  sec->SetAzimuthAngle(phi_zero * TMath::Pi()/180);
  sec->SetTime(t_zero);

  return sec;
}


//____________________________________________________________________________________________________
void BReconstruct::Tauinv(Double_t txf, Double_t tyf, Double_t tzf, Double_t & theta, Double_t & phi) const {
  //inverse function of tau

  const Double_t grad = TMath::RadToDeg();

  theta = TMath::ACos(tzf);

  if(txf == 0) {
    phi = TMath::PiOver2();
  }
  else {
    phi = TMath::ATan2(tyf, txf);
  }

  if (tyf < 0) {
    phi = phi + 2 * TMath::Pi();
  }

  phi *= grad;
  theta *= grad;

  return;
}

//___________________________________________________________________________________________________________
void BReconstruct::To_cos(Double_t the, Double_t phi, Double_t & t1, Double_t & t2, Double_t & f1, Double_t & f2) const {
  const Double_t grad = TMath::RadToDeg();

  Double_t thet = the / grad;
  Double_t phit = phi / grad;

  t1 = TMath::Cos(thet);
  t2 = TMath::Sin(thet);
  f1 = TMath::Cos(phit);
  f2 = TMath::Sin(phit);

  return;
}
