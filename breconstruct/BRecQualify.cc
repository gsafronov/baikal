//////////////////////////////////////////////////////////////////////////////
//
//  BRecQualify
// ===============
//
//
//////////////////////////////////////////////////////////////////////////////
#include "BRecQualify.h"

#include "MLog.h"
#include "MLogManip.h"

#include "MParList.h"
#include "TMath.h"
#include "Math/BrentMinimizer1D.h"
#include "Math/Functor.h"
#include "TTimeStamp.h"

#include <fstream>
#include <iostream>

#include "BSource.h"
#include "BGeomTel.h"
#include "BEvent.h"
#include "BEventMask.h"
#include "BSecInteraction.h"
#include "BRecParameters.h"
#include "BChannelMask.h"
#include <cassert>
#include <cmath>

ClassImp(BRecQualify);

using namespace std;

//______________________________________________________________________
BRecQualify::BRecQualify(const char *name, const char *title) :
  fEventSource(NULL), fEvent(NULL), fEventMask(NULL), fGeomTel(NULL),
  fSecInteraction(NULL), fChannelMask(NULL), fRecParameters(NULL), fMCEventMask(NULL)
{
  fName  = name  ? name  : "BRecQualify";
  fTitle = title ? title : "BRecQualify";

  fFilterEventMaskName = "";

  fVerbose = kFALSE;
  fNchIsRes = -1;
  fQ10FileName = "../config/quant10HE.data";
  fHitCriterionFlag = false;
}

//______________________________________________________________________
Int_t BRecQualify::PreProcess(MParList *pList)
{
  fEventSource = (BSourceEAS*)pList->FindObject("MCEventSource");
  // if (!fEventSource)
  //   {
  //     *fLog << err << AddSerialNumber("MCEventSource") << " not found... aborting." << endl;
  //     return kFALSE;
  //   }

  fEvent = (BEvent*)pList->FindObject("BEvent", "BEvent");
  if (!fEvent)
    {
      *fLog << err << AddSerialNumber("BEvent") << " not found... aborting." << endl;
      return kFALSE;
    }

  fEventMask = (BEventMask*)pList->FindObject(fFilterEventMaskName, "BEventMask");
  if (!fEventMask)
    {
      *fLog << err << AddSerialNumber("BEventMask") << " not found... aborting." << endl;
      return kFALSE;
    }

  fGeomTel = (BGeomTel*)pList->FindObject(AddSerialNumber("BGeomTel"));
  if (!fGeomTel)
    {
      *fLog << err << AddSerialNumber("BGeomTel") << " not found... aborting." << endl;
      return kFALSE;
    }

  fChannelMask = (BChannelMask*)pList->FindObject("BChannelMask");
  if (!fChannelMask)
    {
      *fLog << err << AddSerialNumber("BChannelMask") << " not found... aborting." << endl;
      return kFALSE;
    }

  fSecInteraction = (BMuon*)pList->FindObject(AddSerialNumber("BMuonX0Y0"));
  if(!fSecInteraction) {
    *fLog << err << "Cannot create " << AddSerialNumber("BMuonX0Y0") << endl;
    return kFALSE;
  }

  fRecParameters = (BRecParameters*)pList->FindObject(AddSerialNumber("BRecParameters"));
  if(!fRecParameters) {
    *fLog << err << "Cannot create " << AddSerialNumber("fRecParameters") << endl;
    return kFALSE;
  }

  fMCEventMask = (BEventMask*)pList->FindObject("MCEventMask");
  // if (!fMCEventMask)
  //   {
  //     *fLog << err << AddSerialNumber("MCEventMask") << " not found... aborting." << endl;
  //     return kFALSE;
  //   }

  if (!fHitCriterionFlag) return kTRUE;
  ifstream ifile(fQ10FileName, ios::in);
  if (!ifile) {
    *fLog << err << "file not found: "<< fQ10FileName << endl;
    return kFALSE;
  }
  double rstep , rzero;
  ifile >> rzero;
  ifile >> rstep;
  int kn1;
  ifile >> kn1;
  double* rrrv1 = new double[kn1];
  double* pppv1 = new double[kn1];
  for (int i = 0; i < kn1; i++)
    ifile >> rrrv1[i] >> pppv1[i];
  fInterpolator = new ROOT::Math::Interpolator(kn1, ROOT::Math::Interpolation::kAKIMA);
  fInterpolator->SetData(kn1, rrrv1, pppv1);
  fRhoUpLim = rrrv1[kn1-1];
  delete[] rrrv1;
  delete[] pppv1;

  return kTRUE;
}

//______________________________________________________________________
Int_t BRecQualify::Process()
{
  fXshift = 0;
  fYshift = 0;
  fZshift = 0;


  fTelLength = 0;
  for(int i = 0; i < (int)fGeomTel->GetNumOMs(); i++) {
    for (int j = i+1; j < (int)fGeomTel->GetNumOMs(); j++) {
      double dist = fGeomTel->GetDistance(i,j);
      if (dist > fTelLength) fTelLength = dist;
    }
    fXshift += fGeomTel->At(i)->GetX();
    fYshift += fGeomTel->At(i)->GetY();
    fZshift += fGeomTel->At(i)->GetZ();
  }
  fTelLength += 0.5;
  fXshift = fXshift / fGeomTel->GetNumOMs();
  fYshift = fYshift / fGeomTel->GetNumOMs();
  fZshift = fZshift / fGeomTel->GetNumOMs();

  if (fMCEventMask && fEventSource) {
    fMuonNumber = GetLargestMuonNum();
    MCResiduals();
  }
  if (fHitCriterionFlag) HitCriterion();
  ZDist();
  //MCResidualsX0Y0();
  return kTRUE;
}

//______________________________________________________________________
Int_t BRecQualify::PostProcess()
{
  return kTRUE;
}

//______________________________________________________________________
void BRecQualify::MCResidualsX0Y0()
{
  //find first impulse, passed all filters, and take its magic number
  //BMuon *muon = fEventSource->GetMuonTrack(0);
  //cout << "MCResiduals" << endl;
  const Double_t sigmat = 2; // temp decision

  int muonnum = -1;
  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    BImpulse *imp = fEvent->GetImpulse(i);
    Int_t nch = imp->GetChannelID();
    Int_t channelflag = fChannelMask->GetFlag(nch);

    if(fEventMask->GetOrigin(i)->GetFlag() == 1 && channelflag) {
      int magic = fMCEventMask->GetOrigin(i)->GetFlag();
      //cout << "magic = " << magic << endl;
      if(magic == 0) {
        continue;
      }
      else if(magic < 0) { // muon case
        muonnum = 1000 + magic;
      }
      else {
        int sh = magic/1000;
        muonnum = magic - sh * 1000;
      }
      break;
    }
  }
  //cout << "muonnum = " << muonnum << " numtracks = " << fEventSource->GetNumTracks() << endl;
  BMuon *muon = fEventSource->GetMuonTrack(muonnum - 1);
  if(muon == 0) { // temp dicision because of bug in magic number in MC
    muon = fEventSource->GetMuonTrack(muonnum - 2);
  }

  BMuonX0Y0 muonx0y0;
  muonx0y0.SetX(muon->GetX());
  muonx0y0.SetY(muon->GetY());
  muonx0y0.SetZ(muon->GetZ());
  muonx0y0.SetE(muon->GetE());
  muonx0y0.SetPolarAngle(TMath::Pi() - muon->GetPolarAngle());
  muonx0y0.SetAzimuthAngle(TMath::Pi() + muon->GetAzimuthAngle());

  // shift of coordinate origin along Z
  double zshift = 0;
  int iii = 0;
  int nchcuri = -1;
  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    BImpulse *imp = fEvent->GetImpulse(i);
    Int_t nch = imp->GetChannelID();
    Int_t channelflag = fChannelMask->GetFlag(nch);

    if(nch != nchcuri && fEventMask->GetOrigin(i)->GetFlag() == 1 && channelflag) {
      zshift += fGeomTel->At(nch)->GetZ();
      //cout << nch << " z = " << fGeomTel->At(nch)->GetZ() << endl;
      iii++;
    }
    nchcuri = imp->GetChannelID();
  }
  if(iii != 0) {
    zshift /= iii;
  }
  //cout << "zshift = " << zshift << endl;

  muonx0y0.SetX0(muon->DefineX0(7.8, 16.2, zshift));
  muonx0y0.SetY0(muon->DefineY0(7.8, 16.2, zshift));
  //muonx0y0.SetTime(muon->GetTime() + muon->DefineZ0(7.8, 16.2, zshift) / fGeomTel->GetVelocityVacuum());

  //cout << "muon->GetTime() = " << muon->GetTime() << endl;
  //cout << "muon->DefineZ0(7.8, 16.2, zshift) = " << muon->DefineZ0(7.8, 16.2, zshift) << endl;
  //cout << "x = " << muon->GetX() << " y = " << muon->GetY() << " z = " << muon->GetZ() << " muon->GetPolarAngle() = " << muon->GetPolarAngle() << endl;

  //int kkk;
  //cin >> kkk;

  Int_t      totch = fEvent->GetChannelN();
  Int_t*   nchgeom = new Int_t[totch];
  Double_t* Ttheor = new Double_t[totch];
  Double_t* Texp   = new Double_t[totch];

  // заполнение вектора Ttheor и вычисление T0.
  double T0 = 0;   // начальное время частицы
  int nchcur = -1;
  int ii = 0;
  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    BImpulse *imp = fEvent->GetImpulse(i);
    Int_t nch = imp->GetChannelID();
    Int_t channelflag = fChannelMask->GetFlag(nch);

    if(nch != nchcur && fEventMask->GetOrigin(i)->GetFlag() == 1 && channelflag) {
      nchgeom[ii] = nch;
      Double_t zom = fGeomTel->At(nch)->GetZ() - zshift;

      Ttheor[ii] = muonx0y0.GetStraightLightTime(fGeomTel, 0, 0, zom);
      Texp[ii]   = imp->GetTime();
      T0 += Texp[ii] - Ttheor[ii];
      ii++;
      nchcur = imp->GetChannelID();
    }
  }

  if(ii > 0) {
    T0 /= ii;
  }

  //T0 = muonx0y0.GetTime();

  //cout << endl;
  Double_t chi2 = 0;
  for(int i = 0; i < ii; i++) {
    if(fVerbose) {
      cout << fixed;
      cout.precision(10);
      cout << "Ttheor = " << Ttheor[i] + T0 << " Texp = " << Texp[i] << " dT = " << Ttheor[i] + T0 - Texp[i] << " nchgeom = " << nchgeom[i] << endl;
    }
    Ttheor[i] = Ttheor[i] + T0 - Texp[i];
    chi2 += Ttheor[i] * Ttheor[i];
  }
  chi2 = chi2 / (sigmat * sigmat * (ii - 3));
  //int kkk;
  //cin >> kkk;

  fRecParameters->SetTresMC(ii, Ttheor, nchgeom);

  fRecParameters->SetThetaMC(180 * muon->GetPolarAngle() / TMath::Pi());

  //fRecParameters->SetX0MC(muon->DefineX0(0,0,0));
  //fRecParameters->SetY0MC(muon->DefineY0(0,0,0));



  //fRecParameters->SetX0MC(muon->DefineX0(7.8, 16.2, 244.5));
  //fRecParameters->SetY0MC(muon->DefineY0(7.8, 16.2, 244.5));
  fRecParameters->SetX0MC(muonx0y0.GetX0());
  fRecParameters->SetY0MC(muonx0y0.GetY0());

  // test
  /*Double_t xmuon = muon->GetX();
    Double_t ymuon = muon->GetY();
    Double_t zmuon = muon->GetZ();
    Double_t theta = muon->GetPolarAngle();
    Double_t phi   = muon->GetAzimuthAngle();

    //Double_t x0muon = muon->DefineX0(7.8, 16.2, 244.5);
    //Double_t y0muon = muon->DefineY0(7.8, 16.2, 244.5);
    Double_t x0muon = muon->DefineX0(0, 0, 0);
    Double_t y0muon = muon->DefineY0(0, 0, 0);

    BMuonX0Y0 muonx0y0;

    Double_t par[] = { x0muon, y0muon, zmuon, 0, phi+TMath::Pi(), TMath::Pi() - theta };
    muonx0y0.SetParameters(par);
    muonx0y0.DefineXYZ();

    cout << endl << "xmuon = " << xmuon << " ymuon = " << ymuon << " zmuon = " << zmuon
    << " theta = " << theta << " phi = " << phi << endl;
    cout << "xmuon = " << muonx0y0.GetX() << " ymuon = " << muonx0y0.GetY() << " zmuon = " << muonx0y0.GetZ() <<endl;

    int kkk;
    cin >> kkk;
  */
  fRecParameters->SetFuncValueMC(chi2);

  if(fVerbose) {
    cout << "ThetaMC = " << 180 * muon->GetPolarAngle() / TMath::Pi() << endl;
    cout << "X0MC = " << muon->DefineX0(7.8, 16.2, zshift) << endl;
    cout << "Y0MC = " << muon->DefineY0(7.8, 16.2, zshift) << endl;
    cout << "chi2 = " << chi2 << endl;
  }

  fRecParameters->SetReadyToSave();

  delete [] Ttheor;
  delete [] Texp;
  delete [] nchgeom;

  //cout << "MCResiduals finished" << endl;
}

//______________________________________________________________________
void BRecQualify::MCResiduals()
{
  const Double_t sigmat = 2; // temp decision

  int muonnum = fMuonNumber;

  if(muonnum == -1) { // it happens when there are no impulses from muon or cascade
    fRecParameters->SetTresMC(0, 0, 0);
    fRecParameters->SetTimeMC(0);
    fRecParameters->SetThetaMC(180);
    fRecParameters->SetFuncValueMC(1000);
    fRecParameters->SetX0MC(1000);
    fRecParameters->SetY0MC(1000);
    return;
  }

  BMuon *muon = fEventSource->GetMuonTrack(muonnum - 1); // muonnum begins from 1
  if(muon == 0) {
    cout << "BRecQualify::MCResiduals(): very strange situation muon == 0 muonnum = " << muonnum << endl;
    exit(1);
  }
  //if(muon == 0) { // temp dicision because of bug in magic number in MC
  //muon = fEventSource->GetMuonTrack(muonnum - 2);
  //}

  if(fVerbose) {
    cout << "muon: theta = " << muon->GetPolarAngle() << " phi = " << muon->GetAzimuthAngle() << " x = " << muon->GetX() <<
      " y = " << muon->GetY() << " z = " << muon->GetZ() << " time = " << muon->GetTime() << endl;
  }

  Int_t      totch = fEvent->GetChannelN();
  Int_t*   nchgeom = new Int_t[totch];
  Double_t* Ttheor = new Double_t[totch];
  Double_t* Texp   = new Double_t[totch];

  // заполнение вектора Ttheor и вычисление T0
  double T0 = 0;   // начальное время частицы
  int nchcur = -1;
  int ii = 0;

  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    BImpulse *imp = fEvent->GetImpulse(i);
    Int_t nch = imp->GetChannelID();
    BChannelMask::EMode channelflag = fChannelMask->GetFlag(nch);

    //if(nch != nchcur && fEventMask->GetOrigin(i)->GetFlag() == 1 && channelflag) {
    if(nch != nchcur && fEventMask->GetOrigin(i)->GetFlag() == 1 && (channelflag == BChannelMask::kOn || channelflag == BChannelMask::kBadChargeCalib || nch == fNchIsRes)) {
      nchgeom[ii] = nch;
      Double_t xom = fGeomTel->At(nch)->GetX();
      Double_t yom = fGeomTel->At(nch)->GetY();
      Double_t zom = fGeomTel->At(nch)->GetZ();

      Ttheor[ii] = muon->GetStraightLightTime(fGeomTel, xom, yom, zom);
      Texp[ii]   = imp->GetTime();
      if(fVerbose) {
        cout << "nch = " << nch << " x = " << xom << " y = " << yom << " z = " << zom << " ttheor = " << Ttheor[ii] << " texp = " << Texp[ii] << endl;
      }
      T0 += Texp[ii] - Ttheor[ii];
      ii++;
      nchcur = imp->GetChannelID();
    }
  }

  if(ii > 0) {
    T0 /= ii;
  }
  T0 = 0; //!!!!!! means that it is used time parameter of muon track instead of mean time from response

  Double_t chi2 = 0;
  for(int i = 0; i < ii; i++) {
    if(fVerbose) {
      cout << fixed;
      cout.precision(10);
      cout << "Ttheor = " << Ttheor[i] + T0 << " Texp = " << Texp[i] << " dT = " << Ttheor[i] + T0 - Texp[i] << " nchgeom = " << nchgeom[i] << endl;
    }
    Ttheor[i] = Ttheor[i] + T0 - Texp[i];
    chi2 += Ttheor[i] * Ttheor[i];
  }
  chi2 = chi2 / (sigmat * sigmat * (ii - 3));

  if(fVerbose) {
    cout << "chi2 = " << chi2 << " T0 = " << T0 << endl;
  }


  fRecParameters->SetTresMC(ii, Ttheor, nchgeom);
  fRecParameters->SetTimeMC(muon->GetTime());
  fRecParameters->SetThetaMC(180 * muon->GetPolarAngle() / TMath::Pi());
  fRecParameters->SetPhiMC(180 * muon->GetAzimuthAngle() / TMath::Pi());

  // shift of coordinate origin along Z
  double xshift = 0;
  double yshift = 0;
  double zshift = 0;
  int iii = 0;
  for(int i = 0; i < (int)fGeomTel->GetNumOMs(); i++) {
    //BChannelMask::EMode channelflag = fChannelMask->GetFlag(i);

    //if(nch != nchcuri && fEventMask->GetOrigin(i)->GetFlag() == 1 && channelflag) {
    //if(channelflag == BChannelMask::kOn || channelflag == BChannelMask::kBadChargeCalib || i == fNchIsRes) {
    xshift += fGeomTel->At(i)->GetX();
    yshift += fGeomTel->At(i)->GetY();
    zshift += fGeomTel->At(i)->GetZ();
    if(fVerbose) {
      cout << i << " " << fGeomTel->At(i)->GetZ() << endl;
    }
    iii++;
    //}
  }
  if(iii != 0) {
    xshift /= iii;
    yshift /= iii;
    zshift /= iii;
  }

  //fRecParameters->SetX0MC(muon->DefineX0(7.8, 16.2, 244.5));
  //fRecParameters->SetY0MC(muon->DefineY0(7.8, 16.2, 244.5));
  //fRecParameters->SetX0MC(muon->DefineX0(7.8, 16.2, zshift));
  //fRecParameters->SetY0MC(muon->DefineY0(7.8, 16.2, zshift));
  fRecParameters->SetX0MC(muon->DefineX0(xshift, yshift, zshift));
  fRecParameters->SetY0MC(muon->DefineY0(xshift, yshift, zshift));

  // test
  /*Double_t xmuon = muon->GetX();
    Double_t ymuon = muon->GetY();
    Double_t zmuon = muon->GetZ();
    Double_t theta = muon->GetPolarAngle();
    Double_t phi   = muon->GetAzimuthAngle();

    //Double_t x0muon = muon->DefineX0(7.8, 16.2, 244.5);
    //Double_t y0muon = muon->DefineY0(7.8, 16.2, 244.5);
    Double_t x0muon = muon->DefineX0(0, 0, 0);
    Double_t y0muon = muon->DefineY0(0, 0, 0);

    BMuonX0Y0 muonx0y0;

    Double_t par[] = { x0muon, y0muon, zmuon, 0, phi+TMath::Pi(), TMath::Pi() - theta };
    muonx0y0.SetParameters(par);
    muonx0y0.DefineXYZ();

    cout << endl << "xmuon = " << xmuon << " ymuon = " << ymuon << " zmuon = " << zmuon
    << " theta = " << theta << " phi = " << phi << endl;
    cout << "xmuon = " << muonx0y0.GetX() << " ymuon = " << muonx0y0.GetY() << " zmuon = " << muonx0y0.GetZ() <<endl;

    int kkk;
    cin >> kkk;
  */
  fRecParameters->SetFuncValueMC(chi2);

  if(fVerbose) {
    cout << "zshift = " << zshift << endl << endl;
    cout << "ThetaMC = " << 180 * muon->GetPolarAngle() / TMath::Pi() << endl;
    cout << "X0MC = " << muon->DefineX0(7.8, 16.2, zshift) << endl;
    cout << "Y0MC = " << muon->DefineY0(7.8, 16.2, zshift) << endl;
    cout << "chi2 = " << chi2 << endl;

    cout << "ThetaRec = " << fRecParameters->GetThetaRec() << endl;
    cout << "X0Rec = " << fRecParameters->GetX0Rec() << endl;
    cout << "Y0Rec = " << fRecParameters->GetY0Rec() << endl;
    int kkk;
    cin >> kkk;
  }

  fRecParameters->SetReadyToSave();

  delete [] Ttheor;
  delete [] Texp;
  delete [] nchgeom;
}

//______________________________________________________________________________
double BRecQualify::GetP(double energy) {
  std::vector<int> OMs;
  double hit_p = 1.;
  double no_hit_p = 1.;
  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    int nch = fEvent->GetImpulse(i)->GetChannelID();
    if (!fChannelMask->GetFlag(nch)) continue;
    if (std::find(OMs.begin(), OMs.end(), nch) == OMs.end()) OMs.push_back(nch);
  }
  const double kCAngle = fGeomTel->GetCherenkovAngle();
  //angle between PMT and opposite direction of Cherenkov light
  double angle_below, angle_above;
  if (fTheta < kCAngle) angle_below = kCAngle - fTheta;
  else angle_below = fTheta - kCAngle;
  if (fTheta < M_PI - kCAngle) angle_above = fTheta + kCAngle;
  else angle_above = 2*M_PI - fTheta - kCAngle;
  int num_of_OMs = 0;
  int count = 0;
  for(int i = 0; i < (int)fGeomTel->GetNumOMs(); i++) {
    if (!fChannelMask->GetFlag(i)) continue;
    num_of_OMs++;
    double xOM = fGeomTel->At(i)->GetX() - fXshift;
    double yOM = fGeomTel->At(i)->GetY() - fYshift;
    double zOM = fGeomTel->At(i)->GetZ() - fZshift;
    double x(fX - xOM), y(fY - yOM), z(fZ - zOM);
    // z coordinate of trajectory (or it's projection on a plane parallel to OM
    // axis) and OM axis intersection
    double b;
    if (fTheta == 0 || fTheta == M_PI) b = zOM;
    else b = fZ - x*cos(fPhi)*cos(fTheta)/sin(fTheta)
           - y*sin(fPhi)*cos(fTheta)/sin(fTheta);
    // distance between OM and trajectory
    double rho = sqrt(pow(x*cos(fTheta) - z*cos(fPhi)*sin(fTheta), 2)
                      + pow(sin(fTheta), 2)*pow(y*cos(fPhi) - x*sin(fPhi), 2)
                      + pow(y*cos(fTheta) - z*sin(fTheta)*sin(fPhi), 2));
    double angle;
    if (b < zOM) angle = angle_below;
    else angle = angle_above;
    double pv;
    // cout << scientific << rho << '\t' << energy << '\t' << angle << '\n';
    if (std::find(OMs.begin(), OMs.end(), i) != OMs.end()) {
      double pv = GetProbability(rho, energy, angle, true);
      hit_p *= pv;
      count++;
    }
    else {
      pv = GetProbability(rho, energy, angle, false);
      no_hit_p *= pv;
    }
  }
  double power_hit = 1 / double (count);
  double power_no_hit = 1 / double (num_of_OMs - count);
  assert(!std::isnan(hit_p) && !std::isnan(no_hit_p));
  hit_p = pow(hit_p, power_hit);
  no_hit_p = pow(no_hit_p, power_no_hit);
  return -hit_p*no_hit_p;
}

//______________________________________________________________________________
int BRecQualify::HitCriterion() {
  fTheta = fRecParameters->GetThetaRec();
  fPhi = fRecParameters->GetPhiRec();
  assert(fTheta >= 0 && fTheta <= M_PI);
  assert(fPhi >= 0 && fPhi <= 2*M_PI);
  double X0 = fRecParameters->GetX0Rec();
  double Y0 = fRecParameters->GetY0Rec();
  BMuonX0Y0 muonX0Y0;
  muonX0Y0.SetX0(X0);
  muonX0Y0.SetY0(Y0);
  muonX0Y0.SetPolarAngle(fTheta);
  muonX0Y0.SetAzimuthAngle(fPhi);
  muonX0Y0.DefineXYZ();
  fX = muonX0Y0.GetX() - fXshift;
  fY = muonX0Y0.GetY() - fYshift;
  fZ = muonX0Y0.GetZ() - fZshift;
  ROOT::Math::Functor1D funcP(this, &BRecQualify::GetP);
  ROOT::Math::BrentMinimizer1D bm;
  bm.SetFunction(funcP, 0, 200000);
  bm.Minimize(100, 10, 0);
  fRecParameters->SetP(-bm.FValMinimum());
  fRecParameters->SetE(bm.XMinimum());
  if (fMCEventMask && fEventSource) {
    if(fMuonNumber == -1) {
      fRecParameters->SetTresMC(0, 0, 0);
      fRecParameters->SetTimeMC(0);
      fRecParameters->SetThetaMC(180);
      fRecParameters->SetFuncValueMC(1000);
      fRecParameters->SetX0MC(1000);
      fRecParameters->SetY0MC(1000);
      return kTRUE;
    }
    BMuon *muon = fEventSource->GetMuonTrack(fMuonNumber - 1);
    fTheta = muon->GetPolarAngle();
    fPhi = muon->GetAzimuthAngle();
    // fE = muon->GetE();
    assert(fTheta >= 0 && fTheta <= M_PI);
    assert(fPhi >= 0 && fPhi <= 2*M_PI);
    fX = muon->GetX() - fXshift;
    fY = muon->GetY() - fYshift;
    fZ = muon->GetZ() - fZshift;
    fRecParameters->SetPMC(-GetP(muon->GetE()));
  }
  return kTRUE;
}

//______________________________________________________________________________
double BRecQualify::GetProbability(double rho, double energy, double angle, bool hit) {
  // rho *= 1.e-2; // in main program all distances in cm
  energy *= 1e-3;  // in BSecInteraction energy is in GeV
  if (energy < 0) energy = 0;
  const double Seff_10inch = 506.7074791;  /* 25.4^2*Pi/4  in cm^2*/
  double n_mu = 0.2 + 0.6*energy;
  if (rho > fRhoUpLim) rho = fRhoUpLim;
  double nmean = fInterpolator->Eval(rho);
  /*
    angular dependance of Hamamatsu 10" (lab. measurements 2010)
    x - over Cosine
    this is Jan's approximation
    It contains non zero value at cos=-1
    Lyashuk insists on it (2010)

    c=== BAIKAL 10inchHQ2010 OM sencitivity_coefficients  ===================
    data baicof/0.3082,-0.54192,0.19831,0.04912/
    afectrun=baicof(1)+baicof(2)*ps+baicof(3)*ps**2+baicof(4)*ps**3
  */
  double x = cos(angle);
  double par[]={0.3082,-0.54192,0.19831,0.04912};
  double ham10HQlab = par[0] + x*par[1] + x*x*par[2] + x*x*x*par[3];
  double npe = nmean*n_mu*Seff_10inch*ham10HQlab/rho;
  double prob = exp(-npe);
  // debugging print code
  // if (hit && rho < fRhoUpLim) {
  //   cout.precision(6);
  //   cout << "nmean: " << nmean << "\tn_mu: " << n_mu << "\tham10HQlab: " << ham10HQlab << "\trho: " << rho << "\tX0: " << ((BMuonX0Y0*)fSecInteraction)->GetX0() << "\tY0: " << ((BMuonX0Y0*)fSecInteraction)->GetY0() << "\tprob: " << prob << endl;
  // }
  if (hit) {
    prob = 1 - prob;
    // if (prob > 0.9) cout << scientific << prob << '\t' << rho << '\t' << angle << '\t' << endl;
    // if (prob < 1.e-6) prob = npe;
  }
  // else {
  //   double min = numeric_limits<double>::min();
  //   if (prob < min) prob = min;
  // }
  return prob;
}
//______________________________________________________________________________
double BRecQualify::GetZDist(double theta, double phi) {
  assert(theta >= 0 && theta <= M_PI);
  assert(phi >= 0 && phi <= 2*M_PI);
  double n [] = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
  std::vector<int> channels;
  for (int i=0; i<fEvent->GetTotImpulses(); i++) {
    if (!fEventMask->GetFlag(i)) continue;
    BImpulse *impulse = fEvent->GetImpulse(i);
    int channel = impulse->GetChannelID();
    int channel_flag = fChannelMask->GetFlag(channel);
    if (!channel_flag) continue;
    if (std::find(channels.begin(), channels.end(), channel) == channels.end())
      channels.push_back(channel);
  }
  double z_dist_max = 0;
  for (int i = 0; i < (int) channels.size(); i++) {
    for (int j = i + 1; j < (int) channels.size(); j++) {
      double dx = fGeomTel->At(channels[i])->GetX()
                - fGeomTel->At(channels[j])->GetX();
      double dy = fGeomTel->At(channels[i])->GetY()
                - fGeomTel->At(channels[j])->GetY();
      double dz = fGeomTel->At(channels[i])->GetZ()
                - fGeomTel->At(channels[j])->GetZ();
      double dist = sqrt(dx*dx + dy*dy + dz*dz);
      double dr [] = {dx, dy, dz};
      assert(dist <= fTelLength);
      double asp = 0;
      for (int k = 0; k < 3; k++)
        asp += n[k]*dr[k];
      asp = abs(asp);
      if (asp > z_dist_max) z_dist_max = asp;
    }
  }
  return z_dist_max;
}

//______________________________________________________________________________
void BRecQualify::ZDist() {
  double theta = fRecParameters->GetThetaRec();
  double phi = fRecParameters->GetPhiRec();
  fRecParameters->SetZDist(GetZDist(theta, phi));
  if (!fEventSource || !fMCEventMask) return;
  double thetaMC = fRecParameters->GetThetaMC();
  double phiMC = fRecParameters->GetPhiMC();
  fRecParameters->SetZDistMC(GetZDist(thetaMC, phiMC));
}

//______________________________________________________________________________
Int_t BRecQualify::GetLargestMuonNum() const
{
  // find MC muon number which gave the largest number of filtered impulses
  // if there are no impulses from muons or cascades (pure noise event) return -1

  int muonnum = -1; // output muon number

  // first get number of filtered signal pulses
  int npulses = 0;
  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    BImpulse *imp = fEvent->GetImpulse(i);
    Int_t nch = imp->GetChannelID();
    BChannelMask::EMode channelflag = fChannelMask->GetFlag(nch);

    if(fEventMask->GetOrigin(i)->GetFlag() == 1 && (channelflag == BChannelMask::kOn || channelflag == BChannelMask::kBadChargeCalib)) {
      int magic = fMCEventMask->GetOrigin(i)->GetFlag();
      if(magic != 0) {
        npulses++;
      }
    }
  }

  if(npulses) {
    // second get array for each impulse with muon numbers
    int *muon = new int[npulses];
    int ii = 0;
    for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
      BImpulse *imp = fEvent->GetImpulse(i);
      Int_t nch = imp->GetChannelID();
      BChannelMask::EMode channelflag = fChannelMask->GetFlag(nch);

      if(fEventMask->GetOrigin(i)->GetFlag() == 1 && (channelflag == BChannelMask::kOn || channelflag == BChannelMask::kBadChargeCalib)) {
        int magic = fMCEventMask->GetOrigin(i)->GetFlag();
        if(magic == 0) {
          continue;
        }
        else if(magic < 0) { // muon case
          muonnum = 1000 + magic;
        }
        else {
          int sh = magic/1000;
          muonnum = magic - sh * 1000;
        }
        //cout << i << " nch = " << nch << " magic = " << magic << " muonnum = " << muonnum << endl;
        muon[ii] = muonnum;
        ii++;
      }
    }

    // third sort the array
    int index[npulses];
    TMath::Sort(npulses, muon, index, kFALSE);

    if(fVerbose) {
      for(int i = 0; i < npulses; i++) {
        cout << i << " " << muon[i] << " " << index[i] << endl;
      }

      for(int i = 0; i < npulses; i++) {
        cout << muon[index[i]] << " ";
      }
      cout << endl;
    }

    // fourth the muon number which gave the largest number of impulses.
    // If the number of impulses are equal for some muons take the first one
    int nmax = 1; // number of pulses for maximal muon number
    int nmuonmax = muon[index[0]]; // maximal muon number
    int ncur = 1; // number of pulses for current muon number
    int nmuoncur = muon[index[0]]; // current muon number
    for(int i = 1; i < npulses; i++) {
      if(fVerbose) {
        cout << "nmuoncur = " << nmuoncur << " muon[index[i]] = " << muon[index[i]] << endl;
      }
      if(nmuoncur == muon[index[i]]) {
        ncur++;
      }
      else {
        if(ncur > nmax) {
          nmax = ncur;
          nmuonmax = nmuoncur;
        }
        ncur = 1;
        nmuoncur = muon[index[i]];
      }

      if(fVerbose) {
        cout << "nmax = " << nmax << " nmuonmax = " << nmuonmax << " ncur = " << ncur << " nmuoncur = " << nmuoncur << endl;
        int kkk;
        cin >> kkk;
      }
    }

    muonnum = nmuonmax;
    if(muonnum == 0) {
      for(int i = 0; i < npulses; i++) {
        cout << muon[index[i]] << " ";
      }
      cout << "BRecQualify::GetLargestMuonNum(): strange situation has happened, muon number defined to 0 doesnt exist" << endl;
      exit(1);
    }
  }

  return muonnum;
}
