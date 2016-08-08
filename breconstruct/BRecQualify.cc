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
}

//______________________________________________________________________
Int_t BRecQualify::PreProcess(MParList *pList)
{
  fEventSource = (BSourceEAS*)pList->FindObject("MCEventSource");
  if (!fEventSource)
    {
      *fLog << err << AddSerialNumber("MCEventSource") << " not found... aborting." << endl;
      return kFALSE;
    }

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
  if (!fMCEventMask)
    {
      *fLog << err << AddSerialNumber("MCEventMask") << " not found... aborting." << endl;
      return kFALSE;
    }
  return kTRUE;
}

//______________________________________________________________________
Int_t BRecQualify::Process()
{
  MCResiduals();
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
  //find first impulse, passed all filters, and take its magic number
  const Double_t sigmat = 2; // temp decision

  int muonnum = -1;
  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    BImpulse *imp = fEvent->GetImpulse(i);
    Int_t nch = imp->GetChannelID();
    BChannelMask::EMode channelflag = fChannelMask->GetFlag(nch);

    //if(fEventMask->GetOrigin(i)->GetFlag() == 1 && channelflag) {
    if(fEventMask->GetOrigin(i)->GetFlag() == 1 && (channelflag == BChannelMask::kOn || channelflag == BChannelMask::kBadChargeCalib || nch == fNchIsRes)) {
      int magic = fMCEventMask->GetOrigin(i)->GetFlag();
      if(fVerbose) {
        cout << " magic = " << magic << endl;
      }
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

  if(muonnum == -1) { // it happens when there are no impulses from muon or cascade
    fRecParameters->SetTresMC(0, 0, 0);
    fRecParameters->SetTimeMC(0);
    fRecParameters->SetThetaMC(-180);
    fRecParameters->SetFuncValueMC(1000);
    fRecParameters->SetX0MC(1000);
    fRecParameters->SetY0MC(1000);
    return;
  }

  BMuon *muon = fEventSource->GetMuonTrack(muonnum - 1);
  if(muon == 0) {
    cout << "BRecQualify::MCResiduals(): very strange situation muon == 0" << endl;
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
    BChannelMask::EMode channelflag = fChannelMask->GetFlag(i);

    //if(nch != nchcuri && fEventMask->GetOrigin(i)->GetFlag() == 1 && channelflag) {
    if(channelflag == BChannelMask::kOn || channelflag == BChannelMask::kBadChargeCalib || i == fNchIsRes) {
      xshift += fGeomTel->At(i)->GetX();
      yshift += fGeomTel->At(i)->GetY();
      zshift += fGeomTel->At(i)->GetZ();
      if(fVerbose) {
        cout << i << " " << fGeomTel->At(i)->GetZ() << endl;
      }
      iii++;
    }
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
