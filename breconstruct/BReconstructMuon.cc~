//////////////////////////////////////////////////////////////////////////////
//
//  BReconstructMuon
//
//
//////////////////////////////////////////////////////////////////////////////
#include "BReconstructMuon.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "MLog.h"
#include "MLogManip.h"

#include "MParList.h"
#include "TMath.h"
#include "TTimeStamp.h"

#include <fstream>
#include <iostream>

#include "BGeomTel.h"
#include "BEvent.h"
#include "BEventMask.h"
#include "BSecInteraction.h"
#include "BRecParameters.h"
#include "BChannelMask.h"

#include "TH2F.h"
#include "TFile.h"

ClassImp(BReconstructMuon);

using namespace std;

//
//
//
BReconstructMuon::BReconstructMuon(const char *name, const char *title)
{
  fName  = name  ? name  : "BReconstructMuon";
  fTitle = title ? title : "ReconstructMuon";
  fQ10Filename = "../config/quant10HE.data";
  fHitProbabilityCriterionFlag = false;

  fou=new TFile("recoMonitor.root","RECREATE");
  hMapOfUsedOM=new TH2F("hMapOfUsedOM","hMapOfUsedOM",10,0,10,25,0,25);
}

//______________________________________________________________________
//
//
//
Int_t BReconstructMuon::PostProcess()
{
  fou->cd();
  hMapOfUsedOM->Write();
  fou->Close();    
  return 1;
}


Int_t BReconstructMuon::PreProcess(MParList *pList)
{
  if(!BReconstruct::PreProcess(pList)) {
    return kFALSE;
  }

  fSecInteraction = (BMuon*)pList->FindCreateObj(AddSerialNumber("BMuonX0Y0"));
  if(!fSecInteraction) {
    *fLog << err << "Cannot create " << AddSerialNumber("BMuonX0Y0") << endl;
    return kFALSE;
  }

  fRecParameters = (BRecParameters*)pList->FindCreateObj(AddSerialNumber("BRecParameters"));
  if(!fRecParameters) {
    *fLog << err << "Cannot create " << AddSerialNumber("fRecParameters") << endl;
    return kFALSE;
  }

  // define initial reconstruction parameters
  const Int_t numpar = 6;

  fRecParameters->SetNumPar(numpar);
  /*    for(int i = 0; i < numpar; i++) {
        if(i == 4 && fIsWOAzimuth == kTRUE) {
        fRecParameters->SetIsFixed(i, 1);
        }
        if(i == 4 && fFixAzimuth == kTRUE) {
        fRecParameters->SetIsFixed(i, 1);
        }

        if(i == 2 || i == 3) {
        fRecParameters->SetIsFixed(i, 1);
        }
        else {
        fRecParameters->SetIsFixed(i, 0);
        }
        fRecParameters->SetInitial(i, 0.5);
        fRecParameters->SetMinLimit(i, 0);
        fRecParameters->SetMaxLimit(i, 0);
        fRecParameters->SetStep(i, 1.1);
        }
  */

  fRecParameters->SetIsFixed(0, 0);
  fRecParameters->SetInitial(0, 0.5);
  fRecParameters->SetMinLimit(0, 0);
  fRecParameters->SetMaxLimit(0, 0);
  fRecParameters->SetStep(0, 10.0);

  fRecParameters->SetIsFixed(1, 0);
  fRecParameters->SetInitial(1, 20);
  fRecParameters->SetMinLimit(1, 0);
  fRecParameters->SetMaxLimit(1, 0);
  fRecParameters->SetStep(1, 10.0);

  fRecParameters->SetIsFixed(2, 1);
  fRecParameters->SetInitial(2, 0.5);
  fRecParameters->SetMinLimit(2, 0);
  fRecParameters->SetMaxLimit(2, 0);
  fRecParameters->SetStep(2, 1.0);

  fRecParameters->SetIsFixed(3, 1);
  fRecParameters->SetInitial(3, 0.5);
  fRecParameters->SetMinLimit(3, 0);
  fRecParameters->SetMaxLimit(3, 0);
  fRecParameters->SetStep(3, 1.0);

  //    if(fIsWOAzimuth == kTRUE || fFixAzimuth == kTRUE) {
  if(fFixAzimuth == kTRUE) {
    fRecParameters->SetIsFixed(4, 1);
    fRecParameters->SetInitial(4, 0);
  }
  else {
    fRecParameters->SetIsFixed(4, 0);
    fRecParameters->SetInitial(4, 134*3.1416/180);
  }
  fRecParameters->SetMinLimit(4, -100 * 3.1416);
  fRecParameters->SetMaxLimit(4,  100 * 3.1416);
  fRecParameters->SetStep(4, 0.2);

  fRecParameters->SetIsFixed(5, 0);
  fRecParameters->SetInitial(5, 116*3.1416/180);
  fRecParameters->SetMinLimit(5, -100 * 3.1416);
  fRecParameters->SetMaxLimit(5,  100 * 3.1416);
  fRecParameters->SetStep(5, 0.1);

  if (!fHitProbabilityCriterionFlag) return kTRUE;

  ifstream ifile(fQ10Filename, ios::in);
  if (!ifile) {
    *fLog << err << "file not found: "<< fQ10Filename << endl;
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
  delete[] rrrv1;
  delete[] pppv1;
  fRho_up_lim = rrrv1[kn1-1];

  return kTRUE;
}

//______________________________________________________________________
Int_t BReconstructMuon::Process()
{
  std::cout<<"fff"<<std::endl;
  BReconstruct::Process();

  if (!fHitProbabilityCriterionFlag) return kTRUE;

  ((BMuonX0Y0*)fSecInteraction)->DefineXYZ();
  fSecInteraction->SetReadyToSave();
  HitProbabilityCriterion();

  return kTRUE;
}
//______________________________________________________________________________
int BReconstructMuon::HitProbabilityCriterion() {
  const double kCAngle = fGeomTel->GetCherenkovAngle();

  std::vector<int> OMs;
  std::cout<<"hit probability..."<<std::endl;
  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    int nch = fEvent->GetImpulse(i)->GetChannelID();
    BChannelMask::EMode channelflag = fChannelMask->GetFlag(nch);
    std::cout<<"yoyo"<<std::endl;
    if(fEventMask->GetOrigin(i)->GetFlag() == 1 && (channelflag == BChannelMask::kOn || channelflag == BChannelMask::kBadChargeCalib)) {
      if (std::find(OMs.begin(), OMs.end(), nch) == OMs.end()) {
	OMs.push_back(nch);
	std::cout<<"bloblo"<<std::endl;
	hMapOfUsedOM->Fill(floor(nch/24),(nch)%24+1,1);
      }
    }
  }

  double hit_p(1), no_hit_p(1);
  double theta(M_PI*fRecParameters->GetThetaRec() / 180), phi(M_PI*fRecParameters->GetPhiRec() / 180);
  double angle_below, angle_above;
  if (theta < kCAngle) angle_below = kCAngle - theta;
  else angle_below = theta - kCAngle;
  if (theta < M_PI - kCAngle) angle_above = theta + kCAngle;
  else angle_above = 2*M_PI - theta - kCAngle;
  for(int i = 0; i < (int)fGeomTel->GetNumOMs(); i++) {
    double xOM = fGeomTel->At(i)->GetX() - fXshift;
    double yOM = fGeomTel->At(i)->GetY() - fYshift;
    double zOM = fGeomTel->At(i)->GetZ() - fZshift;
    double X(fSecInteraction->GetX()), Y(fSecInteraction->GetY()), Z(fSecInteraction->GetZ());
    double x(X - xOM), y(Y - yOM), z(Z - zOM);
    double b; // z coordinate of trajectory (or it's projection on a plane parallel to OM axis) and OM axis intersection
    if (theta == 0 || theta == M_PI) b = zOM;
    else b = Z - x*cos(phi)*cos(theta)/sin(theta) - y*sin(phi)*cos(theta)/sin(theta);
    double rho = sqrt(pow(x*cos(theta) - z*cos(phi)*sin(theta), 2) + pow(sin(theta), 2)*pow(y*cos(phi) - x*sin(phi), 2) + pow(y*cos(theta) - z*sin(theta)*sin(phi), 2)); // distance between OM and trajectory
    double angle;
    if (b < zOM) angle = angle_below;
    else angle = angle_above;
    if (std::find(OMs.begin(), OMs.end(), i) != OMs.end()) hit_p *= GetProbability(rho, angle, true);
    else no_hit_p *= GetProbability(rho, angle, false);
  }
  fRecParameters->SetHitP(hit_p);
  fRecParameters->SetNoHitP(no_hit_p);
  return kTRUE;
}

//___________________________________________________________________________________________________________
double BReconstructMuon::GetProbability(double rho, double angle, bool hit) {
  // rho *= 1.e-2; // in main program all distances in cm
  const double Seff_10inch = 506.7074791;  /* 25.4^2*Pi/4  in cm^2*/
  const double log10e = log10(0.4); // 0.4 TeV is the approximate mean energy of the muon
  double n_mu = 0.2 + 0.6*pow(10,log10e);
  if (rho > fRho_up_lim) rho = fRho_up_lim;
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
  double x = - cos(angle);
  double par[]={0.3082,-0.54192,0.19831,0.04912};
  double ham10HQlab = par[0] + x*par[1] + x*x*par[2] + x*x*x*par[3];
  double npe = nmean*n_mu*Seff_10inch*ham10HQlab/rho;
  double prob = exp(-npe);
  // debugging print code
  // if (hit && rho < fRho_up_lim) {
  //   cout.precision(6);
  //   cout << "nmean: " << nmean << "\tn_mu: " << n_mu << "\tham10HQlab: " << ham10HQlab << "\trho: " << rho << "\tX0: " << ((BMuonX0Y0*)fSecInteraction)->GetX0() << "\tY0: " << ((BMuonX0Y0*)fSecInteraction)->GetY0() << "\tprob: " << prob << endl;
  // }
  if (hit) {
    prob = 1 - prob;
    if (prob < 1.e-6) prob = npe;
  }
  else {
    double min = numeric_limits<double>::min();
    if (prob < min) prob = min;
  }
  return prob;
}
