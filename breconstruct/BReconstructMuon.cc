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
#include "BMuonCriterionParameters.h"
#include <assert.h>
// #include "BSource.h"

ClassImp(BReconstructMuon);

using namespace std;

//
//
//
BReconstructMuon::BReconstructMuon(const char *name, const char *title)
{
  fName  = name  ? name  : "BReconstructMuon";
  fTitle = title ? title : "ReconstructMuon";
}

//______________________________________________________________________
//
//
//
Int_t BReconstructMuon::PreProcess(MParList *pList)
{
  if(!BReconstruct::PreProcess(pList)) {
    return kFALSE;
  }

  fZeroSecInteraction = (BMuonX0Y0*)pList->FindCreateObj("BMuonX0Y0", "BZeroMuonX0Y0");
  if(!fZeroSecInteraction) {
    * fLog << err << "Cannot create " << AddSerialNumber("BMuonX0Y0") << endl;
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

  if (fInitialsType == 2) {
    fMuonCriterionParams = (BMuonCriterionParameters*)pList->FindObject(AddSerialNumber("BMuonCriterionParameters"));
    if(!fMuonCriterionParams) {
      *fLog << err << "Cannot find " << AddSerialNumber("BMuonCriterionParameters") << endl;
      return kFALSE;
    }
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

  return kTRUE;
}

//______________________________________________________________________
Int_t BReconstructMuon::Process()
{
  if (InitialValues()) return kCONTINUE;
  return BReconstruct::Process();
}
//______________________________________________________________________________
int BReconstructMuon::InitialValues() {
  // Double_t xmin  = -10;
  // Double_t xstep =   2;
  // Int_t xn =   10;

  // Double_t ymin  = -10;
  // Double_t ystep =   2;
  // Int_t yn =   10;

  //Double_t phimin  = 0;
  //Double_t phistep = 10 * TMath::Pi() / 180;
  //Int_t phin =   36;

  // Double_t thetamin  = 90 * TMath::Pi() / 180;
  // Double_t thetastep = 10 * TMath::Pi() / 180;
  // Int_t thetan =   9;

  // double chimin = 10000000;
  double xinit = 0;
  double yinit = 0;
  double phiinit = 0;
  double thetainit = 0;

  static int counter;

  if (fInitialsType == 2) {
    thetainit = (fMuonCriterionParams->GetMinA() + fMuonCriterionParams->GetMaxA())/2;
    GetX0Y0(xinit, yinit, thetainit, phiinit);
  }
  if (fInitialsType == 0) {
    // plane wave algorithm launches in any case
    BSecInteraction * sec = PlaneWaveSingleMuonDirectionReconstruction();

    if(!sec) {
      if (!fMuonCriterionParams) return 1;
      // degenerated case. So use other zero approximation for zenith. There is no zero approximation for azimuth at the moment
      *fLog << inf << "sec == 0 in " << counter++ << "th time\n";
      thetainit = (fMuonCriterionParams->GetMaxA() + fMuonCriterionParams->GetMinA())/2;
      phiinit = 0;
    }
    else {
      // main case of zero approximation
      fZeroSecInteraction->SetPolarAngle(sec->GetPolarAngle());
      fZeroSecInteraction->SetAzimuthAngle(sec->GetAzimuthAngle());
      fZeroSecInteraction->SetReadyToSave();
      delete sec;
      phiinit = fZeroSecInteraction->GetAzimuthAngle();
      thetainit = fZeroSecInteraction->GetPolarAngle()+0.54;
      if(fVerbose) {
        cout << "sec->GetPolarAngle() = " << sec->GetPolarAngle() << " sec->GetAzimuthAngle() = " << sec->GetAzimuthAngle() << endl;
      }
    }

    GetX0Y0(xinit, yinit, thetainit, phiinit); // zero approximation for X0 and Y0
    ((BMuonX0Y0*)fZeroSecInteraction)->SetX0(xinit);
    ((BMuonX0Y0*)fZeroSecInteraction)->SetY0(yinit);
    //fRecParameters->SetMinLimit(5, fMuonCriterionParams->GetMinA());
    //fRecParameters->SetMaxLimit(5, fMuonCriterionParams->GetMaxA());
  }
  // for(int ix = 0; ix <= xn; ix++) {
  //   double x = xmin + ix * xstep;
  //   for(int iy = 0; iy <= yn; iy++) {
  //     double y = ymin + iy * ystep;
  //     //for(int iphi = 0; iphi < phin; iphi++) {
  //     //double phi = phimin + iphi * phistep;
  //     double phi = 0;
  //     for(int itheta = 0; itheta < thetan; itheta++) {
  //       double theta = thetamin + itheta * thetastep;
  //       double r[] = { x, y, 0, 0, phi, theta };
  //       double chi = Chi2(r);

  //       if(chi < chimin) {
  //         chimin = chi;
  //         xinit = x;
  //         yinit = y;
  //         phiinit = phi;
  //         thetainit = theta;
  //       }
  //     }
  //     //}
  //   }
  // }

  if(fVerbose) {
	  cout << "xinit = " << xinit << " yinit = " << yinit << " phiinit = " << phiinit << " thetainit = " << thetainit << endl;
  }

  fRecParameters->SetInitial(0, xinit);
  fRecParameters->SetInitial(1, yinit);
  fRecParameters->SetInitial(4, phiinit);
  fRecParameters->SetInitial(5, thetainit);
  return 0;
}

//______________________________________________________________________________
////////////////////////////////////////////////////////////////////////////////
// The idea of plane wave model direction reconstruction algorithm is to define
// components of plane's normal vector from a system of equations:
// c*dt1 == nx*dx1 + ny*dy1 + nz*dz1 &&
// c*dt2 == nx*dx2 + ny*dy2 + nz*dz2 &&
// nx*nx + ny*ny + nz*nz == 1,
// where dti is the time difference for i-th pair of channels and the same
// notations for coordinates are used. This system has two sets of solutions:
// {{nx1, ny1, nz1}, {nx2, ny2, nz2}}.
// The algorithm takes all sets of pairs of channel pairs and for each pair
// finds a solution. From all the solutions it picks the one that gives the
// minimum absolute value for difference vec(n)*vec(dr3) - c*dt3, where
// parameters of a third pair, which gives the minimum of the difference for a
// particular set of two pairs, are used.
//
// If solution was not found, alternative method can be applied either after
// previous one or instead of it (skipping Combinations method invocation).
// The method simply takes avarage of velocity vectors v_ij = dr_ij / dt_ij.
// Then normal vector of the direction is calculated.
////////////////////////////////////////////////////////////////////////////////

BSecInteraction * BReconstructMuon::PlaneWaveSingleMuonDirectionReconstruction() {
  int previous_channel_i = -1;
  int previous_channel_j = -1;
  vector<pair> set_of_pairs;
  for (int i = 0; i < fEvent->GetTotImpulses(); i++) {
    int impulse_flag_i = fEventMask->GetOrigin(i)->GetFlag();
    if (!impulse_flag_i) continue;
    BImpulse * impulse_i = fEvent->GetImpulse(i);
    int channel_i = impulse_i->GetChannelID();
    if (channel_i == previous_channel_i) continue;
    BChannelMask::EMode channel_flag_i = fChannelMask->GetFlag(channel_i);
    if (channel_flag_i != BChannelMask::kOn
        && channel_flag_i != BChannelMask::kBadChargeCalib) continue;
    pair a_pair;
    for (int j = i + 1; j < fEvent->GetTotImpulses(); j++) {
      int impulse_flag_j = fEventMask->GetOrigin(j)->GetFlag();
      if (impulse_flag_j != 1) continue;
      BImpulse * impulse_j = fEvent->GetImpulse(j);
      int channel_j = impulse_j->GetChannelID();
      if (channel_j == previous_channel_j) continue;
      BChannelMask::EMode channel_flag_j = fChannelMask->GetFlag(channel_j);
      if (channel_flag_j != BChannelMask::kOn
          && channel_flag_j != BChannelMask::kBadChargeCalib) continue;
      if (channel_j == channel_i) continue;
      a_pair.first = i;
      a_pair.second = j;
      set_of_pairs.push_back(a_pair);
      previous_channel_j = channel_j;
      previous_channel_i = channel_i;
    }
  }

  int need = 2;
  if ((int)set_of_pairs.size() < need) return NULL;
  const double max_double = std::numeric_limits<double>::max();
  double parameters [] = {max_double, 0, 0, 0};
  // Combinations(set_of_pairs.size(), need, 0, 0, set_of_pairs, parameters);

  // debug print
  // if (parameters[0] == max_double) {
  //   cout << "\nNew Event\n";
  //   for (int i = 0; i < fEvent->GetTotImpulses(); i++) {
  //     if( fEventMask->GetOrigin(i)->GetFlag() != 1) continue;
  //     BImpulse * impulse = fEvent->GetImpulse(i);
  //     int ch_id = impulse->GetChannelID();
  //     BGeom * channel = fGeomTel->At(ch_id);
  //     double x = channel->GetX();
  //     double y = channel->GetY();
  //     double z = channel->GetZ();
  //     double t = impulse->GetTime();
  //     cout << i << ":\t" << ch_id << '\t' << x << '\t' << y << '\t' << z << '\t'
  //          << t << endl;
  //   }
  //   int count;
  //   for (auto &p: set_of_pairs)
  //     cout << "pair " << ++count << '\t' << p.first << '\t' << p.second << endl;
  //   parameters[0] = e10d;
  //   Combinations(set_of_pairs.size(), need, 0, 0, set_of_pairs, parameters);
  //   cin >> count;
  // }

  if (parameters[0] == max_double) {
    double velocity [3];
    double zero = 1e-15;
    for (pair &a_pair: set_of_pairs) {
      BImpulse * first_impulse = fEvent->GetImpulse(a_pair.first);
      BImpulse * second_impulse = fEvent->GetImpulse(a_pair.second);
      double t1 = first_impulse->GetTime();
      double t2 = second_impulse->GetTime();
      double dt = t1 - t2;
      if (abs(dt) < zero) continue;
      BGeom * channel_firts = fGeomTel->At(first_impulse->GetChannelID());
      BGeom * channel_second = fGeomTel->At(second_impulse->GetChannelID());
      double dx = channel_firts->GetX() - channel_second->GetX();
      double dy = channel_firts->GetY() - channel_second->GetY();
      double dz = channel_firts->GetZ() - channel_second->GetZ();
      velocity[0] += dx/dt;
      velocity[1] += dy/dt;
      velocity[2] += dz/dt;
    }
    double velocity_squered(1);
    for (int i=0; i<3; i++)
      velocity_squered += pow(velocity[i], 2);
    for (int i=0; i<3; i++) {
      velocity[i] /= sqrt(velocity_squered);
      parameters[i+1] = velocity[i];
    }
  }
  assert(parameters[1] <=1 && parameters[2] <= 1 && parameters[3] <= 1);
  BSecInteraction * sec = new BSecInteraction();
  double theta = acos(parameters[3]);
  assert(theta >= 0 && theta <= M_PI);
  sec->SetPolarAngle(theta);
  double phi = atan(parameters[2]/parameters[1]);
  if (parameters[1] > 0) {
    if (parameters[2] < 0) phi += 2*M_PI;
  }
  else phi += M_PI;
  assert(phi <= 2*M_PI && phi >= 0);
  sec->SetAzimuthAngle(phi);
  return sec;
}

//______________________________________________________________________________
void BReconstructMuon::Combinations(int pool, int need, unsigned long chosen,
                                int at, vector<pair> & set_of_pairs,
                                double * parameters) {
  if (pool < need + at) return;
  unsigned long one = 1;
  double zero = 1e-15;
  if (!need) {
    vector<int> chosen_pairs;
    vector<double> dt;
    vector<double> dx;
    vector<double> dy;
    vector<double> dz;
    for (at = 0; at < pool; at++)
      if (chosen & (one << at)) {
        BImpulse * first_impulse = fEvent->GetImpulse(set_of_pairs[at].first);
        BImpulse * second_impulse = fEvent->GetImpulse(set_of_pairs[at].second);
        double t1 = first_impulse->GetTime();
        double t2 = second_impulse->GetTime();
        dt.push_back((t1-t2) * fGeomTel->GetVelocityVacuum());
        BGeom * channel_firts = fGeomTel->At(first_impulse->GetChannelID());
        BGeom * channel_second = fGeomTel->At(second_impulse->GetChannelID());
        dx.push_back(channel_firts->GetX() - channel_second->GetX());
        dy.push_back(channel_firts->GetY() - channel_second->GetY());
        dz.push_back(channel_firts->GetZ() - channel_second->GetZ());
        chosen_pairs.push_back(at);
      }
    double dxy = dx[0]*dy[1] - dx[1]*dy[0];
    double dxz = dx[0]*dz[1] - dx[1]*dz[0];
    double dzy = dz[0]*dy[1] - dz[1]*dy[0];
    double dty = dt[0]*dy[1] - dt[1]*dy[0];
    double dxt = dx[0]*dt[1] - dx[1]*dt[0];
    double dzt = dz[0]*dt[1] - dz[1]*dt[0];
    double Dxy = 2*dty*dxt*dxz*dzy - dty*dty*(dxy*dxy + dxz*dxz)
      - dxt*dxt*(dxy*dxy + dzy*dzy) + dxy*dxy*(dxy*dxy + dxz*dxz + dzy*dzy);
    double Dzy = - 2*dty*dzt*dxz*dxy - dty*dty*(dzy*dzy + dxz*dxz)
      - dzt*dzt*(dzy*dzy + dxy*dxy) + dzy*dzy*(dzy*dzy + dxz*dxz + dxy*dxy);
    double Dxz = 2*dzt*dxt*dxy*dzy - dzt*dzt*(dxz*dxz + dxy*dxy)
      - dxt*dxt*(dxz*dxz + dzy*dzy) + dxz*dxz*(dxz*dxz + dxy*dxy + dzy*dzy);
    const int two = 2;
    double nx[two], ny[two], nz[two];
    double denominator = dxy*dxy + dxz*dxz + dzy*dzy;
    if ((abs(dxy) > zero) && (Dxy >= 0)) {
      double a = dxt*dxz + dty*dzy;
      nz[0] = (a - sqrt(Dxy)) / denominator;
      nz[1] = (a + sqrt(Dxy)) / denominator;
      for (int i=0; i<two; i++) {
        nx[i] = (dty - dzy*nz[i]) / dxy;
        ny[i] = (dxt - dxz*nz[i]) / dxy;
      }
    }
    else if ((abs(dxz) > zero) && (Dxz >= 0)) {
      double a = dxt*dxy + dzt*dzy;
      ny[0] = (a - sqrt(Dxz)) / denominator;
      ny[1] = (a + sqrt(Dxz)) / denominator;
      for (int i=0; i<two; i++) {
        nx[i] = (- dzt + dzy*ny[i]) / dxz;
        nz[i] = (dxt - dxy*ny[i]) / dxz;
      }
    }
    else if ((abs(dzy) > zero) && (Dzy >= 0)) {
      double a = dty*dxy - dzt*dxz;
      nx[0] = (a - sqrt(Dzy)) / denominator;
      nx[1] = (a + sqrt(Dzy)) / denominator;
      for (int i=0; i<two; i++) {
        ny[i] = (dzt + dxz*nx[i]) / dzy;
        nz[i] = (dty - dxy*nx[i]) / dzy;
      }
    }
    else return;
    for (int i=0; i<(int)set_of_pairs.size(); i++) {
      if (i == chosen_pairs[0] || i == chosen_pairs[1]) continue;
      BImpulse * first_impulse = fEvent->GetImpulse(set_of_pairs[i].first);
      BImpulse * second_impulse = fEvent->GetImpulse(set_of_pairs[i].second);
      double t1 = first_impulse->GetTime();
      double t2 = second_impulse->GetTime();
      double dt_third = (t1-t2) * fGeomTel->GetVelocityVacuum();
      BGeom * channel_firts = fGeomTel->At(first_impulse->GetChannelID());
      BGeom * channel_second = fGeomTel->At(second_impulse->GetChannelID());
      double dx_third = channel_firts->GetX() - channel_second->GetX();
      double dy_third = channel_firts->GetY() - channel_second->GetY();
      double dz_third = channel_firts->GetZ() - channel_second->GetZ();
      double res [two];
      for (int j=0; j<two; j++) {
        res[j] = abs(nx[j]*dx_third + ny[j]*dy_third + nz[j]*dz_third - dt_third);
        if (res[j] < parameters[0]) {
          parameters[0] = res[j];
          parameters[1] = nx[j];
          parameters[2] = ny[j];
          parameters[3] = nz[j];
        }
      }
    }
    return;
  }
  Combinations(pool, need-1, chosen | (one<<at), at+1, set_of_pairs, parameters);
  Combinations(pool, need, chosen, at+1, set_of_pairs, parameters);
}

//___________________________________________________________________________________________________________
void BReconstructMuon::GetX0Y0(Double_t &x0, Double_t &y0, Double_t theta, Double_t phi) const
{
  // zero approximation for x0, y0 if theta and phi have already been defined

  x0 = 0;
  y0 = 0;
  int n = 0;
  for(int i = 0; i < fEvent->GetTotImpulses(); i++) {
    BImpulse * imp = fEvent->GetImpulse(i);
    Int_t nch = imp->GetChannelID();
    BChannelMask::EMode channelflag = fChannelMask->GetFlag(nch);

    if(fEventMask->GetOrigin(i)->GetFlag() == 1 && (channelflag == BChannelMask::kOn || channelflag == BChannelMask::kBadChargeCalib)) {
      Double_t xom = fGeomTel->At(nch)->GetX() - fXshift;
      Double_t yom = fGeomTel->At(nch)->GetY() - fYshift;
      Double_t zom = fGeomTel->At(nch)->GetZ() - fZshift;

      BMuon muon;
      Double_t p[6];
      p[0] = xom;
      p[1] = yom;
      p[2] = zom;
      p[3] = 0;
      p[4] = phi;
      p[5] = theta;
      muon.SetParameters(p);

      x0 += muon.DefineX0(0, 0, 0);
      y0 += muon.DefineY0(0, 0, 0);
      n++;
    }
  }
  x0 /= n;
  y0 /= n;
}
