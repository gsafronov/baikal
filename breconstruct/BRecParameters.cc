#include "BRecParameters.h"
#include <iostream>
#include "TVector3.h" 
#include "TMath.h" 

ClassImp(BRecParameters)

//______________________________________________________________________
BRecParameters::BRecParameters(const char* name, const char* title)
{
    // Initialize name and title for MARS container
    
    fName = name ? name : "BRecParameters";
    fTitle = title ? title : "Reconstruct Parameters";

	fNumPar   = 0;
	fIsFixed  = 0;
	fInitial  = 0;
	fMinLimit = 0;
	fMaxLimit = 0;
	fStep     = 0;
	fErrors   = 0;
	
	fNhit = 0;
	fTres = 0;
	fNchGeom = 0;
	
	fNhitMC = 0;
	fTresMC = 0;
	fNchGeomMC = 0;
	fFuncValueMC = 0;
	
	fTimeRec = 0;
	fTimeMC = 0;
}

//______________________________________________________________________
void BRecParameters::SetNumPar(Int_t value)
{
	Clear();
	
	if(value > 0) {
		fNumPar = value;
		fIsFixed = new Int_t[value];
		fInitial = new Double_t[value];
		fMinLimit = new Double_t[value];
		fMaxLimit = new Double_t[value];
		fStep = new Double_t[value];
		fErrors = new Double_t[value];
	}	
}

//______________________________________________________________________
void BRecParameters::Clear(Option_t* o)
{
	if(fIsFixed) {
		delete [] fIsFixed;
		fIsFixed = 0;
	}
	if(fInitial) {
		delete [] fInitial;
		fInitial = 0;
	}	
	if(fMinLimit) {
		delete [] fMinLimit;
		fMinLimit = 0;
	}	
	if(fMaxLimit) {
		delete [] fMaxLimit;
		fMaxLimit = 0;
	}	
	if(fStep) {
		delete [] fStep;
		fStep = 0;
	}	
	if(fErrors) {
		delete [] fErrors;
		fErrors = 0;
	}	
	
	if(fTres) {
		delete [] fTres;
		fTres = 0;
	}	
	
	if(fNchGeom) {
		delete [] fNchGeom;
		fNchGeom = 0;
	}	
	
	if(fTresMC) {
		delete [] fTresMC;
		fTresMC = 0;
	}	
	
	if(fNchGeomMC) {
		delete [] fNchGeomMC;
		fNchGeomMC = 0;
	}		
	
	fNhit   = 0;
	fNhitMC   = 0;
	fNumPar = 0;
}

//______________________________________________________________________
BRecParameters::~BRecParameters() 
{
	Clear();
}

//______________________________________________________________________
void BRecParameters::SetErrors(const Double_t *value)
{
	for(int i = 0; i < fNumPar; i++) {
		fErrors[i] = value[i];
	}
}

//______________________________________________________________________
void BRecParameters::SetNhit(Int_t nhit) 
{ 
	fNhit = nhit; 
	if(fNhit > 0) {
		fTres = new Double_t[nhit];
		fNchGeom = new Int_t[nhit];
	}
}

//______________________________________________________________________
void BRecParameters::SetTresMC(Int_t nch, Double_t *tresmc, Int_t *nchgeom) 
{ 
	if(nch == 0) {
		return;
	}
	
	if(fNhitMC > 0) {
		if(fTresMC) {
			delete [] fTresMC;
			fTresMC = 0;
		}	
	
		if(fNchGeomMC) {
			delete [] fNchGeomMC;
			fNchGeomMC = 0;
		}		
	}
	
	fNhitMC = nch;
	fTresMC = new Double_t[fNhitMC];
	fNchGeomMC = new Int_t[fNhitMC];
	
	for(int i = 0; i < fNhitMC; i++) {
		fTresMC[i] = tresmc[i];
		fNchGeomMC[i] = nchgeom[i];
	}
}

//______________________________________________________________________
Double_t BRecParameters::GetTresNchGeom(Bool_t ismc, Int_t nchgeom) const
{
	for(int i = 0; i < fNhit; i++) {
		if(fNchGeom[i] == nchgeom) {
			if(ismc) {
				return fTresMC[i];			
			}
			else {
				return fTres[i];			
			}
		}
	}
	return -100000;
}

//______________________________________________________________________
Double_t BRecParameters::GetPsi() const
{
	Double_t c = TMath::Pi() / 180;
	
	TVector3 r1;
	r1.SetMagThetaPhi(1, fThetaRec * c, fPhiRec * c);
	
	TVector3 r2;
	r2.SetMagThetaPhi(1, fThetaMC * c, fPhiMC * c);
	
	return r1.Angle(r2);
}
