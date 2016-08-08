#include "BRecQuality.h"

#include "TVector3.h"
#include "TMath.h"

ClassImp(BRecQuality);

BRecQuality::BRecQuality(const char *name, const char *title)
{
    fName  = name  ? name  : "BRecQuality";
    fTitle = title ? title : "Rec quality";

	fSourceAzimuthAngle = -1;
    fSourcePolarAngle = -1;

    fRecAzimuthAngle = -1;
    fRecPolarAngle = -1;
    
    fPsi = -1;
}
 
void BRecQuality::SetReadyToSave(Bool_t flag)
{
    MParContainer::SetReadyToSave(flag);
    
    CalcPsi();
}

void BRecQuality::CalcPsi()
{
	TVector3 source;
	//source.SetMagThetaPhi(1, fSourceAzimuthAngle, fSourcePolarAngle);
	source.SetMagThetaPhi(1, fSourcePolarAngle, fSourceAzimuthAngle);
	
	TVector3 rec;
	//rec.SetMagThetaPhi(1, fRecAzimuthAngle, fRecPolarAngle);
	rec.SetMagThetaPhi(1, TMath::Pi() - fRecPolarAngle, fRecAzimuthAngle + TMath::Pi());
	
	fPsi = source.Angle(rec);
}

void BRecQuality::SetRecAzimuthAngle(Float_t ang)
{ 
	Double_t pi2 = 2*TMath::Pi();
	
	if(ang < 0) {
		Long_t n = ang / pi2;
		ang += TMath::Abs(n) * pi2;
		ang = pi2 - TMath::Abs(ang);
	}
	
	if(ang >= 2*TMath::Pi()) {
		Long_t n = ang / pi2;
		ang -= n * pi2;
	}

	fRecAzimuthAngle = ang; 	
}

void BRecQuality::SetRecPolarAngle(Float_t ang)
{ 
	Double_t pi = TMath::Pi();
	
	ang = TMath::Abs(ang);
	
	if(ang > pi) {
		Long_t n = ang / pi;
		ang -= n * pi;
		ang = pi - ang;
	}

	fRecPolarAngle = ang; 	
}
