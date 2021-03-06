#include "BTimeClusters.h"
#include "BGeomTel.h"
#include "BEvent.h"
#include "BMCEvent.h"
#include "BExtractedImpulse.h"
#include "BExtractedImpulseTel.h"
#include "BJoinExtractedImpulseTel.h"
#include "BEventMask.h"
#include "BChannelMask.h"
#include "BExtractedHeader.h"

#include "BRecParameters.h"

#include "MParList.h"
#include "MLog.h"
#include "MLogManip.h"
 
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"

#include <math.h>
#include <vector>

ClassImp(BTimeClusters);

BTimeClusters::BTimeClusters(const char *name, const char *title)
{
  // chanIDs=chans;
  // fMaskNoise=maskNoise;
  // fInputMaskName ="";
  // fOutputMaskName="";

  fName  = name ? name : "BTimeClusters";
  fTitle = title ? title : "TimeClusterFilterMask";

  cVacuum=0.3;
  cWater=cVacuum/1.33;

  //define WP
  fSignalCut_hotspot1=2.5;
  fSignalCut_hotspot2=-1.5;
  fSafetyWindow=20;

  //gen-level cuts
  fGen_rhoCut=40;
  fGen_minAngle=10;
  fGen_maxAngle=90;

  fSignalCut_gen=0.1;
  fSeedSignalCut_gen=0.2;

  fVerbose=false;
  
  fEventCounter=0;
  fNoiseOMs=0;
  fSignalOMs=0;
  fSignalStrings=0;
  
  //  fOutputMaskName = "AmplitudeFilterMask";

  fOUT=new TFile("timeClustersDebug.root","RECREATE");
  h_1mu_hits_per_string_2pe=new TH1F("h_1mu_hits_per_string_2pe","h_1mu_hits_per_string_2pe",25,0,25);
  h_1mu_hits_2pe=new TH1F("h_1mu_hits_2pe","h_1mu_hits_2pe",200,0,200);
  h_1mu_fired_strings_2pe=new TH1F("h_1mu_fired_strings_2pe","h_1mu_fired_strings_2pe",9,0,9);
  hreco_hits_per_string=new TH1F("hreco_hits_per_string","hreco_hits_per_string",25,0,25);
  hreco_hits=new TH1F("hreco_hits","hreco_hits",200,0,200);
  hreco_fired_strings=new TH1F("hreco_fired_strings","hreco_fired_strings",9,0,9);
  h_strings_reco_vs_gen=new TH2F("h_strings_reco_vs_gen","h_strings_reco_vs_gen",9,0,9,9,0,9);
  h_strings_bevt_vs_gen=new TH2F("h_strings_bevt_vs_gen","h_strings_bevt_vs_gen",9,0,9,9,0,9);

  h_hitsPerString_reco_vs_gen=new TH2F("h_hitsPerString_reco_vs_gen","h_hitsPerString_reco_vs_gen",25,0,25,25,0,25);
  h_hitsPerString_bevt_vs_gen=new TH2F("h_hitsPerString_bevt_vs_gen","h_hitsPerString_bevt_vs_gen",25,0,25,25,0,25);

  h_smuons_polar_vs_rho=new TH2F("h_smuons_polar_vs_rho","h_smuons_polar_vs_rho",10,0,100,25,0,250);
  h_rcand_polar_vs_rho=new TH2F("h_rcand_polar_vs_rho","h_rcand_polar_vs_rho",10,0,100,25,0,250);
  h_strings_polar_vs_rho=new TProfile2D("h_strings_polar_vs_rho","h_strings_polar_vs_rho",10,0,100,25,0,250,0,10000);

  float energy_bins[37]={0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,2000,2200,2400,2600,2800,3000,3250,3500,3750,4000,4500,5000,5500,6000,7000,8000,9000,10000};
    
  h_muon_energy=new TH1F("h_muon_energy","h_muon_energy",36,energy_bins);
  h_muon_energy_rcand_gen=new TH1F("h_muon_energy_rcand_gen","h_muon_energy_rcand_gen",36,energy_bins);
  h_muon_energy_rcand=new TH1F("h_muon_energy_rcand","h_muon_energy_rcand",36,energy_bins);
  
  h_hitSignal_reco_vs_gen=new TH2F("h_hitSignal_reco_vs_gen","h_hitSignal_reco_vs_gen",22,-1,10,22,-1,10);

  h_noiseFrac_clustered=new TH2F("h_noiseFrac_clustered","h_noiseFrac_clustered",20,0.5,10.5,20,0,200);
  h_signalFrac_clustered=new TH2F("h_signalFrac_clustered","h_signalFrac_clustered",20,0.5,10.5,20,0,200);
  h_stringSignalFrac_clustered=new TH2F("h_stringSignalFrac_clustered","h_stringSignalFrac_clustered",20,0.5,10.5,20,0,200);
  
  h_clusteredFrac_noise=new TH2F("h_clusteredFrac_noise","h_clusteredFrac_noise",20,0.5,10.5,20,0,200);
  h_clusteredFrac_signal=new TH2F("h_clusteredFrac_signal","h_clusteredFrac_signal",20,0.5,10.5,20,0,200);

  h_effpur_mult=new TH2F("h_effpur_mult","h_effpur_mult",20,0.5,10.5,20,0,200);
  
  h_clustered=new TH2F("h_clustered","h_clustered",20,0.5,10.5,20,0,200);
  h_clustered_highMult=new TH2F("h_clustered_highMult","h_clustered_highMult",20,0.5,10.5,20,0,200);
  
  hWhatIsOM=new TH1F("hWhatIsOM","hWhatIsOM",3,0,3);

  hitMap_gen=new TH2F("hitMap_gen","hitMap_gen",9,0,9,25,0,25);
  hitMap_det=new TH2F("hitMap_det","hitMap_det",9,0,9,25,0,25);
}


BTimeClusters::~BTimeClusters()
{
}

Int_t BTimeClusters::PreProcess(MParList * pList)
{
  fMCEvent=(BMCEvent*)pList->FindObject("BMCEvent","BMCEvent");
  
  fEvent=(BEvent*)pList->FindObject("BEvent","BEvent");
  if (!fEvent){
    * fLog << err << AddSerialNumber("BEvent") << " not found... aborting." << endl;
    return kFALSE;
  }
  
  fChannelMask = (BChannelMask*)pList->FindObject("BChannelMask");
  if (!fChannelMask)
    {
      * fLog << err << AddSerialNumber("BChannelMask") << " not found... aborting." << endl;
      return kFALSE;
    }
  
  fInputEventMask = (BEventMask*)pList->FindObject(fInputMaskName);
  if (!fInputEventMask)
    {
      * fLog << err << AddSerialNumber(fInputMaskName) << " not found... aborting." << endl;
      return kFALSE;
    }
  
  fGeomTel = (BGeomTel*)pList->FindObject(AddSerialNumber("BGeomTel"));

  /*
  std::cout<<"geometry in filter:    "<<std::endl;
  for (int i=0; i<192; i++){
    float z=fGeomTel->At(i)->GetZ();
    std::cout<<"chan :"<<i<<"     Z: "<<z<<std::endl;
  }
  */
  /*
  std::cout<<"BChanSetter: exclude calibrated channel from reconstruction"<<std::endl;
  for (int i=0; i<chanIDs.size(); i++){
    if (fChannelMask->GetFlag(chanIDs[i])==1) fChannelMask->SetFlag(chanIDs[i], BChannelMask::EMode::kOff);
  }
  */
  if (!BFilter::PreProcess(pList)) {
    return kFALSE;
  }

 
  
  return kTRUE;
}


Bool_t BTimeClusters::Filter()
{
  fEventCounter++;
  if (fEventCounter%10000==0) std::cout<<fEventCounter<<std::endl;
  int n_impulse_channels=fMCEvent->GetChannelN();
  /*
  for (int i=0; i<fEvent->GetTotImpulses(); i++){
    bool isNoise=false;

    for (int j=0; j<n_impulse; j++){
      BMCHitChannel* fHitChan=fMCEvent->GetHitChannel(j);
      for (int k=0; k<fHitChan->GetPulseN(); k++){
	if (fEvent->GetImpulse(i)->GetTime()==fHitChan->GetPulse(k)->GetTime()){
	  if (fHitChan->GetPulse(k)->GetMagic()==1) isNoise=true;
	  fOutputEventMask->SetFlag(i,0);
	}
      }
    }
    if (isNoise) fOutputEventMask->SetFlag(i,0);
    else fOutputEventMask->SetFlag(i,fInputEventMask->GetFlag(i));
  }
    
  */
  
  ///////////////////////////////////////////////////////////////////
  //check gen-level, fill histos
  //  std::cout<<"eprst  "<<n_impulse<<std::endl;
  //find direction and rho of the muon
  float primThetaRad=M_PI*fMCEvent->GetPrimaryParticlePolar()/180;
  float primPhiRad=M_PI*fMCEvent->GetPrimaryParticleAzimuth()/180;

  TVector3 genVec(sin(primThetaRad)*cos(primPhiRad),
		  sin(primThetaRad)*sin(primPhiRad),
		  cos(primThetaRad));

  TVector3 inPoTrack(fMCEvent->GetTrack(0)->GetX()-500*genVec.X(), fMCEvent->GetTrack(0)->GetY()-500*genVec.Y(), fMCEvent->GetTrack(0)->GetZ()-500*genVec.Z());

  TVector3 zero(0,0,0);
  float rho=getTrackDistanceToOM(inPoTrack, genVec, zero);

  
  //set the event mask
  int impulse_n=fEvent->GetTotImpulses();
  for (int iPulse=0; iPulse<impulse_n; iPulse++){
    fOutputEventMask->SetFlag(iPulse,0);
  }
  
  if (fMCEvent->GetResponseMuonsN()!=1) return kFALSE;
  if (rho>fGen_rhoCut) return kFALSE;
  if (180-fMCEvent->GetPrimaryParticlePolar()>fGen_maxAngle) return kFALSE;
  if (180-fMCEvent->GetPrimaryParticlePolar()<fGen_minAngle) return kFALSE;
  //  if (

  //  if (fMCEvent->GetTrack(0)->GetMuonEnergy()<1000||
  //  if (rho>20) return kFALSE;
  /////////////////////////

  //find trigger string and exclude it from analysis
  int trigger_string=findTriggerString(4,1.5);
  //int trigger_string=-1;
  
  std::vector<int> string_impulses_gen[8];
  std::vector<int> hasGaps;
  for (int i=0; i<8; i++) hasGaps.push_back(false);
  
  int noisePulseID=-100;
  int nHits=0;
  bool hasLargeHit=false;
  std::vector<bool> stringHasLargeHit;
  for (int i=0; i<8; i++) stringHasLargeHit.push_back(false);
  
  float maxAmplitude_gen=-1;
  int maxAmplitude_gen_chanID=-1;

  std::vector<float> singlePulse_ampls;
  std::vector<int> singlePulse_ids;

  std::vector<bool> noiseOMs_gen;
  std::vector<bool> signalOMs_gen;
  std::vector<bool> signalStrings_gen;
  std::vector<bool> isMatched;
  for (int iPulse=0; iPulse<impulse_n; iPulse++) isMatched.push_back(false);

  for (int i=0; i<192; i++) {
    noiseOMs_gen.push_back(false);
    signalOMs_gen.push_back(false);
  }

  for (int i=0; i<8; i++){
    signalStrings_gen.push_back(false);
  }    
  
  for (int i=0; i<n_impulse_channels; i++) {
    //    std::cout<<i<<std::endl;
    BMCHitChannel* fHitChan=fMCEvent->GetHitChannel(i);
    int idch=fHitChan->GetChannelID()-1;
    idch=24*floor(idch/24)+(24-idch%24);
    idch=idch-1;
    int iString=floor(idch/24);
    //exclude the trigger string
    if (iString==trigger_string) continue;
    hitMap_gen->Fill(iString, idch%24,1);
    float signal=0;
    float time=-1;
    bool isNoise=false;
    bool isSignal=false;
    for (int k=0; k<fHitChan->GetPulseN(); k++){
      if (fHitChan->GetPulse(k)->GetMagic()!=1) {
	signal+=fHitChan->GetPulse(k)->GetAmplitude();
	time=fHitChan->GetPulse(k)->GetTime();
	isSignal=true;
      }
      if (fHitChan->GetPulse(k)->GetMagic()==1) isNoise=true;
    }

    if ((isNoise&&!isSignal)) {
      noiseOMs_gen[idch]=true;
      fNoiseOMs++;
    }

    if (isSignal) {
      if (signal>fSignalCut_gen) signalOMs_gen[idch]=true;
      for (int k=0; k<impulse_n; k++){
	if (fabs(time-fEvent->GetImpulse(k)->GetTime())<30) {
	  isMatched[k]=true;;
	}
      }
      //fSignalOMs++;
    }
    
    if (isNoise&&!isSignal) hWhatIsOM->Fill(0.1,1);
    if (isNoise&&isSignal) hWhatIsOM->Fill(1.1,1);
    if (!isNoise&&isSignal) hWhatIsOM->Fill(2.1,1);
    
    //for gen-reco amplitude comparison
    //    if (fHitChan->GetPulseN()==1&&fHitChan->GetPulse(0)->GetMagic()!=1){
    if (signal!=0){
      singlePulse_ampls.push_back(fHitChan->GetPulse(0)->GetAmplitude());
      singlePulse_ids.push_back(idch);
    }
    //    std::cout<<"impulse "<<i<<"   iString: "<<iString<<"   idch "<<idch<<"   bb  "<<fHitChan->GetChannelID()<<std::endl;
    if (signal>=fSignalCut_gen) string_impulses_gen[iString].push_back(i);
    
    
    if (signal>=fSeedSignalCut_gen) {
      hasLargeHit=true;
      stringHasLargeHit[iString]=true;
    }

    if (signal>maxAmplitude_gen){
      maxAmplitude_gen=signal;
      maxAmplitude_gen_chanID=idch;
    }
  }

  /*
  //SET MASK TO MATCHED IMPULSES/////////////////////////
  int nMatched=0;
  int nStr=0;
  std::vector<int> stringHasHits;
  for (int i=0; i<8; i++) stringHasHits.push_back(0);
  for (int i=0; i<impulse_n; i++){
    if (isMatched[i]) {
      nMatched++;
      int ch=fEvent->GetImpulse(i)->GetChannelID();
      stringHasHits[int(floor(ch/24))]++;
    }
  }
  int nStrings=0;
  for (int i=0; i<8; i++){
    if (stringHasHits[i]>1) nStr++;
  }
  if (nStr>=3&&nMatched>=6){
    for (int i=0; i<impulse_n; i++){
      if (isMatched[i]) fOutputEventMask->SetFlag(i,1);
    }
    return kTRUE;
  }
  else return kFALSE;  
  
  ////////////////////////////////////////////////////////
  */
  
  for (int iString=0; iString<8; iString++) {
    int nGaps=countMCStringGaps(string_impulses_gen[iString]);
    if (nGaps>0) {
      hasGaps[iString]=true;
      string_impulses_gen[iString].clear();
    }
    if (!hasGaps[iString]&&stringHasLargeHit[iString]) nHits+=string_impulses_gen[iString].size();
  }

 
  //  if (!hasLargeHit) return kFALSE;
  
  h_muon_energy->Fill(fMCEvent->GetTrack(0)->GetMuonEnergy(),1);
  
  if (nHits>=5) h_rcand_polar_vs_rho->Fill(180-fMCEvent->GetPrimaryParticlePolar(),rho,1);
  h_smuons_polar_vs_rho->Fill(180-fMCEvent->GetPrimaryParticlePolar(),rho,1);
  
  //  std::cout<<"popo"<<std::endl;
  h_1mu_hits_2pe->Fill(nHits,1);

  
  int firedStrings=0;
  for (int iString=0; iString<8; iString++){
    if (hasGaps[iString]) continue;
    h_1mu_hits_per_string_2pe->Fill(string_impulses_gen[iString].size(),1);
    if (string_impulses_gen[iString].size()>1&&stringHasLargeHit[iString]) {
      fSignalStrings++;
      firedStrings++;
      signalStrings_gen[iString]=true;
      if (string_impulses_gen[iString].size()>=3) {
	//signalStrings_gen[iString]=true;
	fSignalOMs+=string_impulses_gen[iString].size();
      }
    }
    //    if string_impulses_gen[iString].size()std::cout<<"gen string:  "<<iString<<":   "<<string_impulses_gen[iString].size()<<" hits"<<std::endl;
       
  }

  if (nHits>=5&&firedStrings>=2) h_muon_energy_rcand_gen->Fill(fMCEvent->GetTrack(0)->GetMuonEnergy(),1);
  
  //  if (firedStrings<2) return kFALSE;
  
  h_1mu_fired_strings_2pe->Fill(firedStrings,1);
  
  h_strings_polar_vs_rho->Fill(180-fMCEvent->GetPrimaryParticlePolar(),rho,firedStrings,1);

  //  std::cout<<"zoob"<<std::endl;
  //detector-level impulses, combine in strings

  ////////////////////////////////////////////////////////////////ACTUAL RECO
  
  std::vector<int> string_impulses[8];

  //  std::cout<<"aa"<<std::endl;

  float maxAmplitude_rec;
  int maxAmplitude_rec_chanID;
  
  int nPulses_initial=0;
  
  for (int iPulse=0; iPulse<impulse_n; iPulse++){
    if (fInputEventMask->GetFlag(iPulse)==0) continue;
    int iChannel=fEvent->GetImpulse(iPulse)->GetChannelID();
    //exclude trigger string:
    int iString=int(floor(iChannel/24));
    if (iString==trigger_string) continue;
    
    nPulses_initial++;
    
    hitMap_det->Fill(iString, iChannel%24,1);
    
    string_impulses[iString].push_back(iPulse);
    
    for (int iGen=0; iGen<singlePulse_ids.size(); iGen++){
      if (iChannel==singlePulse_ids[iGen]&&stringHasLargeHit[iString]) h_hitSignal_reco_vs_gen->Fill(singlePulse_ampls[iGen], fEvent->GetImpulse(iPulse)->GetAmplitude(),1); 
    }
    if (fEvent->GetImpulse(iPulse)->GetAmplitude()>maxAmplitude_rec){
      maxAmplitude_rec=fEvent->GetImpulse(iPulse)->GetAmplitude();
      maxAmplitude_rec_chanID=iChannel;
    }
    
  }

  // if (maxAmplitude_rec_chanID==maxAmplitude_gen_chanID) h_hitSignal_reco_vs_gen->Fill(maxAmplitude_gen,  maxAmplitude_rec,1);
  //else h_hitSignal_reco_vs_gen->Fill(-1,-1,1);
  
  //std::cout<<"bb"<<std::endl;
  std::vector<std::vector<int> > stringClusters[8];
  int hitStrings=0;
  if (fVerbose) std::cout<<">>>>>>>>>>>>>>>>>>>>>>>>>NEXT EVENT"<<std::endl;
  for (int iString=0; iString<8; iString++){
    
    //do not do anything is there is a gap in MC signal
    if (hasGaps[iString]||iString==trigger_string) continue;
    /////////
    
    if (string_impulses[iString].size()>0) hitStrings++;
    if (fVerbose){
      if (string_impulses_gen[iString].size()>0&&firedStrings>=0){
	std::cout<<"gen string:  "<<iString<<":   "<<string_impulses_gen[iString].size()<<" hits"<<std::endl;
	for (int ip=0; ip<string_impulses_gen[iString].size(); ip++){
	  int  idch=fMCEvent->GetHitChannel(string_impulses_gen[iString][ip])->GetChannelID()-1;
	  idch=24*floor(idch/24)+(24-idch%24);
	  idch=idch-1;
	  std::cout<<"    impulse idch: "<<idch<<std::endl;
	}
	stringClusters[iString]=buildStringClusters(iString, string_impulses[iString]);
      }
    }
    else stringClusters[iString]=buildStringClusters(iString, string_impulses[iString]);
  }

  //  std::cout<<"cc"<<std::endl;
  
  //build global cluster [SIMPLE]
  int usedStrings=0;
  
  std::vector<int> globalCluster=buildGlobalCluster(stringClusters);

  for (int iString=0; iString<8; iString++){
    if (hasGaps[iString]) continue;
    if (stringClusters[iString].size()>0) usedStrings++;
    hreco_hits_per_string->Fill(stringClusters[iString].size(),1);
    if (stringHasLargeHit[iString]) h_hitsPerString_reco_vs_gen->Fill(string_impulses_gen[iString].size(),stringClusters[iString].size(),1);
  }
  if (fVerbose) std::cout<<"firedStrings: "<<firedStrings<<"     recoStrings: "<<usedStrings<<std::endl;
  if (usedStrings>=2&&globalCluster.size()>=5) h_muon_energy_rcand->Fill(fMCEvent->GetTrack(0)->GetMuonEnergy(),1);
  
  hreco_fired_strings->Fill(usedStrings,1);

  if (globalCluster.size()>=5&&usedStrings>=3){
    //    std::cout<<"GOTCHA"<<std::endl;
    for (int iClustered=0; iClustered<globalCluster.size(); iClustered++){
      fOutputEventMask->SetFlag(globalCluster[iClustered],1);
    }
    for (int i=0; i<impulse_n; i++){
      std::cout<<"impulse:  "<<i<<"   flag: "<<fOutputEventMask->GetFlag(i)<<std::endl;
    }
    
  }
  
  h_strings_reco_vs_gen->Fill(firedStrings, usedStrings,1);
  h_strings_bevt_vs_gen->Fill(firedStrings, hitStrings,1);
  hreco_hits->Fill(globalCluster.size(),1);
    //std::cout<<"OLD METH STRINGS: "<<usedStrings<<std::endl;
  //  std::cout<<"initial pulses: "<<nPulses_initial<<"   after clusterisation: "<<globalCluster.size()<<" used strings:   "<<usedStrings<<"     mc response muons:   "<<fMCEvent->GetResponseMuonsN()<<std::endl;

  //  runWPscan(string_impulses, noiseOMs_gen, signalOMs_gen, signalStrings_gen);
  
  if (globalCluster.size()>=5&&usedStrings>=3) return kTRUE;
  else return kFALSE;
}

std::vector<std::vector<int> > BTimeClusters::buildStringClusters(int iString, std::vector<int> string_impulses)
{
  // std::cout<<"BUILD"<<std::endl;
  std::vector<int> stringCluster;
  std::vector<int> isClustered;
  std::vector<std::vector<int> > result;
  
  if (string_impulses.size()==0) return result;				   
  
  ///////////////////TEMPORARY
  float primThetaRad=M_PI*fMCEvent->GetPrimaryParticlePolar()/180;
  float primPhiRad=M_PI*fMCEvent->GetPrimaryParticleAzimuth()/180;

  TVector3 genVec(sin(primThetaRad)*cos(primPhiRad),
		  sin(primThetaRad)*sin(primPhiRad),
		  cos(primThetaRad));

  TVector3 inPoTrack(fMCEvent->GetTrack(0)->GetX()-500*genVec.X(), fMCEvent->GetTrack(0)->GetY()-500*genVec.Y(), fMCEvent->GetTrack(0)->GetZ()-500*genVec.Z());

  TVector3 zero(0,0,0);
  float rho=getTrackDistanceToOM(inPoTrack, genVec, zero);
  
   std::vector<std::vector<int>> hot_spots=findHotSpots(string_impulses);
  
  //  if (fMCEvent->GetResponseMuonsN()==1&&rho<fGen_rhoCut){
  //    std::cout<<"string #"<<iString+1<<"   spots found: "<<hot_spots.size()<<std::endl;
  //    for (int i=0; i<hot_spots.size(); i++){
  //      std::cout<<"       spot #"<<i<<"  size "<<hot_spots[i].size()<<std::endl;
  //    }
  // }
    
  std::vector<int> max_size_cluster;
  int max_cluster_size=0;

  if (fVerbose) std::cout<<"string: "<<iString<<"    "<<hot_spots.size()<<" hotspots found,      size of string: "<<string_impulses.size()<<std::endl;

  for (int i=0; i<hot_spots.size(); i++){
    std::vector<int> stringCluster=addImpulses(hot_spots[i], string_impulses);

    /*
    //remove clustered impulses from string_impulses and hot_spots
    for (int j=0; j<string_impulses.size(); j++){
      bool drop=false;
      for (int k=0; k<stringCluster.size(); k++){
	if (string_impulses[j]==stringCluster[k]) drop=true;
      }
      if (drop) string_impulses.erase(string_impulses.begin()+j);
    }
    */
    //remove clustered impulses from hotspots
    for (int j=0; j<hot_spots.size(); j++){
      bool drop=false;
      for (int q=0; q<stringCluster.size(); q++){
	if (string_impulses[hot_spots[j][0]]==stringCluster[q]) drop=true;
	if (string_impulses[hot_spots[j][1]]==stringCluster[q]) drop=true;
      }
      if (drop) hot_spots.erase(hot_spots.begin()+j);
    }
    
    result.push_back(stringCluster);
    
    //stop if hotspots have been clustered
    if (i>=hot_spots.size()) break;

    /*
    //check if the cluster is the biggest in the event
    if (stringCluster.size()>max_cluster_size){
      max_cluster_size=stringCluster.size();
      max_size_cluster=stringCluster;
    }
    */
  }
  return result;
  
}

std::vector<int> BTimeClusters::addImpulses(std::vector<int> hotspot, std::vector<int> string_impulses)
//						   int max_ampl_id, int adjacent_id, std::vector<int> string_impulses, float window)
{
  std::vector<int> stringCluster;
  std::vector<int> storeys_pulse_id;
  std::vector<float> storeys_pulse_time;
  
  for (int i=0; i<24; i++){
    storeys_pulse_id.push_back(-1);
    storeys_pulse_time.push_back(-1);
  }

  stringCluster.push_back(string_impulses[hotspot[0]]);
  stringCluster.push_back(string_impulses[hotspot[1]]);
  
  for (int i=0; i<stringCluster.size(); i++){
    storeys_pulse_id[(fEvent->GetImpulse(stringCluster[i])->GetChannelID())%24]=stringCluster[i];
    storeys_pulse_time[(fEvent->GetImpulse(stringCluster[i])->GetChannelID())%24]=fEvent->GetImpulse(stringCluster[i])->GetTime();
  }
  
  /*
  stringCluster.push_back(string_impulses[max_ampl_id]);
  stringCluster.push_back(string_impulses[adjacent_id]);
  fOutputEventMask->SetFlag(string_impulses[max_ampl_id],10);
  fOutputEventMask->SetFlag(string_impulses[adjacent_id],10);
  */
  for (int iSeed=0; iSeed<hotspot.size(); iSeed++){

    float seed_time=fEvent->GetImpulse(string_impulses[hotspot[iSeed]])->GetTime();
    int seed_channel_id=fEvent->GetImpulse(string_impulses[hotspot[iSeed]])->GetChannelID();
    int seed_storey=(seed_channel_id)%24;
    float adjacent_time=fEvent->GetImpulse(string_impulses[hotspot[1-iSeed]])->GetTime();
    int adjacent_channel_id=fEvent->GetImpulse(string_impulses[hotspot[1-iSeed]])->GetChannelID();
    
    float deltaZ=fabs(fGeomTel->At(seed_channel_id)->GetZ()-fGeomTel->At(adjacent_channel_id)->GetZ());
    float deltaT=seed_time-adjacent_time;
    int id_increment=seed_channel_id-adjacent_channel_id;

    //add +-1 channels
    float time_early=seed_time+deltaT-fSafetyWindow; //experimental
    float time_late=max(seed_time+deltaZ/cWater+fSafetyWindow, seed_time+deltaT+fSafetyWindow);

    if (fVerbose) std::cout<<" deltaT: "<<deltaT<<"   seed_time: "<<seed_time<<"    seed chan: "<<seed_channel_id<<"   t_early: "<<time_early<<"   t_late: "<<time_late<<std::endl;
    
    for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
      if (iPulse==hotspot[iSeed]||iPulse==hotspot[1-iSeed]) continue;
      int chanID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
      float chan_time=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
      if (fVerbose) std::cout<<"time: "<<chan_time<<"       chanID: "<<chanID<<std::endl;
      if (chanID==seed_channel_id+id_increment){
	if (chan_time>time_early&&chan_time<time_late){
	  stringCluster.push_back(string_impulses[iPulse]);
	  storeys_pulse_id[chanID%24]=string_impulses[iPulse];
	  storeys_pulse_time[chanID%24]=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
	}
      }
    }
  
    if (fVerbose){
      if (stringCluster.size()==2) std::cout<<"fail"<<std::endl;
      if (stringCluster.size()==0) std::cout<<"super fail"<<std::endl;
      if (stringCluster.size()>2) std::cout<<"ok"<<std::endl;
    }
    
    
    //add up to 7 channels
    
    int id[]={1,3};
    for (int k=0; k<2; k++){
      int i=id_increment*id[k];
      if (storeys_pulse_id[seed_storey+i]==-1) continue;
      for (int j=1; j<3; j++){
	//find min and max time for new channel
	if (seed_storey+i+sgn(i)*j<0||seed_storey+i+sgn(i)*j>23) continue;
	float tim0=storeys_pulse_time[seed_storey+i]+j*(storeys_pulse_time[seed_storey+i]-storeys_pulse_time[seed_storey+i-sgn(i)]);
	float tim1=storeys_pulse_time[seed_storey+i]+j*(storeys_pulse_time[seed_storey+i]-storeys_pulse_time[seed_storey+i-2*sgn(i)])/2;
	float tim2=storeys_pulse_time[seed_storey+i-sgn(i)]+(j+1)*(storeys_pulse_time[seed_storey+i-sgn(i)]-storeys_pulse_time[seed_storey+i-2*sgn(i)]);
	
	float min_estimate=min(tim0,tim1);
	min_estimate=min(min_estimate,tim2);
	
	float max_estimate=max(tim0,tim1);
	max_estimate=max(max_estimate,tim2);
	//loop over string pulses
	for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
	  //if (fOutputEventMask->GetFlag(string_impulses[iPulse])==10) continue;
	  int chanID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
	  if (chanID!=seed_channel_id+i+sgn(i)*j) continue;
	  float time=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
	  if (time<max_estimate+fSafetyWindow&&time>min_estimate-fSafetyWindow){
	    stringCluster.push_back(string_impulses[iPulse]);
	    storeys_pulse_id[seed_storey+i+sgn(i)*j]=string_impulses[iPulse];
	    storeys_pulse_time[seed_storey+i+sgn(i)*j]=string_impulses[iPulse];
	    //fOutputEventMask->SetFlag(string_impulses[iPulse],10);
	  } 
	}
      }
    }
    
  }

  if (fVerbose){
    std::cout<<"size of assembled cluster: "<<stringCluster.size()<<std::endl;
    for (int i=0; i<stringCluster.size(); i++){
      std::cout<<"         chanID: "<<fEvent->GetImpulse(stringCluster[i])->GetChannelID()<<std::endl;
    }
  }
  return stringCluster;  
}


Int_t BTimeClusters::PostProcess()
{
  //normalise 2d strings response matrix
  for (int i=0; i<9; i++){
    double normFactor=h_1mu_fired_strings_2pe->GetBinContent(i+1);
    if (normFactor==0) continue;
    for (int j=0; j<9; j++){
      double bico=(double)h_strings_reco_vs_gen->GetBinContent(i+1,j+1)/h_1mu_fired_strings_2pe->GetBinContent(i+1);
      //      std::cout<<i<<"   "<<j<<"   "<<"   "<<(double)normFactor<<"   "<<"   "<<(double)bico<<"   "<<(double)pow((double)normFactor,-1)<<std::endl;
      h_strings_reco_vs_gen->SetBinContent(i+1, j+1, (double)bico);
    }      
  }

    //normalise 2d strings response matrix
  for (int i=0; i<9; i++){
    double normFactor=h_1mu_fired_strings_2pe->GetBinContent(i+1);
    if (normFactor==0) continue;
    for (int j=0; j<9; j++){
      double bico=(double)h_strings_bevt_vs_gen->GetBinContent(i+1,j+1)/h_1mu_fired_strings_2pe->GetBinContent(i+1);
      //      std::cout<<i<<"   "<<j<<"   "<<"   "<<(double)normFactor<<"   "<<"   "<<(double)bico<<"   "<<(double)pow((double)normFactor,-1)<<std::endl;
      h_strings_bevt_vs_gen->SetBinContent(i+1, j+1, (double)bico);
    }      
  }

  for (int i=0; i<25; i++){
    double normFactor=h_hitsPerString_reco_vs_gen->ProjectionX()->GetBinContent(i+1);
    if (normFactor==0) continue;
    for (int j=0; j<25; j++){
      double bico=(double)h_hitsPerString_reco_vs_gen->GetBinContent(i+1,j+1)/normFactor;
      //      std::cout<<i<<"   "<<j<<"   "<<"   "<<(double)normFactor<<"   "<<"   "<<(double)bico<<"   "<<(double)pow((double)normFactor,-1)<<std::endl;
      h_hitsPerString_reco_vs_gen->SetBinContent(i+1, j+1, (double)bico);
    }      
  }

  for (int i=0; i<22; i++){
    double normFactor=h_hitSignal_reco_vs_gen->ProjectionX()->GetBinContent(i+1);
    if (normFactor==0) continue;
    for (int j=0; j<22; j++){
      double bico=(double)h_hitSignal_reco_vs_gen->GetBinContent(i+1,j+1)/normFactor;
      //      std::cout<<i<<"   "<<j<<"   "<<"   "<<(double)normFactor<<"   "<<"   "<<(double)bico<<"   "<<(double)pow((double)normFactor,-1)<<std::endl;
      h_hitSignal_reco_vs_gen->SetBinContent(i+1, j+1, (double)bico);
    }      
  }

  h_noiseFrac_clustered->Add(h_clusteredFrac_noise,1);
  h_noiseFrac_clustered->Scale(pow(fNoiseOMs,-1));

  h_signalFrac_clustered->Add(h_clusteredFrac_signal,1);
  h_signalFrac_clustered->Scale(pow(fSignalOMs,-1));

  h_clusteredFrac_noise->Divide(h_clusteredFrac_noise,h_clustered,1,1);
  h_clusteredFrac_signal->Divide(h_clusteredFrac_signal,h_clustered_highMult,1,1);
  
  h_effpur_mult->Add(h_clusteredFrac_signal,1);
  h_effpur_mult->Multiply(h_signalFrac_clustered);

  h_stringSignalFrac_clustered->Scale(pow(fSignalStrings,-1));
  
  h_rcand_polar_vs_rho->Divide(h_rcand_polar_vs_rho,h_smuons_polar_vs_rho,1,1);

  h_muon_energy->TH1F::Sumw2();
  h_muon_energy_rcand->TH1F::Sumw2();
  h_muon_energy_rcand->Divide(h_muon_energy_rcand, h_muon_energy,1,1,"B");

  h_muon_energy_rcand_gen->TH1F::Sumw2();
  h_muon_energy_rcand_gen->Divide(h_muon_energy_rcand_gen, h_muon_energy,1,1,"B");
  
  fOUT->cd();
  h_clusteredFrac_noise->Write();
  h_clusteredFrac_signal->Write();
  h_noiseFrac_clustered->Write();
  h_signalFrac_clustered->Write();
  h_stringSignalFrac_clustered->Write();
  h_effpur_mult->Write();
  h_1mu_hits_per_string_2pe->Write();
  h_1mu_hits_2pe->Write();
  h_1mu_fired_strings_2pe->Write();
  hreco_hits_per_string->Write();
  hreco_fired_strings->Write();
  hreco_hits->Write();
  h_strings_reco_vs_gen->Write();
  h_strings_bevt_vs_gen->Write();
  h_smuons_polar_vs_rho->Write();
  h_rcand_polar_vs_rho->Write();
  h_strings_polar_vs_rho->Write();
  h_hitsPerString_reco_vs_gen->Write();
  h_muon_energy->Write();
  h_muon_energy_rcand->Write();
  h_muon_energy_rcand_gen->Write();
  h_hitSignal_reco_vs_gen->Write();
  hWhatIsOM->Write();
  hitMap_gen->Write();
  hitMap_det->Write();
  fOUT->Close();
  return kTRUE;
}

float BTimeClusters::getTrackDistanceToOM(TVector3 initialPoint, TVector3 direction, TVector3 xyzOM)
{
  TVector3 diff(xyzOM.X()-initialPoint.X(),xyzOM.Y()-initialPoint.Y(),xyzOM.Z()-initialPoint.Z());

  float absDir=direction.Mag();

  float scalarProd=diff.Dot(direction);

  float cosAlpha=scalarProd/(diff.Mag()*absDir);

  //cosAlpha=round(cosAlpha*10000000)*pow(10000000,-1); 

  if (cosAlpha*cosAlpha>1) cosAlpha=1;
  
  float dist=diff.Mag()*sqrt(1-cosAlpha*cosAlpha);

  return dist;
}

int BTimeClusters::runWPscan(std::vector<int>* string_impulses, std::vector<bool> noiseOMs_gen, std::vector<bool> signalOMs_gen, std::vector<bool> signalStrings_gen)
{
  std::vector<std::vector<int> > stringClusters[8];
  int hitStrings=0;
  float fSignalCut_hotspot1_buf=fSignalCut_hotspot1;
  float fSafetyWindow_buf=fSafetyWindow;
  for (int i=0; i<20; i++){
    fSignalCut_hotspot1=0.5+0.5*i;
    for (int j=0; j<20; j++){
      fSafetyWindow=10*j;
      std::vector<int> clusteredOMs;
      std::vector<bool> omIsClustered;
      for (int q=0; q<192; q++) omIsClustered.push_back(false);
      int nSignalOMs[8]={0};
      for (int iString=0; iString<8; iString++){
	//if (!signalStrings_gen
	//find how many signal OMs
	//	nSignalOMs[iString]=0;
	for (int q=0; q<signalOMs_gen.size(); q++){
	  if (signalOMs_gen[q]&&int(floor(q/24))==iString&&signalStrings_gen[iString]) nSignalOMs[iString]++;
	}

	//	if (string_impulses[iString].size()>0) hitStrings++;
	stringClusters[iString]=buildStringClusters(iString, string_impulses[iString]);

	if (stringClusters[iString].size()>0&&signalStrings_gen[iString]) h_stringSignalFrac_clustered->Fill(fSignalCut_hotspot1, fSafetyWindow,1);

	for (int k=0; k<stringClusters[iString].size(); k++){
	  for (int p=0; p<stringClusters[iString][p].size(); p++){
	    int omID=fEvent->GetImpulse(stringClusters[iString][p][k])->GetChannelID();
	    //	  std::cout<<omID<<std::endl;
	    omIsClustered[omID]=true;
	  }
	}

	if (stringClusters[iString].size()>1) hitStrings++;
	
      }
      for (int iOM=0; iOM<192; iOM++){
	if (hitStrings>0&&omIsClustered[iOM]){
	  h_clustered->Fill(fSignalCut_hotspot1, fSafetyWindow,1);
	  if (noiseOMs_gen[iOM]) h_clusteredFrac_noise->Fill(fSignalCut_hotspot1, fSafetyWindow,1);
	}
	if (signalStrings_gen[int(floor(iOM/24))]&&nSignalOMs[int(floor(iOM/24))]>2&&omIsClustered[iOM]){
	  h_clustered_highMult->Fill(fSignalCut_hotspot1, fSafetyWindow,1);
	  if (signalOMs_gen[iOM]) h_clusteredFrac_signal->Fill(fSignalCut_hotspot1, fSafetyWindow,1);
	}
      }
    }
  }
  fSafetyWindow=fSafetyWindow_buf;
  fSignalCut_hotspot1=fSignalCut_hotspot1_buf;
  return 1;
}

int BTimeClusters::sgn(float x)
{
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}


int BTimeClusters::findTriggerString(float high_threshold, float low_threshold)
{
  //find the string with trigger condition:
  //two adjacent modules have signals >high and >low in 100ns window
  int triggerString=-1;
  int n_impulses=fEvent->GetTotImpulses();

  //vector with trigger spot impulses;
  std::vector<int> trigger_impulses;
  //find pairs with trigger conditions
  for (int i=0; i<n_impulses; i++){
    int chID=fEvent->GetImpulse(i)->GetChannelID();
    if (fEvent->GetImpulse(i)->GetAmplitude()>=high_threshold){
      int secID=fGeomTel->At(chID)->GetSecNum();
      for (int j=0; j<n_impulses; j++){
	int chIDpair=fEvent->GetImpulse(j)->GetChannelID();
	int secIDpair=fGeomTel->At(chIDpair)->GetSecNum();
	if (fabs(chID-chIDpair)!=1||secID!=secIDpair) continue;
	if (fEvent->GetImpulse(i)->GetAmplitude()>=low_threshold){
	  trigger_impulses.push_back(i);
	  trigger_impulses.push_back(j);
	}
      }
    }
  }

  //find trigger impulse with max amplitude and pick corresponding string
  float ampl_max=0;
  int chIDmax=-1;
  for (int i=0; i<trigger_impulses.size(); i++){
    if (fEvent->GetImpulse(trigger_impulses[i])->GetAmplitude()>ampl_max){
      ampl_max=fEvent->GetImpulse(trigger_impulses[i])->GetAmplitude();
      chIDmax=fEvent->GetImpulse(trigger_impulses[i])->GetChannelID();
    }
  }

  int triggerSection=-1;
  if (chIDmax==-1) triggerString=6;
  else {
    triggerString=int(floor(chIDmax/24));
    triggerSection=fGeomTel->At(chIDmax)->GetSecNum();
  }
  
  return triggerString;
}


int BTimeClusters::countMCStringGaps(std::vector<int> string_impulses_gen)
{
  std::vector<bool> chanID;
  for (int i=0; i<192; i++) chanID.push_back(false);
  for (int i=0; i<string_impulses_gen.size(); i++){
    int idch=fMCEvent->GetHitChannel(string_impulses_gen[i])->GetChannelID()-1;
    idch=24*floor(idch/24)+(24-idch%24);
    idch=idch-1;
    chanID[idch]=true;
  }

  int nGaps=0;
  int nOMinGap=0;
  bool start=false;
  for (int i=0; i<chanID.size(); i++){
    if (chanID[i]) start=true;
    if (start&&(!chanID[i])) {
      nGaps++;
      nOMinGap++;
    }
    if (start&&chanID[i]) nOMinGap=0;
    if (start&&nOMinGap>=10) {
      nGaps=nGaps-nOMinGap;
      start=false;
    }
    if (start&&i==chanID.size()-1) {
      nGaps=nGaps-nOMinGap;
      start=false;
    }
  }
  //  std::cout<<"nGaps:   "<<nGaps<<"     nOMinGap: "<<nOMinGap<<std::endl;
  return nGaps;  
}



std::vector<std::vector<int>> BTimeClusters::findHotSpots(std::vector<int> string_impulses)
{
  //: PRODUCE ALL POSSIBLE CASUAL CONNECTED PAIRS OF IMPULSES AT NEIGHBOURING MODULES 
  //: ORDER PAIRS IN SUM OF IMPULSE AMPLITUDE
  
  std::vector<std::vector<int> > hot_spots;
  std::vector<bool> isClustered;
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++) isClustered.push_back(false);  
  
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
    if (fEvent->GetImpulse(string_impulses[iPulse])->GetAmplitude()<fSignalCut_hotspot1) continue;
    int chID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
    float timePulse=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
    for (int iCand=0; iCand<string_impulses.size(); iCand++){
      int id_increment=fEvent->GetImpulse(string_impulses[iCand])->GetChannelID()-chID;
      if (abs(id_increment)==1&&
	  fEvent->GetImpulse(string_impulses[iCand])->GetAmplitude()>fSignalCut_hotspot2){
	float timeCand=fEvent->GetImpulse(string_impulses[iCand])->GetTime();
	if (fabs(timeCand-timePulse)<(abs(id_increment)*15/cWater+fSafetyWindow)&&(!isClustered[iCand])){
	  std::vector<int> hotspot;
	  hotspot.push_back(iCand);
	  hotspot.push_back(iPulse);
	  hot_spots.push_back(hotspot);
	  isClustered[iPulse]=true;
	}
      }
    }
  }
  
  //ORDER PAIRS IN SUM IMPULSE
  std::vector<std::vector<int>> result;
  for (int i=0; i<hot_spots.size(); i++){
    std::vector<int> dummy;
    result.push_back(dummy);
  }
  
  for (int iHotspot=0; iHotspot<hot_spots.size(); iHotspot++){
    float sumAmpl=fEvent->GetImpulse(string_impulses[hot_spots[iHotspot][0]])->GetAmplitude()+
      fEvent->GetImpulse(string_impulses[hot_spots[iHotspot][1]])->GetAmplitude();
    int nLarger=0;
    for (int i=0; i<hot_spots.size(); i++){
      float probeSumAmpl=fEvent->GetImpulse(string_impulses[hot_spots[i][0]])->GetAmplitude()+
      fEvent->GetImpulse(string_impulses[hot_spots[i][1]])->GetAmplitude();
      if (sumAmpl<probeSumAmpl) nLarger++;
    }
    result[nLarger]=hot_spots[iHotspot];
  }

  return result;
 
  // return hot_spots;
}


std::pair<float,float> BTimeClusters::getPolarEstimate_string(std::vector<int> stringCluster)
{
  //work with fEventMask
  //estimate polar angle from all clustered hits
  //find cluster with largest amplitude
  //connect clusters which satisfy the reqirement:
  //deltaZ ~ deltaXY/tg(Theta)

  //dt=fEvent->GetImpulse(i)->GetTime();
  
  const double kCAngle = fGeomTel->GetCherenkovAngle();
  const double kCAngleCos = cos(kCAngle);
  const double kC = fGeomTel->GetVelocityVacuum();

  std::vector<float> thetaMin;
  std::vector<float> thetaMax;

  //  std::cout<<"string: "<<stringCluster.size()<<std::endl;
  
  for (int i = 0; i < stringCluster.size(); i++) {
    int chID=fEvent->GetImpulse(stringCluster[i])->GetChannelID();
    for (int j = i+1; j < stringCluster.size(); j++) {
      int chIDnext=fEvent->GetImpulse(stringCluster[j])->GetChannelID();
      if (abs(chID-chIDnext)==0) continue;
      float time=fEvent->GetImpulse(stringCluster[i])->GetTime();
      float timeNext=fEvent->GetImpulse(stringCluster[j])->GetTime();
      float dt = time-timeNext;
      // chPair.dt = dt; // for debugging
      float dz = fGeomTel->At(chID)->GetZ() - fGeomTel->At(chIDnext)->GetZ();
      //      std::cout<<" dt: "<<dt<<"   dz/kC "<<dz/kC<<std::endl;
      if (dz < 0) {
	dt = -dt;
	dz = -dz;
      }
      if (dt >= dz / kC && dt <= dz/(kCAngleCos * kC) + fSafetyWindow) {
	if (dt - fSafetyWindow < dz / kC) thetaMin.push_back(0);
	else thetaMin.push_back(kCAngle - acos(kC * (dt - fSafetyWindow) * kCAngleCos / dz));
	thetaMax.push_back(kCAngle + acos(kC * (dt -fSafetyWindow) * kCAngleCos / dz));
      }
      else if (dt < dz / kC && dt > - dz / kC) {
	if (dt + fSafetyWindow > dz / kC) thetaMin.push_back(0);
	else thetaMin.push_back(- kCAngle + acos(kC * (dt + fSafetyWindow) * kCAngleCos / dz));
	if (dt - fSafetyWindow < - dz / kC) thetaMax.push_back(M_PI);
	else thetaMax.push_back(kCAngle + acos(kC * (dt - fSafetyWindow) * kCAngleCos / dz));
      }
      else if (dt <= - dz / kC && dt >= - dz/(kCAngleCos * kC) - fSafetyWindow) {
	thetaMin.push_back(- kCAngle + acos(kC * (dt + fSafetyWindow) * kCAngleCos / dz));
	if (dt + fSafetyWindow > - dz / kC) thetaMax.push_back(M_PI);
	else thetaMax.push_back(2 * M_PI - kCAngle - acos(kC * (dt + fSafetyWindow) * kCAngleCos / dz));
      }
      else continue;

      
      
      /*  
      try {
	strings.at(fGeomTel->GetStringNum(ChID[i])).spair.push_back(chPair);
      } catch (const std::out_of_range & ) {
	OMString str;
	str.spair.push_back(chPair);
	strings.insert(std::pair<int, OMString>(fGeomTel->GetStringNum(ChID[i]), str));
      }
      strings.at(fGeomTel->GetStringNum(ChID[i])).channel.insert(std::pair<int, bool>(ChID[i], false));
      strings.at(fGeomTel->GetStringNum(ChID[j])).channel.insert(std::pair<int, bool>(ChID[j], false));
      */
    }
  }
  
  //find maximum thetaMin and minimum thetaMax
  
  float thetaMin_max=-1000;
  float thetaMax_min=1000;
  //  std::cout<<thetaMin.size()<<std::endl;
  for (int i=0; i<thetaMin.size(); i++){
    if (thetaMin[i]>thetaMin_max) thetaMin_max=thetaMin[i];
    if (thetaMax[i]<thetaMax_min) thetaMax_min=thetaMax[i];
  }

  return (std::make_pair(thetaMin_max, thetaMax_min));
    
}


std::pair<float,float> BTimeClusters::getClusterCenter(std::vector<int> stringCluster)
{
  float weightedSum=0;
  float amplSum=0;
  for (int i=0; i<stringCluster.size(); i++){
    float zPulse=fGeomTel->At(fEvent->GetImpulse(stringCluster[i])->GetChannelID())->GetZ();
    float amplPulse=fEvent->GetImpulse(stringCluster[i])->GetAmplitude();
    weightedSum+=zPulse*amplPulse;
    amplSum+=amplPulse;
  }

  if (weightedSum==0&&amplSum==0) return std::make_pair(0,0);
  float zCenter=weightedSum/amplSum;

  //initialise time with preliminary estimate from central pulse, then find closest to the center;
  float tCenter=fEvent->GetImpulse(stringCluster[int(floor(stringCluster.size()/2))])->GetTime();
  float zMin=1000;
  for (int i=0; i<stringCluster.size(); i++){
    float zPulse=fGeomTel->At(fEvent->GetImpulse(stringCluster[i])->GetChannelID())->GetZ();
    if (fabs(zPulse-zCenter)<zMin){
      zMin=fabs(zPulse-zCenter);
      tCenter=fEvent->GetImpulse(stringCluster[i])->GetTime();
    }
  }
  return std::make_pair(zCenter, tCenter);
}


std::vector<int> BTimeClusters::buildGlobalCluster(std::vector<std::vector<int> >* stringClusters)
{
  std::vector<int> globalCluster;
  
  //TRY TO LINK CLUSTERS ON DIFFERENT STRINGS:
  //FIND CLUSTER WITH LARGEST AMPLITUDE HOTSPOT (WITH LARGEST AMPLITUDE FOR THE MOMENT)
  float maxAmpl=0;
  int stringID=-1;
  int clusterID=-1;

  float minTheta=-1000;
  float maxTheta=1000;
  
  for (int iString=0; iString<8; iString++){
    for (int iCluster=0; iCluster<stringClusters[iString].size(); iCluster++){
      {
	//find min and max theta:
	if (stringClusters[iString][iCluster].size()>2){
	  std::pair<float,float> minmax=getPolarEstimate_string(stringClusters[iString][iCluster]);
	  if (minmax.first>minTheta) minTheta=minmax.first;
	  if (minmax.second<maxTheta) maxTheta=minmax.second;
	}
	float sum=0;
	for (int j=0; j<stringClusters[iString][iCluster].size(); j++){
	  sum+=fEvent->GetImpulse(stringClusters[iString][iCluster][j])->GetAmplitude();
	}
	if (sum>maxAmpl&&stringClusters[iString][iCluster].size()>2){
	  maxAmpl=sum;
	  stringID=iString;
	  clusterID=iCluster;
	}
      }
    }
  }

  //  float minTheta=1000;
  //  float maxTheta=-1000;
  
  if (stringID!=-1) {
    for (int i=0; i<stringClusters[stringID].size(); i++){
      globalCluster.push_back(stringClusters[stringID][clusterID][i]);
    }
    //  for (int iString=0; iString<8; iString++){
    //  if (stringClusters[iString].size()==0) continue;
    std::pair<float,float> center=getClusterCenter(stringClusters[stringID][clusterID]);
    /*
    if (usedStrings>1&&firedStrings>1&&stringID!=-1){
      std::cout<<"string: "<<stringID<<"     center z, t:   "<<center.first<<", "<<center.second<<std::endl;
      for (int k=0; k<stringClusters[stringID].size(); k++){
	std::cout<<"         hit: "<<k<<"   z, t: "<<
	  fGeomTel->At(fEvent->GetImpulse(stringClusters[stringID][k])->GetChannelID())->GetZ()<<", "<<
	  fEvent->GetImpulse(stringClusters[stringID][k])->GetTime()<<
	  std::endl;
      }
    }
    */
    //check how muon criterion works
    std::pair <float,float> minmax=getPolarEstimate_string(stringClusters[stringID][clusterID]);
    
    std::cout<<"seed cluster: "<<stringID<<"    size: "<<stringClusters[stringID][clusterID].size()<<"    "<<minTheta<<"   "<<maxTheta<<std::endl;
    
    //check other strings
    //link if deltaZ*deltaT<0
    //
    for (int iStringLink=0; iStringLink<8; iStringLink++){
      if (iStringLink==stringID) continue;
      for (int iClusterLink=0; iClusterLink<stringClusters[iStringLink].size(); iClusterLink++){
	std::pair<float,float> centerLink=getClusterCenter(stringClusters[iStringLink][iClusterLink]);
	if ((centerLink.first-center.first)*(center.second-centerLink.second)>0) continue;
	for (int j=0; j<stringClusters[iStringLink].size(); j++){
	  globalCluster.push_back(stringClusters[iStringLink][iClusterLink][j]);
	}
      }
    }
  }
  return globalCluster;
}
