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
  fSignalCut_hotspot=1.5;
  fTimeMargin=80;

  //gen-level cuts
  fGen_rhoCut=40;
  fGen_minAngle=45;
  fGen_maxAngle=90;

  fSignalCut_gen=1.0;
  fSeedSignalCut_gen=4;

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

  h_noiseFrac_clustered=new TH2F("h_noiseFrac_clustered","h_noiseFrac_clustered",20,0.5,10.5,20,0,400);
  h_signalFrac_clustered=new TH2F("h_signalFrac_clustered","h_signalFrac_clustered",20,0.5,10.5,20,0,400);
  h_stringSignalFrac_clustered=new TH2F("h_stringSignalFrac_clustered","h_stringSignalFrac_clustered",20,0.5,10.5,20,0,400);
  
  h_clusteredFrac_noise=new TH2F("h_clusteredFrac_noise","h_clusteredFrac_noise",20,0.5,10.5,20,0,400);
  h_clusteredFrac_signal=new TH2F("h_clusteredFrac_signal","h_clusteredFrac_signal",20,0.5,10.5,20,0,400);

  h_effpur_mult=new TH2F("h_effpur_mult","h_effpur_mult",20,0.5,10.5,20,0,400);
  
  h_clustered=new TH2F("h_clustered","h_clustered",20,0.5,10.5,20,0,400);
    
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
  int n_impulse=fMCEvent->GetChannelN();
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
  
  if (fMCEvent->GetResponseMuonsN()!=1) return kFALSE;
  if (rho>fGen_rhoCut) return kFALSE;
  if (180-fMCEvent->GetPrimaryParticlePolar()>fGen_maxAngle) return kFALSE;
  if (180-fMCEvent->GetPrimaryParticlePolar()<fGen_minAngle) return kFALSE;
  //  if (

  //  if (fMCEvent->GetTrack(0)->GetMuonEnergy()<1000||
  //  if (rho>20) return kFALSE;
  /////////////////////////
  
  std::vector<int> string_impulses_gen[8];
  
  int noisePulseID=-100;
  int nHits=0;
  bool hasLargeHit=false;

  float maxAmplitude_gen=-1;
  int maxAmplitude_gen_chanID=-1;

  std::vector<float> singlePulse_ampls;
  std::vector<int> singlePulse_ids;

  std::vector<bool> noiseOMs_gen;
  std::vector<bool> signalOMs_gen;
  std::vector<bool> signalStrings_gen;

  for (int i=0; i<192; i++) {
    noiseOMs_gen.push_back(false);
    signalOMs_gen.push_back(false);
  }

  for (int i=0; i<8; i++){
    signalStrings_gen.push_back(false);
  }    
  
  for (int i=0; i<n_impulse; i++) {
    //    std::cout<<i<<std::endl;
    BMCHitChannel* fHitChan=fMCEvent->GetHitChannel(i);
    int idch=fHitChan->GetChannelID()-1;
    idch=24*floor(idch/24)+(24-idch%24);
    idch=idch-1;
    int iString=floor(idch/24);
    hitMap_gen->Fill(iString, idch%24,1);
    float signal=0;
    bool isNoise=false;
    bool isSignal=false;
    for (int k=0; k<fHitChan->GetPulseN(); k++){
      if (fHitChan->GetPulse(k)->GetMagic()!=1) {
	signal+=fHitChan->GetPulse(k)->GetAmplitude();
	isSignal=true;
      }
      if (fHitChan->GetPulse(k)->GetMagic()==1) isNoise=true;
    }

    if ((isNoise&&!isSignal)) {
      noiseOMs_gen[idch]=true;
      fNoiseOMs++;
    }
    if (isSignal) {
      signalOMs_gen[idch]=true;
      fSignalOMs++;
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
    if (signal>fSignalCut_gen) {
      string_impulses_gen[iString].push_back(i);
      nHits++;
    }
    
    if (signal>fSeedSignalCut_gen) hasLargeHit=true;

    if (signal>maxAmplitude_gen){
      maxAmplitude_gen=signal;
      maxAmplitude_gen_chanID=idch;
    }
  }

  if (!hasLargeHit) return kFALSE;
  
  h_muon_energy->Fill(fMCEvent->GetTrack(0)->GetMuonEnergy(),1);
  if (nHits>=5&&hasLargeHit) h_muon_energy_rcand_gen->Fill(fMCEvent->GetTrack(0)->GetMuonEnergy(),1);
  
  if (nHits>=5) h_rcand_polar_vs_rho->Fill(180-fMCEvent->GetPrimaryParticlePolar(),rho,1);
  h_smuons_polar_vs_rho->Fill(180-fMCEvent->GetPrimaryParticlePolar(),rho,1);
  
  //  std::cout<<"popo"<<std::endl;
  h_1mu_hits_2pe->Fill(nHits,1);

  
  int firedStrings=0;
  for (int iString=0; iString<8; iString++){
    h_1mu_hits_per_string_2pe->Fill(string_impulses_gen[iString].size(),1);
    if (string_impulses_gen[iString].size()>1) {
      fSignalStrings++;
      firedStrings++;
      signalStrings_gen[iString]=true;
    }
    //    if string_impulses_gen[iString].size()std::cout<<"gen string:  "<<iString<<":   "<<string_impulses_gen[iString].size()<<" hits"<<std::endl;
       
  }
  
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
  int impulse_n=fEvent->GetTotImpulses();
  for (int iPulse=0; iPulse<impulse_n; iPulse++){
    if (fInputEventMask->GetFlag(iPulse)==0) continue;
    nPulses_initial++;
    int iChannel=fEvent->GetImpulse(iPulse)->GetChannelID();

    hitMap_det->Fill(int(floor(iChannel/24)), iChannel%24,1);
    
    string_impulses[int(floor(iChannel/24))].push_back(iPulse);
    
    for (int iGen=0; iGen<singlePulse_ids.size(); iGen++){
      if (iChannel==singlePulse_ids[iGen]) h_hitSignal_reco_vs_gen->Fill(singlePulse_ampls[iGen], fEvent->GetImpulse(iPulse)->GetAmplitude(),1); 
    }
    if (fEvent->GetImpulse(iPulse)->GetAmplitude()>maxAmplitude_rec){
      maxAmplitude_rec=fEvent->GetImpulse(iPulse)->GetAmplitude();
      maxAmplitude_rec_chanID=iChannel;
    }
    
  }

  // if (maxAmplitude_rec_chanID==maxAmplitude_gen_chanID) h_hitSignal_reco_vs_gen->Fill(maxAmplitude_gen,  maxAmplitude_rec,1);
  //else h_hitSignal_reco_vs_gen->Fill(-1,-1,1);
  
  //std::cout<<"bb"<<std::endl;
  std::vector<int> stringClusters[8];
  int hitStrings=0;
  if (fVerbose) std::cout<<">>>>>>>>>>>>>>>>>>>>>>>>>NEXT EVENT"<<std::endl;
  for (int iString=0; iString<8; iString++){
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
	stringClusters[iString]=BuildStringCluster(iString, string_impulses[iString]);
      }
    }
    else stringClusters[iString]=BuildStringCluster(iString, string_impulses[iString]);
  }

  //  std::cout<<"cc"<<std::endl;
  
  //build global cluster [SIMPLE]
  int usedStrings=0;
  
  std::vector<int> globalCluster;
  for (int iString=0; iString<8; iString++){
    if (stringClusters[iString].size()>0) usedStrings++;
    hreco_hits_per_string->Fill(stringClusters[iString].size(),1);
    for (int iPulse=0; iPulse<stringClusters[iString].size(); iPulse++){
      globalCluster.push_back(stringClusters[iString][iPulse]);
    }
    h_hitsPerString_reco_vs_gen->Fill(string_impulses_gen[iString].size(),stringClusters[iString].size(),1);
  }
  if (fVerbose) std::cout<<"firedStrings: "<<firedStrings<<"     recoStrings: "<<usedStrings<<std::endl;
  if (usedStrings>=2&&globalCluster.size()>=5) h_muon_energy_rcand->Fill(fMCEvent->GetTrack(0)->GetMuonEnergy(),1);
  
  hreco_fired_strings->Fill(usedStrings,1);

  //set the event mask
  for (int iPulse=0; iPulse<impulse_n; iPulse++){
    fOutputEventMask->SetFlag(iPulse,0);
  }

  for (int iClustered=0; iClustered<globalCluster.size(); iClustered++){
    fOutputEventMask->SetFlag(globalCluster[iClustered],1);
  }

  h_strings_reco_vs_gen->Fill(firedStrings, usedStrings,1);
  h_strings_bevt_vs_gen->Fill(firedStrings, hitStrings,1);
  hreco_hits->Fill(globalCluster.size(),1);
    //std::cout<<"OLD METH STRINGS: "<<usedStrings<<std::endl;
  //  std::cout<<"initial pulses: "<<nPulses_initial<<"   after clusterisation: "<<globalCluster.size()<<" used strings:   "<<usedStrings<<"     mc response muons:   "<<fMCEvent->GetResponseMuonsN()<<std::endl;

  runWPscan(string_impulses, noiseOMs_gen, signalOMs_gen, signalStrings_gen);
  
  return kTRUE;
}

std::vector<int> BTimeClusters::BuildStringCluster(int iString, std::vector<int> string_impulses)
{
  // std::cout<<"BUILD"<<std::endl;
  std::vector<int> stringCluster;

  if (string_impulses.size()==0) return stringCluster;
  
  ///////////////////TEMPORARY
  float primThetaRad=M_PI*fMCEvent->GetPrimaryParticlePolar()/180;
  float primPhiRad=M_PI*fMCEvent->GetPrimaryParticleAzimuth()/180;

  TVector3 genVec(sin(primThetaRad)*cos(primPhiRad),
		  sin(primThetaRad)*sin(primPhiRad),
		  cos(primThetaRad));

  TVector3 inPoTrack(fMCEvent->GetTrack(0)->GetX()-500*genVec.X(), fMCEvent->GetTrack(0)->GetY()-500*genVec.Y(), fMCEvent->GetTrack(0)->GetZ()-500*genVec.Z());

  TVector3 zero(0,0,0);
  float rho=getTrackDistanceToOM(inPoTrack, genVec, zero);
  
  
  ///////////////?TEST: PRODUCE ALL POSSIBLE CASUAL CONNECTED PAIRS OF IMPULSES AT NEIGHBOURING MODULES 
  std::vector<std::vector<int>> hot_spots;
  std::vector<bool> isClustered;
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++) isClustered.push_back(false);  
  
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
    if (fEvent->GetImpulse(string_impulses[iPulse])->GetAmplitude()<fSignalCut_hotspot||isClustered[iPulse]) continue;
    std::vector<int> hotspot;
    hotspot.push_back(iPulse);
    isClustered[iPulse]=true;
    int chID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
    float timePulse=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
    for (int iCand=iPulse+1; iCand<string_impulses.size(); iCand++){
      //if (fEvent->GetImpulse(string_impulses[iCand])->GetAmplitude()<fSignalCut_hotspot||
      if (isClustered[iCand]) continue;
      if (fabs(chID-fEvent->GetImpulse(string_impulses[iCand])->GetChannelID())!=1) continue;
      float timeCand=fEvent->GetImpulse(string_impulses[iCand])->GetTime();
      //      std::cout<<timeCand<<"      "<<timePulse<<"       "<<fabs(timePulse-timeCand)<<std::endl;
      if (fabs(timeCand-timePulse)<(15/cWater)+30){
	isClustered[iCand]=true;
	hotspot.push_back(iCand);
	break;
      }
    }
    if (hotspot.size()>1) hot_spots.push_back(hotspot);
  }

  
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
    
    //find max and next to max hits from the hotspot [WHAT FOR??]
    /*
    int max_id=-1;
    int submax_id=-1;
    for (int j=0; j<hot_spots[i].size(); j++){
      int nLargerAmpl=0;
      for (int k=0; k<hot_spots[i].size(); k++){
	if(k==j) continue;     
	if (fEvent->GetImpulse(string_impulses[hot_spots[i][j]])->GetAmplitude()<
	    fEvent->GetImpulse(string_impulses[hot_spots[i][k]])->GetAmplitude()) nLargerAmpl++;
      }
      if (nLargerAmpl==0) max_id=hot_spots[i][j];
      if (nLargerAmpl==1) submax_id=hot_spots[i][j];
    }
    */
    //std::cout<<max_id<<"  "<<submax_id<<std::endl;
    //they should be on adjacent channels (??????)
    //    if (fabs(fEvent->GetImpulse(string_impulses[max_id])->GetChannelID()-
    //	     fEvent->GetImpulse(string_impulses[submax_id])->GetChannelID())!=1) continue;
    //if ok, add further hits to the hotspot
    std::vector<int> stringCluster=AddClusterImpulses(hot_spots[i], string_impulses, fTimeMargin);
    //std::cout<<"check"<<std::endl;
    //check if the cluster is the biggest in the event
    if (stringCluster.size()>max_cluster_size){
      max_cluster_size=stringCluster.size();
      max_size_cluster=stringCluster;
    }
  }
  
  /*    
  ////////////
  //find seed channel - the one with maximum signal
  
  float max_ampl=-1;
  int max_ampl_id=-1;
  
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
    if (max_ampl<fEvent->GetImpulse(string_impulses[iPulse])->GetAmplitude()){
      max_ampl=fEvent->GetImpulse(string_impulses[iPulse])->GetAmplitude();
      max_ampl_id=iPulse;
    }
  }

  if (string_impulses.size()==0) return stringCluster;
  
  // std::cout<<"bbbbb1:  "<<string_impulses.size()<<"  "<<max_ampl<<"  "<<max_ampl_id<<std::endl;
  
  //check signal in seed neighbours
  float adjacent_ampl=0;
  int adjacent_id=-1;
  int seed_channel_id=fEvent->GetImpulse(string_impulses[max_ampl_id])->GetChannelID();
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
    int channel_id=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
    if (!(channel_id==seed_channel_id+1||channel_id==seed_channel_id-1)) continue;
    if (adjacent_ampl<fEvent->GetImpulse(string_impulses[iPulse])->GetAmplitude()){
      adjacent_ampl=fEvent->GetImpulse(string_impulses[iPulse])->GetAmplitude();
      adjacent_id=iPulse;
    }
  }

  //fail if no adjacent channel with ampl > 2 pe. But need to check other seeds then. Possibly build pairs of adjacent channels with A > 2pe and treat the one with higher signal as seed.
  
  //  std::cout<<"bbbbb2"<<std::endl;

  if (adjacent_ampl<fSignalCut_hotspot) return stringCluster;

  //hot spot is found, run cluster buiding procedure
  float window = fTimeMargin;
  stringCluster=AddClusterImpulses(max_ampl_id, adjacent_id, string_impulses, window);
  */
  //  std::cout<<"bloblo"<<std::endl;
  return max_size_cluster;
  //return stringCluster;
}

std::vector<int> BTimeClusters::AddClusterImpulses(std::vector<int> hotspot, std::vector<int> string_impulses, float window)
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
  //  std::cout<<"bbbbb3"<<std::endl;

  for (int iSeed=0; iSeed<hotspot.size(); iSeed++){
    
    float seed_time=fEvent->GetImpulse(string_impulses[hotspot[iSeed]])->GetTime();
    int seed_channel_id=fEvent->GetImpulse(string_impulses[hotspot[iSeed]])->GetChannelID();
    int seed_storey=(seed_channel_id)%24;
    float adjacent_time=fEvent->GetImpulse(string_impulses[hotspot[1-iSeed]])->GetTime();
    int adjacent_channel_id=fEvent->GetImpulse(string_impulses[hotspot[1-iSeed]])->GetChannelID();
    
    float deltaZ=fGeomTel->At(seed_channel_id)->GetZ()-fGeomTel->At(adjacent_channel_id)->GetZ();
    float deltaT=seed_time-adjacent_time;
    int id_increment=seed_channel_id-adjacent_channel_id;

    //  std::cout<<"bbbbb4"<<std::endl;

    //add +-1 channels
    float time_early=min(seed_time+sgn(deltaT)*fabs(deltaZ)/cWater-window, seed_time+deltaT-window); //experimental
    float time_late=max(seed_time+sgn(deltaT)*fabs(deltaZ)/cWater+window, seed_time+deltaT+window);

    if (fVerbose) std::cout<<" deltaT: "<<deltaT<<"   seed_time: "<<seed_time<<"    seed chan: "<<seed_channel_id<<"   t_early: "<<time_early<<"   t_late: "<<time_late<<std::endl;
    
    for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
      if (iPulse==hotspot[iSeed]||iPulse==hotspot[1-iSeed]) continue;
      int chanID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
      float chan_time=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
      if (fVerbose) std::cout<<"time: "<<chan_time<<"       chanID: "<<chanID<<std::endl;
      if (chanID==seed_channel_id+id_increment){
	if (chan_time>time_early&&chan_time<time_late){
	  stringCluster.push_back(string_impulses[iPulse]);
	  fOutputEventMask->SetFlag(string_impulses[iPulse],10);
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
    /*
      if (stringCluster.size()==2) {
      //std::cout<<"size 2"<<std::endl;
      return stringCluster;
      }
    */
    
    //add up to 7 channels
    
    int id[]={1,3};
    for (int k=0; k<2; k++){
      int i=id_increment*id[k];
      if (storeys_pulse_id[seed_storey+i]==-1) continue;
      for (int j=1; j<3; j++){
	//     std::cout<<"kkk"<<std::endl;
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
	  if (fOutputEventMask->GetFlag(string_impulses[iPulse])==10) continue;
	  int chanID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
	  if (chanID!=seed_channel_id+i+sgn(i)*j) continue;
	  float time=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
	  if (time<max_estimate+fTimeMargin&&time>min_estimate-fTimeMargin){
	    stringCluster.push_back(string_impulses[iPulse]);
	    storeys_pulse_id[seed_storey+i+sgn(i)*j]=string_impulses[iPulse];
	    storeys_pulse_time[seed_storey+i+sgn(i)*j]=string_impulses[iPulse];
	    fOutputEventMask->SetFlag(string_impulses[iPulse],10);
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
    double normFactor=h_1mu_hits_per_string_2pe->GetBinContent(i+1);
    if (normFactor==0) continue;
    for (int j=0; j<25; j++){
      double bico=(double)h_hitsPerString_reco_vs_gen->GetBinContent(i+1,j+1)/h_1mu_hits_per_string_2pe->GetBinContent(i+1);
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
  h_clusteredFrac_signal->Divide(h_clusteredFrac_signal,h_clustered,1,1);
  
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
  std::vector<int> stringClusters[8];
  int hitStrings=0;
  for (int i=0; i<20; i++){
    fSignalCut_hotspot=0.5+0.5*i;
    for (int j=0; j<20; j++){
      fTimeMargin=20*j;
      std::vector<int> clusteredOMs;
      for (int iString=0; iString<8; iString++){
	if (string_impulses[iString].size()>0) hitStrings++;
	stringClusters[iString]=BuildStringCluster(iString, string_impulses[iString]);

	if (stringClusters[iString].size()>0&&signalStrings_gen[iString]) h_stringSignalFrac_clustered->Fill(fSignalCut_hotspot, fTimeMargin,1);

	for (int k=0; k<stringClusters[iString].size(); k++){
	  int omID=fEvent->GetImpulse(stringClusters[iString][k])->GetChannelID();
	  clusteredOMs.push_back(omID);
	}
      }
      for (int iOM=0; iOM<clusteredOMs.size(); iOM++){
	h_clustered->Fill(fSignalCut_hotspot, fTimeMargin,1);
	if (noiseOMs_gen[clusteredOMs[iOM]]) h_clusteredFrac_noise->Fill(fSignalCut_hotspot, fTimeMargin,1);
	if (signalOMs_gen[clusteredOMs[iOM]]) h_clusteredFrac_signal->Fill(fSignalCut_hotspot, fTimeMargin,1);
      }
    }
  }
  
  return 1;
}

int BTimeClusters::sgn(float x)
{
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}
