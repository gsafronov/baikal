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
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"

#include <math.h>
#include <vector>

ClassImp(BTimeClusters);

BTimeClusters::BTimeClusters(const char *name, const char *title)
{
  //chanIDs=chans;
  // fMaskNoise=maskNoise;
  //fInputMaskName ="";
  // fOutputMaskName="";

  fName  = name ? name : "BTimeClusters";
  fTitle = title ? title : "TimeClusterFilterMask";

  cVacuum=0.3;
  cWater=cVacuum/1.33;

  fSignalCut_gen=2;
  fSignalCut_hotspot=1;
  fGen_rhoCut=20;
  fTimeMargin=50;
  
  //  fOutputMaskName = "AmplitudeFilterMask";

  fOUT=new TFile("timeClustersDebug.root","RECREATE");
  h_hits_per_string_2pe=new TH1F("h_hits_per_string_2pe","h_hits_per_string_2pe",25,0,25);
  h_hits_2pe=new TH1F("h_hits_2pe","h_hits_2pe",200,0,200);
  h_fired_strings_2pe=new TH1F("h_fired_strings_2pe","h_fired_strings_2pe",9,0,9);
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
  //if (fMCEvent->GetTrack(0)->GetMuonEnergy()<1000||
  if (rho>20) return kFALSE;
  /////////////////////////
  
  std::vector<int> string_impulses_gen[8];
  int noisePulseID=-100;
  int nHits=0;
  for (int i=0; i<n_impulse; i++) {
    //    std::cout<<i<<std::endl;
    BMCHitChannel* fHitChan=fMCEvent->GetHitChannel(i);
    int idch=fHitChan->GetChannelID()-1;
    idch=24*floor(idch/24)+(24-idch%24);
    idch=idch-1;
    int iString=floor(idch/24);
    float signal=0;
    for (int k=0; k<fHitChan->GetPulseN(); k++){
      if (fHitChan->GetPulse(k)->GetMagic()!=1) signal+=fHitChan->GetPulse(k)->GetAmplitude();
    }
    //    std::cout<<"impulse "<<i<<"   iString: "<<iString<<"   idch "<<idch<<"   bb  "<<fHitChan->GetChannelID()<<std::endl;
    if (signal>fSignalCut_gen) {
      string_impulses_gen[iString].push_back(i);
      nHits++;
    }
  }
  //  std::cout<<"popo"<<std::endl;
  h_hits_2pe->Fill(nHits,1);
  if (fMCEvent->GetResponseMuonsN()==1&&rho<fGen_rhoCut) h_1mu_hits_2pe->Fill(nHits,1);

  int firedStrings=0;
  for (int iString=0; iString<8; iString++){
    h_hits_per_string_2pe->Fill(string_impulses_gen[iString].size(),1);
    if (fMCEvent->GetResponseMuonsN()==1&&rho<fGen_rhoCut) h_1mu_hits_per_string_2pe->Fill(string_impulses_gen[iString].size(),1);
    if (string_impulses_gen[iString].size()>1) firedStrings++;
  }
  h_fired_strings_2pe->Fill(firedStrings,1);
  
  if (fMCEvent->GetResponseMuonsN()==1&&rho<fGen_rhoCut){
    // std::cout<<"NEXT EVENT     firedStrings: "<<firedStrings<<std::endl;
    h_1mu_fired_strings_2pe->Fill(firedStrings,1);
  }

  //  std::cout<<"zoob"<<std::endl;
  //detector-level impulses, combine in strings

  ////////////////////////////////////////////////////////////////ACTUAL RECO
  
  std::vector<int> string_impulses[8];

  //  std::cout<<"aa"<<std::endl;

  int nPulses_initial=0;
  int impulse_n=fEvent->GetTotImpulses();
  for (int iPulse=0; iPulse<impulse_n; iPulse++){
    if (fInputEventMask->GetFlag(iPulse)==0) continue;
    nPulses_initial++;
    int iChannel=fEvent->GetImpulse(iPulse)->GetChannelID();
    string_impulses[int(floor(iChannel/24))].push_back(iPulse);    
  }

  // std::cout<<"bb"<<std::endl;
  int hitStrings=0;
  std::vector<int> stringClusters[8];
  for (int iString=0; iString<8; iString++){
    if (string_impulses[iString].size()>0) hitStrings++;
    stringClusters[iString]=BuildStringCluster(iString, string_impulses[iString]);
  }

  //  std::cout<<"cc"<<std::endl;
  
  //build global cluster [SIMPLE]
  int usedStrings=0;
  
  std::vector<int> globalCluster;
  for (int iString=0; iString<8; iString++){
    if (stringClusters[iString].size()>0) usedStrings++;
    if (fMCEvent->GetResponseMuonsN()==1&&rho<fGen_rhoCut) hreco_hits_per_string->Fill(stringClusters[iString].size(),1);
    for (int iPulse=0; iPulse<stringClusters[iString].size(); iPulse++){
      globalCluster.push_back(stringClusters[iString][iPulse]);
    }
    h_hitsPerString_reco_vs_gen->Fill(string_impulses_gen[iString].size(),stringClusters[iString].size(),1);
  }

  if (fMCEvent->GetResponseMuonsN()==1&&rho<fGen_rhoCut) hreco_fired_strings->Fill(usedStrings,1);

  //set the event mask
  for (int iPulse=0; iPulse<impulse_n; iPulse++){
    fOutputEventMask->SetFlag(iPulse,0);
  }

  for (int iClustered=0; iClustered<globalCluster.size(); iClustered++){
    fOutputEventMask->SetFlag(globalCluster[iClustered],1);
  }

  if (fMCEvent->GetResponseMuonsN()==1&&rho<fGen_rhoCut) {
    h_strings_reco_vs_gen->Fill(firedStrings, usedStrings,1);
    h_strings_bevt_vs_gen->Fill(firedStrings, hitStrings,1);
    hreco_hits->Fill(globalCluster.size(),1);
    //std::cout<<"OLD METH STRINGS: "<<usedStrings<<std::endl;
  }
  //  std::cout<<"initial pulses: "<<nPulses_initial<<"   after clusterisation: "<<globalCluster.size()<<" used strings:   "<<usedStrings<<"     mc response muons:   "<<fMCEvent->GetResponseMuonsN()<<std::endl;
  
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
    int chID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
    float timePulse=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
    for (int iCand=iPulse+1; iCand<string_impulses.size(); iCand++){
      //if (fEvent->GetImpulse(string_impulses[iCand])->GetAmplitude()<fSignalCut_hotspot||
      if (isClustered[iCand]) continue;
      if (fabs(chID-fEvent->GetImpulse(string_impulses[iCand])->GetChannelID())!=1) continue;
      float timeCand=fEvent->GetImpulse(string_impulses[iCand])->GetTime();
      // std::cout<<timeCand<<"      "<<timePulse<<"       "<<fabs(timePulse-timeCand)<<std::endl;
      if (fabs(timeCand-timePulse)<(15/cWater)+fTimeMargin){
	isClustered[iCand]=true;
	hotspot.push_back(iCand);
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
  
  for (int i=0; i<hot_spots.size(); i++){
    //find max and next to max hits from the hotspot
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
    
    //std::cout<<max_id<<"  "<<submax_id<<std::endl;
    //they should be on adjacent channels
    if (fabs(fEvent->GetImpulse(string_impulses[max_id])->GetChannelID()-
	     fEvent->GetImpulse(string_impulses[submax_id])->GetChannelID())!=1) continue;
    //if ok, add further hits to the hotspot
    std::vector<int> stringCluster=AddClusterImpulses(max_id, submax_id, string_impulses, fTimeMargin);
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

std::vector<int> BTimeClusters::AddClusterImpulses(int max_ampl_id, int adjacent_id, std::vector<int> string_impulses, float window)
{
  std::vector<int> stringCluster;
  std::vector<int> storeys_pulse_id;
  std::vector<float> storeys_pulse_time;
  
  for (int i=0; i<24; i++){
    storeys_pulse_id.push_back(-1);
    storeys_pulse_time.push_back(-1);
  }
  
  stringCluster.push_back(string_impulses[max_ampl_id]);
  stringCluster.push_back(string_impulses[adjacent_id]);
  fOutputEventMask->SetFlag(string_impulses[max_ampl_id],10);
  fOutputEventMask->SetFlag(string_impulses[adjacent_id],10);
  
  //  std::cout<<"bbbbb3"<<std::endl;
  
  float seed_time=fEvent->GetImpulse(string_impulses[max_ampl_id])->GetTime();
  int seed_channel_id=fEvent->GetImpulse(string_impulses[max_ampl_id])->GetChannelID();
  int seed_storey=(seed_channel_id)%24;
  float adjacent_time=fEvent->GetImpulse(string_impulses[adjacent_id])->GetTime();
  int adjacent_channel_id=fEvent->GetImpulse(string_impulses[adjacent_id])->GetChannelID();

  float deltaZ=fGeomTel->At(adjacent_channel_id)->GetZ()-fGeomTel->At(seed_channel_id)->GetZ();
  float deltaT=adjacent_time-seed_time;
  int id_increment=adjacent_channel_id-seed_channel_id;

  //  std::cout<<"bbbbb4"<<std::endl;
  
  //add +-1 channels
  float time_early=min(seed_time-(deltaT/fabs(deltaT))*fabs(deltaZ)/cWater-window, seed_time-deltaT-window); //experimental
  float time_late=max(seed_time-(deltaT/fabs(deltaT))*fabs(deltaZ)/cWater+window, seed_time-deltaT+window);
  
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
    if (iPulse==max_ampl_id||iPulse==adjacent_id) continue;
    int chanID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
    if (chanID==seed_channel_id-id_increment){
      float chan_time=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
      if (chan_time>time_early&&chan_time<time_late){
	stringCluster.push_back(string_impulses[iPulse]);
	fOutputEventMask->SetFlag(string_impulses[iPulse],10);
      }
    }
  }

  if (stringCluster.size()==2) {
    //std::cout<<"size 2"<<std::endl;
    return stringCluster;
  }

  //add up to 7 channels
  //fill storeys vector
  for (int i=0; i<stringCluster.size(); i++){
    storeys_pulse_id[(fEvent->GetImpulse(stringCluster[i])->GetChannelID())%24]=stringCluster[i];
    storeys_pulse_time[(fEvent->GetImpulse(stringCluster[i])->GetChannelID())%24]=fEvent->GetImpulse(stringCluster[i])->GetTime();
  }

  for (int i=-1; i<3; i+=2){
    for (int j=1; j<3; j++){
      //     std::cout<<"kkk"<<std::endl;
      //find min and max time for new channel
      if (seed_storey+i+i*j<0||seed_storey+i+i*j>23) continue;
      float tim0=storeys_pulse_time[seed_storey+i]+j*(storeys_pulse_time[seed_storey+i]-storeys_pulse_time[seed_storey]);
      float tim1=storeys_pulse_time[seed_storey+i]+j*(storeys_pulse_time[seed_storey+i]-storeys_pulse_time[seed_storey-i])/2;
      float tim2=storeys_pulse_time[seed_storey]+(j+1)*(storeys_pulse_time[seed_storey]-storeys_pulse_time[seed_storey-i]);

      float min_estimate=min(tim0,tim1);
      min_estimate=min(min_estimate,tim2);

      float max_estimate=max(tim0,tim1);
      max_estimate=max(max_estimate,tim2);

      //loop over string pulses
      for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
	if (fOutputEventMask->GetFlag(string_impulses[iPulse])==10) continue;
	int chanID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
	if (chanID!=seed_channel_id+i+i*j) continue;
	float time=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
	if (time<max_estimate+fTimeMargin&&time>min_estimate-fTimeMargin){
	  stringCluster.push_back(string_impulses[iPulse]);
	  storeys_pulse_id[seed_storey+i+i*j]=string_impulses[iPulse];
	  storeys_pulse_time[seed_storey+i+i*j]=string_impulses[iPulse];
	  fOutputEventMask->SetFlag(string_impulses[iPulse],10);
	} 
      }
    }
  }
  //  std::cout<<stringCluster.size()<<std::endl;
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
    
  fOUT->cd();
  h_hits_per_string_2pe->Write();
  h_hits_2pe->Write();
  h_fired_strings_2pe->Write();
  h_1mu_hits_per_string_2pe->Write();
  h_1mu_hits_2pe->Write();
  h_1mu_fired_strings_2pe->Write();
  hreco_hits_per_string->Write();
  hreco_fired_strings->Write();
  hreco_hits->Write();
  h_strings_reco_vs_gen->Write();
  h_strings_bevt_vs_gen->Write();

  h_hitsPerString_reco_vs_gen->Write();
  h_hitsPerString_bevt_vs_gen->Write();
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
