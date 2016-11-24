#include "BClusterProducer.h"
#include "BGeomTel.h"
#include "BEvent.h"
#include "BMCEvent.h"
#include "BExtractedImpulse.h"
#include "BExtractedImpulseTel.h"
#include "BJoinExtractedImpulseTel.h"
#include "BEventMask.h"
#include "BChannelMask.h"
#include "BExtractedHeader.h"

//#include "BStringCluster.h"

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

ClassImp(BClusterProducer);

BClusterProducer::BClusterProducer(const char *name, const char *title)
{
  // chanIDs=chans;
  // fMaskNoise=maskNoise;
  // fInputMaskName ="";
  // fOutputMaskName="";

  fName  = name ? name : "BClusterProducer";
  fTitle = title ? title : "TimeClusterFilterMask";

  cVacuum=0.3;
  cWater=cVacuum/1.33;

  //define WP
  fSignalCut_hotspot1=2.5;
  fSignalCut_hotspot2=-1.5;
  fSafetyWindow=10;

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

  fOUT=new TFile("timeClusterProducerDebug.root","RECREATE");

  h_ntracks=new TH1F("h_ntracks","h_ntracks",100,0,100);
  h_ntracks_6hits3strings=new TH1F("h_ntracks_6hits3strings","h_ntracks_6hits3strings",100,0,100);
  h_debug=new TH1F("h_debug","h_debug",100, 0,1);
  
  h_strClu_ntracks=new TH1F("h_strClu_ntracks","h_strClu_ntracks",20,0,20);
  h_gloClu_ntracks=new TH1F("h_gloClu_ntracks","h_gloClu_ntracks",20,0,20);
  h_gloClu_vs_firedStrings=new TH2F("h_gloClu_vs_firedStrings","h_gloClu_vs_firedStrings", 10,0,10,10,0,10);
				    
  h_pulse_ntracks=new TH1F("h_pulse_ntracks","h_pulse_ntracks",20,0,20);

  h_strClu_noiseFrac=new TH1F("h_strClu_noiseFrac","h_strClu_noiseFrac",100,0,1);
}



BClusterProducer::~BClusterProducer()
{
}

Int_t BClusterProducer::PreProcess(MParList * pList)
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


Bool_t BClusterProducer::Filter()
{
  fEventCounter++;
  if (fEventCounter%10000==0) std::cout<<fEventCounter<<std::endl;
  int n_impulse_channels=fMCEvent->GetChannelN();

  h_ntracks->Fill(fMCEvent->GetResponseMuonsN(),1);
  if (fMCEvent->GetResponseMuonsN()<1) return kFALSE;

  //check multiplicity in event with 6 hits at 3 strings
  int firedString[8]={0};
  int signalCh[194]={0};
  int nHitChans=fMCEvent->GetChannelN();
  int nHits=0;
  int nStrings=0;
  for (int i=0; i<nHitChans; i++){
    BMCHitChannel* ch=fMCEvent->GetHitChannel(i);
    int idch=ch->GetChannelID()-1;
    idch=24*floor(idch/24)+(24-idch%24);
    idch=idch-1;
    int iString=floor(idch/24);
    bool isSignal=false;
    float signal=0;
    for (int j=0; j<ch->GetPulseN(); j++){
      if (ch->GetPulse(j)->GetMagic()!=1) {
	signal+=ch->GetPulse(j)->GetAmplitude();
      }
    }
    signalCh[idch+1]=signal;
  }

  for (int i=1; i<194; i++){
    if (signalCh[i+1]>2.5&&
	((signalCh[i]>0&&int(floor(i-1)/24)==int(floor(i)/24))
	 ||(signalCh[i+2]>0.1&&int(floor(i)/24)==int(floor(i+1)/24)))) firedString[int(floor(i-1)/24)]++ ;
    if (signalCh[i+1]>0) nHits++;
  }

  for (int i=0; i<8; i++){
    if (firedString[i]>0) nStrings++;
  }
  
  if (nStrings>2&&nHits>5) h_ntracks_6hits3strings->Fill(fMCEvent->GetResponseMuonsN(),1);
  //////////////////////////////////////////////////////////////
  
  std::vector<int> string_impulses[8];
  int impulse_n=fEvent->GetTotImpulses();
  
  for (int iPulse=0; iPulse<impulse_n; iPulse++){
    if (fInputEventMask->GetFlag(iPulse)==0) continue;
    int iChannel=fEvent->GetImpulse(iPulse)->GetChannelID();
    //exclude trigger string:
    int iString=int(floor(iChannel/24));
    string_impulses[iString].push_back(iPulse);
   }

  //build string clusters
  std::vector<BStringCluster> stringClusters[8];
  for (int iString=0; iString<8; iString++){
    stringClusters[iString]=buildStringClusters(iString, string_impulses[iString]);
    //plot number of tracks in max hotspot ampl cluster
    float max=-1;
    BStringCluster leadingCluster;
    for (int k=0; k<stringClusters[iString].size(); k++){
      if (stringClusters[iString][k].GetSize()>4) h_strClu_ntracks->Fill(stringClusters[iString][k].GetNSignalTracks(),1);
      if (max<stringClusters[iString][k].GetHotSpotAmpl()){
	max=stringClusters[iString][k].GetHotSpotAmpl();
	leadingCluster=stringClusters[iString][k];
      }
      for (int iPu=0; iPu<stringClusters[iString][k].GetSize(); iPu++){
	h_pulse_ntracks->Fill(stringClusters[iString][k].GetNSignalTracksHit(iPu),1);
      }
      //      h_strClu_ntracks->Fill(stringClusters[iString][k].GetNSignalTracks(),1);
    }
    //if (leadingCluster.GetHotSpotAmpl()/leadingCluster.GetSumAmpl() > 0.8)
    // if (leadingCluster.GetSize()>3) h_strClu_ntracks->Fill(leadingCluster.GetNSignalTracks(),1);
    //    h_strClu_noiseFrac->Fill(
  }
  
  //cout string clusters
  int nStr=0;
  for (int i=0; i<8; i++){
    //  std::cout<<"size of string clusters at string "<<i<<": "<<stringClusters[i].size()<<std::endl;
    if (stringClusters[i].size()>0) nStr++;
    for (int j=0; j<stringClusters[i].size(); j++){
      // std::cout<<"       z of cluster "<<j<<":  "<<stringClusters[i][j].GetCenterZ()<<std::endl;
    }
  }
  
  //  if (!(fMCEvent->GetResponseMuonsN()==1&&nStrings>2&&nHits>5&&firedString[7]>1)) return kFALSE;

  //  std::cout<<">>>>>>>>>>>>>NEXT EVENT"<<std::endl;
  
  std::vector<BStringCluster> globalCluster=buildGlobalCluster(stringClusters);

  int nTracksInGlobalCluster=0;
  int trackCounter[100]={0};
  for (int i=0; i<globalCluster.size(); i++){
    std::vector<int> tracks=globalCluster[i].GetTracks();
    for (int j=0; j<tracks.size(); j++) trackCounter[tracks[j]]++;

    for (int j=0; j<globalCluster[i].GetSize(); j++){
      fOutputEventMask->SetFlag(globalCluster[i].GetImpulseID(j),1);
    }
  }
  for (int i=0; i<100; i++) {if (trackCounter[i]>0) nTracksInGlobalCluster++;}
  if (globalCluster.size()>2) {
    h_gloClu_ntracks->Fill( nTracksInGlobalCluster,1);
    //    std::cout<<"ntacks in cluster: "<<nTracksInGlobalCluster<<"    response tracks: "<<fMCEvent->GetResponseMuonsN()<<"     ratio: "<<(float)nTracksInGlobalCluster/fMCEvent->GetResponseMuonsN()<<std::endl;
    h_debug->Fill((float)nTracksInGlobalCluster/fMCEvent->GetResponseMuonsN(),1);
    //    std::cout<<"size of cluster:  "<<globalCluster.size()<<" ;  clustered strings:   ";
    //    for (int i=0; i<globalCluster.size(); i++) std::cout<<globalCluster[i].GetStringID()<<"   ";
    //    std::cout<<std::endl;
  }
  
  //  if (fMCEvent->GetResponseMuonsN()==1&&nStrings>2&&nHits>5&&firedString[7]>1) std::cout<<"fired strings: "<<nStrings<<"    cluster: "<<globalCluster.size()<<"   "<<nStr<<std::endl;
  
  return kTRUE;
}



std::vector<BStringCluster> BClusterProducer::buildStringClusters(int iString, std::vector<int> string_impulses)
{
  std::vector<BStringCluster> result;
  
  if (string_impulses.size()==0) return result;				   
  
  ///////////////////TEMPORARY
  /*
  float primThetaRad=M_PI*fMCEvent->GetPrimaryParticlePolar()/180;
  float primPhiRad=M_PI*fMCEvent->GetPrimaryParticleAzimuth()/180;

  TVector3 genVec(sin(primThetaRad)*cos(primPhiRad),
		  sin(primThetaRad)*sin(primPhiRad),
		  cos(primThetaRad));

  TVector3 inPoTrack(fMCEvent->GetTrack(0)->GetX()-500*genVec.X(), fMCEvent->GetTrack(0)->GetY()-500*genVec.Y(), fMCEvent->GetTrack(0)->GetZ()-500*genVec.Z());

  TVector3 zero(0,0,0);
  float rho=getTrackDistanceToOM(inPoTrack, genVec, zero);
  */

  std::vector<BStringCluster> hot_spots=findHotSpots(iString, string_impulses);
  
  if (fVerbose) std::cout<<"string: "<<iString<<"    "<<hot_spots.size()<<" hotspots found,      size of string: "<<string_impulses.size()<<std::endl;

  for (int i=0; i<hot_spots.size(); i++){

    addImpulses(&hot_spots[i], string_impulses);

    //remove clustered impulses from hotspots
    for (int j=i+1; j<hot_spots.size(); j++){
      bool drop=false;
      for (int q=0; q<hot_spots[j].GetSize(); q++){
	if (hot_spots[j].GetImpulseID(0)==hot_spots[i].GetImpulseID(q)) drop=true;
	if (hot_spots[j].GetImpulseID(1)==hot_spots[i].GetImpulseID(q)) drop=true;
      }
      if (drop) hot_spots.erase(hot_spots.begin()+j);
    }

    result.push_back(hot_spots[i]);
    
    //stop if hotspots have been clustered
    if (i>=hot_spots.size()) break;
  }

  return result;
}

std::vector<BStringCluster> BClusterProducer::findHotSpots(int iString, std::vector<int> string_impulses)
{
  //PRODUCE ALL POSSIBLE CASUAL CONNECTED PAIRS OF IMPULSES AT NEIGHBOURING MODULES 
  //ORDER PAIRS IN SUM OF IMPULSE AMPLITUDE
  
  std::vector<BStringCluster> hot_spots;
  std::vector<bool> isClustered;
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++) isClustered.push_back(false);  
  
  for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
    if (fEvent->GetImpulse(string_impulses[iPulse])->GetAmplitude()<fSignalCut_hotspot1) continue;
    int chID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
    float timePulse=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
    for (int iCand=0; iCand<string_impulses.size(); iCand++){
      int id_increment=fEvent->GetImpulse(string_impulses[iCand])->GetChannelID()-chID;
      if (fVerbose) std::cout<<"            increment: "<<id_increment<<"    "<<chID<<"   "<<fEvent->GetImpulse(string_impulses[iCand])->GetChannelID()<<std::endl;
      if (abs(id_increment)==1&&
	  fEvent->GetImpulse(string_impulses[iCand])->GetAmplitude()>fSignalCut_hotspot2){
	float timeCand=fEvent->GetImpulse(string_impulses[iCand])->GetTime();
	if (fabs(timeCand-timePulse)<(abs(id_increment)*15/cWater+fSafetyWindow)&&(!isClustered[iCand])){
	  std::vector<int> spot;
	  spot.push_back(string_impulses[iCand]);
	  spot.push_back(string_impulses[iPulse]);
	  BStringCluster strClu(iString, spot, fEvent, fGeomTel, true, fMCEvent);
	  hot_spots.push_back(strClu);
	  isClustered[iPulse]=true;
	}
      }
    }
  }
  
  //ORDER PAIRS IN SUM IMPULSE
  BStringCluster buf[hot_spots.size()];
  
  for (int iHotspot=0; iHotspot<hot_spots.size(); iHotspot++){
    if (fVerbose) std::cout<<"       hotspot: "<<hot_spots[iHotspot].GetImpulseID(0)<<"    "<<hot_spots[iHotspot].GetImpulseID(1)<<"    "<<hot_spots[iHotspot].GetConstituent(0)->GetChannelID()<<"   "<<hot_spots[iHotspot].GetConstituent(1)->GetChannelID()<<"    ampl: "<<hot_spots[iHotspot].GetHotSpotAmpl()<<std::endl;
    float sumAmpl=hot_spots[iHotspot].GetSumAmpl();
    int nLarger=0;
    for (int i=0; i<hot_spots.size(); i++){
      float probeSumAmpl=hot_spots[i].GetSumAmpl();
      if (sumAmpl<probeSumAmpl) nLarger++;
    }
    buf[nLarger]=hot_spots[iHotspot];
  }
  
  std::vector<BStringCluster> result;
  for (int i=0; i<hot_spots.size(); i++) result.push_back(buf[i]);  
  
  return result;
 
}


int BClusterProducer::addImpulses(BStringCluster* hotspot, std::vector<int> string_impulses)
{
  std::vector<BImpulse*> storeys_pulse;

  for (int i=0; i<24; i++) storeys_pulse.push_back(NULL);
  
  for (int i=0; i<hotspot->GetSize(); i++) storeys_pulse[(hotspot->GetConstituent(i)->GetChannelID())%24]=hotspot->GetConstituent(i);
  
  for (int iSeed=0; iSeed<2; iSeed++){
    bool added=false;
    float seed_time=hotspot->GetConstituent(iSeed)->GetTime();
    int seed_channel_id=hotspot->GetConstituent(iSeed)->GetChannelID();
    int seed_storey=(seed_channel_id)%24;
    float adjacent_time=hotspot->GetConstituent(1-iSeed)->GetTime();
    int adjacent_channel_id=hotspot->GetConstituent(1-iSeed)->GetChannelID();
    
    float deltaZ=fabs(fGeomTel->At(seed_channel_id)->GetZ()-fGeomTel->At(adjacent_channel_id)->GetZ());
    float deltaT=seed_time-adjacent_time;
    int id_increment=seed_channel_id-adjacent_channel_id;

    //add +-1 channels
    float time_early=seed_time+deltaT-fSafetyWindow; 
    float time_late=max(seed_time+deltaZ/cWater+fSafetyWindow, seed_time+deltaT+fSafetyWindow);

    if (fVerbose) std::cout<<" deltaT: "<<deltaT<<"   seed_time: "<<seed_time<<"    seed chan: "<<seed_channel_id<<"   t_early: "<<time_early<<"   t_late: "<<time_late<<std::endl;
    
    for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
      if (string_impulses[iPulse]==hotspot->GetImpulseID(iSeed)||
	  string_impulses[iPulse]==hotspot->GetImpulseID(1-iSeed)) continue;
      int chanID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
      float chan_time=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
      if (fVerbose) std::cout<<"time: "<<chan_time<<"       chanID: "<<chanID<<std::endl;
      if (chanID==seed_channel_id+id_increment){
	if (chan_time>time_early&&chan_time<time_late){
 	  hotspot->AddImpulse(string_impulses[iPulse]);
	  storeys_pulse[chanID%24]=fEvent->GetImpulse(string_impulses[iPulse]);
	}
      }
    }
  
    if (fVerbose){
      if (hotspot->GetSize()==2) std::cout<<"fail"<<std::endl;
      if (hotspot->GetSize()==0) std::cout<<"super fail"<<std::endl;
      if (hotspot->GetSize()>2) std::cout<<"ok"<<std::endl;
    }

    if (!added) continue;
    
    //add up to 7 channels
    
    int id[]={1,3};
    for (int k=0; k<2; k++){
      int i=id_increment*id[k];
      if (storeys_pulse[seed_storey+i]==NULL) continue;
      for (int j=1; j<3; j++){
	//find min and max time for new channel
	if (seed_storey+i+sgn(i)*j<0||seed_storey+i+sgn(i)*j>23) continue;
	float tim0=storeys_pulse[seed_storey+i]->GetTime()+j*(storeys_pulse[seed_storey+i]->GetTime()-storeys_pulse[seed_storey+i-sgn(i)]->GetTime());
	float tim1=storeys_pulse[seed_storey+i]->GetTime()+j*(storeys_pulse[seed_storey+i]->GetTime()-storeys_pulse[seed_storey+i-2*sgn(i)]->GetTime())/2;
	float tim2=storeys_pulse[seed_storey+i-sgn(i)]->GetTime()+(j+1)*(storeys_pulse[seed_storey+i-sgn(i)]->GetTime()-storeys_pulse[seed_storey+i-2*sgn(i)]->GetTime());
	
	float min_estimate=min(tim0,tim1);
	min_estimate=min(min_estimate,tim2);
	
	float max_estimate=max(tim0,tim1);
	max_estimate=max(max_estimate,tim2);
	
	//loop over string pulses
	for (int iPulse=0; iPulse<string_impulses.size(); iPulse++){
	  int chanID=fEvent->GetImpulse(string_impulses[iPulse])->GetChannelID();
	  if (chanID!=seed_channel_id+i+sgn(i)*j) continue;
	  float time=fEvent->GetImpulse(string_impulses[iPulse])->GetTime();
	  if (time<max_estimate+fSafetyWindow&&time>min_estimate-fSafetyWindow){
	    hotspot->AddImpulse(string_impulses[iPulse]);
	    storeys_pulse[seed_storey+i+sgn(i)*j]=fEvent->GetImpulse(string_impulses[iPulse]);
	  } 
	}
      }
    }
    
  }

  if (fVerbose){
    std::cout<<"size of assembled cluster: "<<hotspot->GetSize()<<std::endl;
    for (int i=0; i<hotspot->GetSize(); i++){
      std::cout<<"         chanID: "<<hotspot->GetConstituent(i)->GetChannelID()<<std::endl;
    }
  }
  return 1;  
}


std::vector<BStringCluster> BClusterProducer::buildGlobalCluster(std::vector<BStringCluster>* stringClusters)
{
  std::vector<BStringCluster> result;
  
  //TRY TO LINK CLUSTERS ON DIFFERENT STRINGS:
  //FIND CLUSTER WITH LARGEST AMPLITUDE HOTSPOT (WITH LARGEST AMPLITUDE FOR THE MOMENT)
  float maxAmpl=0;
  int stringID=-1;
  int clusterID=-1;

  float minTheta=-1000;
  float maxTheta=1000;
  
   //select events which have clusters on the central string
  //consider all clusters at the central string select global cluster with maximum size
  int iString=7;
  stringID=7;
  for (int iCluster=0; iCluster<stringClusters[stringID].size(); iCluster++){
    stringID=7;
    //find min and max theta:
    /*
    if (stringClusters[stringID][iCluster].GetSize()>4){
      std::pair<float,float> minmax=getPolarEstimate_string( stringClusters[stringID][iCluster]);
      if (minmax.first>minTheta) minTheta=minmax.first;
      if (minmax.second<maxTheta) maxTheta=minmax.second;
    }
    */
    clusterID=iCluster;
    //    if (stringClusters[iString][iCluster].GetHotSpotAmpl()>maxAmpl/*&&stringClusters[iString][iCluster].GetSize()>5*/){
    //      maxAmpl=stringClusters[iString][iCluster].GetHotSpotAmpl();
    //      clusterID=iCluster;
    //    }
    
    
  //  float minTheta=1000;
  //  float maxTheta=-1000;
    std::vector<BStringCluster> globalCluster;
    globalCluster.push_back(stringClusters[stringID][clusterID]);
    //  for (int iString=0; iString<8; iString++){
    //  if (stringClusters[iString].size()==0) continue;
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
    
    if (fVerbose) std::cout<<"seed cluster: "<<stringID<<"    size: "<<stringClusters[stringID][clusterID].GetSize()<<"    "<<minTheta<<"   "<<maxTheta<<"       gen angle: "<<M_PI*(180-fMCEvent->GetPrimaryParticlePolar())/180<<std::endl;
    
    //check other strings
    //link if deltaT~deltaS/cVacuum
    //
    float seedTime=stringClusters[stringID][clusterID].GetCenterTime();
    float seedZ=stringClusters[stringID][clusterID].GetCenterZ();
    //get other coord. of the cluster
    float seedX=stringClusters[stringID][clusterID].GetX();
    float seedY=stringClusters[stringID][clusterID].GetY();
    
    for (int iStringLink=0; iStringLink<8; iStringLink++){
      if (iStringLink==stringID) continue;
      for (int iClusterLink=0; iClusterLink<stringClusters[iStringLink].size(); iClusterLink++){
	float linkTime=stringClusters[iStringLink][iClusterLink].GetCenterTime();
	float linkZ=stringClusters[iStringLink][iClusterLink].GetCenterZ();
	float linkX=stringClusters[iStringLink][iClusterLink].GetX();
	float linkY=stringClusters[iStringLink][iClusterLink].GetY();
	float deltaR=sqrt(pow(seedX-linkX,2)+pow(seedY-linkY,2)+pow(seedZ-linkZ,2));
	
	int direction=(linkZ-seedZ)/fabs(linkZ-seedZ);
	//	std::cout<<"linking: "<<fabs(linkTime-seedTime)<<"  times: "<<linkTime<<"   "<<seedTime<<"   "<<deltaR<<"   "<<deltaR/cVacuum<<std::endl;
	
	if (fabs(linkTime-seedTime)<deltaR/cVacuum - 200 || fabs(linkTime-seedTime)>deltaR/cVacuum + 200) continue;
	//	std::cout<<"linking: "<<fabs(linkTime-seedTime)<<"  times: "<<linkTime<<"   "<<seedTime<<"   "<<deltaR<<"   "<<deltaR/cVacuum<<std::endl;
	globalCluster.push_back(stringClusters[iStringLink][iClusterLink]);
      }
    }
    //    std::cout<<"got cluster: "<<globalCluster.size()<<std::endl;
    if (globalCluster.size()>2){
      //check if chosen clusters are on the same straight line.
      //      std::cout<<"cluster bigger than 2: "<<globalCluster.size()<<"   seedZ: "<<seedZ<<std::endl;
      for (int iClu=1; iClu<globalCluster.size(); iClu++){
	//	std::cout<<"      vec0: "<<globalCluster[iClu].GetStringID()<<"   xyz: "<<globalCluster[iClu].GetX()<<" "<<globalCluster[iClu].GetY()<<" "<<globalCluster[iClu].GetCenterZ()<<std::endl;
	TVector3 vec0(seedX-globalCluster[iClu].GetX(), seedY-globalCluster[iClu].GetY(), seedZ-globalCluster[iClu].GetCenterZ());
	for (int iCluConn=iClu+1; iCluConn<globalCluster.size(); iCluConn++){
	  TVector3 vec1(seedX-globalCluster[iCluConn].GetX(), seedY-globalCluster[iCluConn].GetY(), seedZ-globalCluster[iCluConn].GetCenterZ());
	  float angle=vec0.Angle(vec1);
	  //	  std::cout<<"              vec1: "<<globalCluster[iCluConn].GetStringID()<<"   xyz: "<<globalCluster[iCluConn].GetX()<<" "<<globalCluster[iCluConn].GetY()<<" "<<globalCluster[iCluConn].GetCenterZ()<<"     angle: "<<180*angle/M_PI<<std::endl;
	  if (angle < M_PI*0.8) {
	    //	    globalCluster.erase(globalCluster.begin()+iCluConn);
	  }
	  //	  else std::cout<<"take.   angle between vectors: "<<180*angle/M_PI<<std::endl;
	}
      }
    }
    if (globalCluster.size()>2) {
      result=globalCluster;
      break;
    }
  }
  //  std::cout<<"ghkhjk  "<<globalCluster.size()<<std::endl;
  return result;
}


Int_t BClusterProducer::PostProcess()
{
  fOUT->cd();
  h_ntracks->Write();
  h_ntracks_6hits3strings->Write();
  h_debug->Write();
  h_strClu_ntracks->Write();
  h_gloClu_ntracks->Write();
  h_pulse_ntracks->Write();
  h_strClu_noiseFrac->Write();
  fOUT->Close();
  return kTRUE;
}

float BClusterProducer::getTrackDistanceToOM(TVector3 initialPoint, TVector3 direction, TVector3 xyzOM)
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
/*
int BClusterProducer::runWPscan(std::vector<int>* string_impulses, std::vector<bool> noiseOMs_gen, std::vector<bool> signalOMs_gen, std::vector<bool> signalStrings_gen)
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
*/
int BClusterProducer::sgn(float x)
{
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}


int BClusterProducer::findTriggerString(float high_threshold, float low_threshold)
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


int BClusterProducer::countMCStringGaps(std::vector<int> string_impulses_gen)
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

std::pair<float,float> BClusterProducer::getPolarEstimate_string(BStringCluster stringCluster)
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
  
  for (int i = 0; i < stringCluster.GetSize(); i++) {
    int chID=stringCluster.GetConstituent(i)->GetChannelID();
    for (int j = i+1; j < stringCluster.GetSize(); j++) {
      int chIDnext=stringCluster.GetConstituent(i)->GetChannelID();
      if (abs(chID-chIDnext)!=1) continue;
      float time=stringCluster.GetConstituent(i)->GetTime();
      float timeNext=stringCluster.GetConstituent(i)->GetTime();
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
