#include "MDirIter.h"
#include "MEvtLoop.h"
#include "MParList.h"
#include "MTaskList.h"
#include "TSystem.h"

#include "BReadTree.h"
#include "BTelConfigReader.h"
#include "BEvoGeomApply.h"
#include "BChannelMaskApply.h"
#include "BReadCalibPlain.h"
#include "BCalibrateAmpCharge.h"
#include "BCalibrateTime.h"
#include "BCausality.h"
#include "BMuonCriterion.h"
#include "BOfflineTriggerCut.h"
#include "BFQualify.h"
#include "BRecQualify.h"
#include "BReconstructMuon.h"
#include "MWriteRootFile.h"
#include "BExtractMCEvent.h"
#include "BAmplitudeFilter.h"
#include "BMCRead.h"
#include "BReadMCGeomWout.h"


void timeCalibstage1_mc(string mc_infile="/home/local1/work/baikal/mc/cors-n8x10m-2sec-tr112-thr4.0-1.5.wout", string root_infile="OUT/outout.root",  int calibChanID=10, Bool_t isrec = kTRUE)
{
	// Script is used to get MC events after filtrations (amplitude, causality, muon) or even more after muon reconstruction
	//
	// Usual execution: root -l 'rec/stage1_mc.C+(fnamein, fnameout, 1)'. Tested execution only from macros/ directory
	//
	// Data are available at baikalsnr.jinr.ru:/data/mc
	//
	// Steps:
	// - get MC events (wout file),
	// - attach static MC geometry, 
	// - attach channel mask
	// - generates common calibrated events (BEvent) 
	// - apply amplitude filter
	// - apply causality criterion 
	// - apply muon criterion 
	// - apply offline trigger cut 
	// - gather parameters after filtration
	// - if isrec flag is activated apply muon reconstruction
	// - if isrec flag is activated gather MC parameters after reconstruction

	gSystem->Load("../libmars");
		
    MParList plist;
    MTaskList tasks; 
    plist.AddToList(&tasks);
    
    // Read MC-events from binary file
    BMCRead reader(mc_infile.c_str());
    tasks.AddToList(&reader);
    
    //Initialize geometry
    //BGeomApply geom;
    //tasks.AddToList(&geom);
    //geom.SetVerbose(1);
    //geom.SetGeometry("BGeomTel2015MC");
    
    BReadMCGeomWout geom;
    geom.SetNSectionsInString(2);
    geom.SetVerbose(kTRUE);
    tasks.AddToList(&geom);
    
    TString config_path = "../config/";
    BChannelMaskApply chmask;
    chmask.SetMaskFileName(config_path + "2015_channelmask_vs_e0556");
    chmask.SetVerbose(kTRUE);
    tasks.AddToList(&chmask);

    BMyRecoChanSetter chSet(calibChanID);
    tasks.AddToList(&chSet);
    
    //Unify mc-events 
    BExtractMCEventSum eventExtractor;
    tasks.AddToList(&eventExtractor);
    
    BAmplitudeFilter ampfilt;
    ampfilt.SetAmpThreshold(3);
    tasks.AddToList(&ampfilt);
    
    BCausality causFilter1("CAUS1");
    causFilter1.SetCausType(3);
    causFilter1.SetCausLimit(50);
    causFilter1.SetInputMaskName("AmplitudeFilterMask");
    causFilter1.SetOutputMaskName("CausalityFilterMask");
    //causFilter1.SetVerbose(kTRUE);
    tasks.AddToList(&causFilter1);
    
    BMuonCriterion muonCriterion("muonCriterion", "Muon Criterion");
    muonCriterion.SetInputMaskName("CausalityFilterMask");
    muonCriterion.SetOutputMaskName("MuonCriterionFilterMask");	
    muonCriterion.FilterSinglePulses(kFALSE);
    //muonCriterion.SetFilterSinglePulsesOnString(kTRUE);
    tasks.AddToList(&muonCriterion);
    
    BOfflineTriggerCut trigcut;
    trigcut.SetMinStringNum(2);
    trigcut.SetMinChannelNum(5);
    trigcut.SetMaskName("MuonCriterionFilterMask");
    tasks.AddToList(&trigcut);
    
    BFQualify qual;
    qual.SetFilterEventMaskName("MuonCriterionFilterMask");
    tasks.AddToList(&qual);
    
    BReconstructMuon muonrec;
    muonrec.SetEventMaskName("MuonCriterionFilterMask");
    muonrec.SetInitialsType(0);
    //muonrec.SetVerbose(kTRUE);
    
    BRecQualify recqual;
    recqual.SetHitCriterionFlag();
    recqual.SetQ10Path(gSystem->ExpandPathName("$BARSSYS/config/quant10HE.data"));
    recqual.SetFilterEventMaskName("MuonCriterionFilterMask");
    
    if(isrec) {
      tasks.AddToList(&muonrec);
      tasks.AddToList(&recqual);
    }
    
    BMyRecoReco myrecoTest("foutMCreco.root", true, calibChanID);
    tasks.AddToList(&myrecoTest);
    
    
    MWriteRootFile writer(root_infile.c_str(), "RECREATE", "Magic root-file", 2); 
    tasks.AddToList(&writer);
    writer.AddContainer("BMCArrayConfig", "ArrayConfig");    
    writer.AddContainer("BMCEvent", "Events");
    writer.AddContainer("BEvent", "Events");
    writer.AddContainer("MCEventSource", "Events");
    writer.AddContainer("MCEventMask", "Events");
    writer.AddContainer("CausalityFilterMask", "Events");
    writer.AddContainer("MuonCriterionFilterMask", "Events");
    writer.AddContainer("BMuonCriterionParameters", "Events");
    writer.AddContainer("BFQualityEvent", "Events");
    writer.AddContainer("BGeomTel", "Events");
    if(isrec) {
      writer.AddContainer("BRecParameters", "Events");
    }
    
    //Init and run event loop
    MEvtLoop eventLoop;
    eventLoop.SetParList(&plist);
    eventLoop.Eventloop();
}
