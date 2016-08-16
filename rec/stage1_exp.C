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
#include "BReconstructMuon.h"
#include "MWriteRootFile.h"
#include "BExtractExperimentalEvent.h"
#include "BAmplitudeFilter.h"
#include "TEnv.h"

//void stage1_exp(const char *fnamein = "/data/joint/e0558_joint.root", const char *fnameout="OUT/e0558_amp3_caus50_type3_muon_crit_rec_off52.root", Bool_t isrec = kTRUE)
void stage1_exp(int nrun, const char *pathin = "/data", const char *pathout = "/data/rec", int season = 2015, Bool_t isrec = kTRUE)
{
	// Script is used to get experimental events after filtrations (amplitude, causality, muon) or even more after muon reconstruction
	//
	// Usual execution: root -l 'rec/stage1_exp.C+(fnamein, fnameout, 1)'. Tested execution only from macros/ directory
	//
	// Data are available at baikalsnr.jinr.ru:/data/joint
	//
	// Steps:
	// - get joint events (BExtractedImpulseTel),
	// - attach dynamic geometry, 
	// - attach charge and time calibration, 
	// - attach channel mask
	// - generates common calibrated events (BEvent) 
	// - apply amplitude filter
	// - apply causality criterion 
	// - apply muon criterion 
	// - apply offline trigger cut 
	// - gather parameters after filtration
	// - if isrec flag is activated apply muon reconstruction	
	// PS channel mask could be modified in charge and time calib tasks
	
	TEnv env(gSystem->ExpandPathName("$BARSSYS/config/rec.rc"));	
	cout << "gSystem->ExpandPathName = " << gSystem->ExpandPathName("$BARSSYS/config/rec.rc") << endl;
	const Float_t ampthr      = env.GetValue("AmplitudeFilter.AmpThreshold", 3);     // amplitude threshold in p.e.
	if(!env.Defined("AmplitudeFilter.AmpThreshold")) {
		cout << "Parameter AmplitudeFilter.AmpThreshold is not defined in environment" << endl;
		return;
	}
	const Float_t causlimit   = env.GetValue("Causality.CausLimit", 50); // causality window in ns
    const Int_t   caustype    = env.GetValue("Causality.CausType", 3);     // causality type
    const Int_t   minstrnum   = env.GetValue("OfflineTriggerCut.MinStringNum", 3);
    const Int_t   minchnum    = env.GetValue("OfflineTriggerCut.MinChannelNum", 5);	
	cout << env.Defined("AmplitudeFilter.AmpThreshold") << endl;
	return ;
	
	const char *prefix;
	if(season == 2015) {
		prefix = "e";
	}
	else if(season == 2016) {
		prefix = "f";		
	}
	else {
		cout << "wrong season. Abort." << endl;
		exit(1);
	}
	
	TString fnamein = Form("%s/joint/%s%04d_joint.root",  pathin, prefix, nrun);
	
	gSystem->Load("../libmars");
	MParList plist;
    MTaskList tasks;
    plist.AddToList(&tasks);
	
	BReadTree reader("Events", fnamein.Data()); // read BExtractedImpulseTel
    tasks.AddToList(&reader);

	BTelConfigReader telreader("TelescopeConfigurationReader");
	TString confname = Form("%s/processed/%s%04d.root",  pathin, prefix, nrun);
	telreader.SetSource(confname.Data(), "Configuration"); // file with MainInfo
	tasks.AddToList(&telreader);
	
	// dynamic geometry 
	BEvoGeomApplyLinear geom;
	TString geomname = Form("%s/acoustic/%d/%04d/%s%04d.evo.beacons.aligned.guessed.root",  pathin, season, nrun, prefix, nrun);
	geom.SetGeometryData(geomname.Data());
	geom.SetStringTemplate("../config/stringtemplate_2sec");
	//geom.SetVerbose(kTRUE);
	tasks.AddToList(&geom);
	
	// fixed geometry 
	//BGeomApply geom;
	//tasks.AddToList(&geom);
	//geom.SetGeometry("BGeomTel2015");

	// Channel Masking
	TString config_path = "../config/";
	BChannelMaskApply chmask;
	//chmask.SetMaskFileName(config_path + "2015_channelmask_vs_str8");
	chmask.SetMaskFileName(config_path + "2015_channelmask_vs_e0556");
	chmask.SetVerbose(kTRUE);
   	tasks.AddToList(&chmask);

	BReadCalibPlain read_tcalib("TCalibReader", "TCalibReader", config_path + "2015_tcalib_ver2");
	BReadCalibPlain read_acalib("ACalibReader", "ACalibReader", config_path + "2015_acalib");
	BReadCalibPlain read_pcalib("PCalibReader", "PCalibReader", config_path + "2015_pcalib");
	BReadCalibPlain read_qcalib("QCalibReader", "QCalibReader", config_path + "2015_qcalib_ver2");
	BReadCalibPlain read_offcalib("OffsetCalibReader", "OffsetCalibReader", config_path + "2015_offsets_556");
	read_tcalib.SetCalibTelName("TCalibTel");
	read_acalib.SetCalibTelName("ACalibTel");
	read_pcalib.SetCalibTelName("PCalibTel");
	read_qcalib.SetCalibTelName("QCalibTel");
	read_offcalib.SetCalibTelName("OffCalibTel");
	tasks.AddToList(&read_tcalib);
	tasks.AddToList(&read_acalib);
	tasks.AddToList(&read_pcalib);
	tasks.AddToList(&read_qcalib);	
	tasks.AddToList(&read_offcalib);
	
	BCalibrateAmpCharge aq_calib;	
	aq_calib.SetNameAmplitudeTel("ACalibTel");
	aq_calib.SetNameChargeTel("QCalibTel");
	aq_calib.SetNamePedestalTel("PCalibTel");
	tasks.AddToList(&aq_calib);

	BCalibrateTime t_calib;
	t_calib.SetNameTimeTel("TCalibTel");
	t_calib.SetNameOffsetTel("OffCalibTel");
	tasks.AddToList(&t_calib);

	BExtractExperimentalEvent beventxtractor; 
	//beventxtractor.SetVerbose(kTRUE);
	tasks.AddToList(&beventxtractor);
		   	   
	BAmplitudeFilter ampfilt;
    ampfilt.SetAmpThreshold(ampthr);
	tasks.AddToList(&ampfilt);	
	
    //  
    BCausality causFilter1("CAUS");
    causFilter1.SetInputMaskName("AmplitudeFilterMask");
    causFilter1.SetOutputMaskName("CausalityFilterMask");
    causFilter1.SetCausLimit(causlimit);
    causFilter1.SetCausType(caustype);
    //causFilter1.SetVerbose(1);
    tasks.AddToList(&causFilter1);
    
    BMuonCriterion muoncrit("muonCriterion", "Muon Criterion");	
	muoncrit.SetInputMaskName("CausalityFilterMask");
	muoncrit.SetOutputMaskName("MuonCriterionFilterMask");
	muoncrit.SetFilterSinglePulsesOnString(kFALSE);
	tasks.AddToList(&muoncrit);
    
    BOfflineTriggerCut trigcut;
    trigcut.SetMinStringNum(minstrnum);
    trigcut.SetMinChannelNum(minchnum);
    trigcut.SetMaskName("MuonCriterionFilterMask");
    tasks.AddToList(&trigcut);

    BFQualify qual;
    qual.SetFilterEventMaskName("MuonCriterionFilterMask");
    tasks.AddToList(&qual);
    
	BReconstructMuon muonrec;
	muonrec.SetEventMaskName("MuonCriterionFilterMask");
	muonrec.SetInitialsType(0);
	//muonrec.SetVerbose(1);
	
	if(isrec) {
		tasks.AddToList(&muonrec);
	}
	
	// building a name for output file
	const char *suffix;
	if(isrec) {
		suffix = "rec";
	}
	else {
		suffix = "filt";
	}		
	TString fnameout = Form("%s/%s%04d_%s_amp%d_caus%d_type%d_off%d%d.root",  pathout, prefix, nrun, suffix, int(ampthr), int(causlimit), caustype, minstrnum, minchnum);
	
    MWriteRootFile writer(fnameout.Data(), "RECREATE", "Magic root-file", 2); 
    tasks.AddToList(&writer);
    writer.AddContainer("BConfigTel", "ConfigTel");
    writer.AddContainer("BEvent", "Events");
    writer.AddContainer("CausalityFilterMask", "Events");
    writer.AddContainer("MuonCriterionFilterMask", "Events");
    writer.AddContainer("BJointExtractedHeader", "Events");
    writer.AddContainer("BFQualityEvent", "Events");
    writer.AddContainer("BMuonCriterionParameters", "Events");
    writer.AddContainer("AmplitudeFilterMask", "Events");
    writer.AddContainer("BGeomTel", "Events");
    if(isrec) {
		writer.AddContainer("BRecParameters", "Events");
	}
    
	// event loop
    MEvtLoop magic;
    magic.SetParList(&plist);
    magic.Eventloop();
    
    env.Write();
}
