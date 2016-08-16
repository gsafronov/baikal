void rec_mc() {
  const char *path = "/data/mc/root";
  const char *file_in = "cors-n8x10m-2sec-tr112-thr4.0-1.5_caus100_type3_AT3PE";
  const int min_strings = 3;
  const int min_channels = 6;
  const char *config_path = gSystem->ExpandPathName("$BARSSYS/config");
  const char *channel_mask_path = Form("%s/2015_channelmask_vs_e0556",
                                       config_path);
  const char *file_out_path = Form("%s/%s_off%d%d_rec.root", path, file_in,
                              min_channels, min_strings);
  const char *file_in_path = Form("%s/%s_off52.root", path, file_in);
  const char *filter_mask_name = "CausalityFilterMask";
  TFile *file = TFile::Open(file_in_path);

  TTree *geom_tree = (TTree*)file->Get("GeomTel");
  BGeomTel *geom = 0;
  geom_tree->SetBranchAddress("BGeomTel.", &geom);
  geom_tree->GetEntry(0);
  geom->SetName("BGeomTel");
  geom->SetReadyToSave();

  TTree *config_tree = (TTree*)file->Get("ArrayConfig");
  TBranch *config_branch = config_tree->GetBranch("BMCArrayConfig.");
  BMCArrayConfig *config = 0;
  config_branch->SetAddress(&config);
  config_branch->GetEntry(0);
  config->SetReadyToSave();

  file->Close();

  MEvtLoop magic;
  MParList plist;
  MTaskList tasks;
  BReadTree reader("Events", file_in_path);
  BChannelMaskApply chmask;
  chmask.SetMaskFileName(channel_mask_path);
  chmask.SetVerbose(kTRUE);
  BOfflineTriggerCut trigcut;
  trigcut.SetMinStringNum(min_strings);
  trigcut.SetMinChannelNum(min_channels);
  trigcut.SetMaskName(filter_mask_name);
	BReconstructMuon muonrec;
	muonrec.SetEventMaskName(filter_mask_name);
	muonrec.SetInitialsType(0);
  BRecQualify rec_qual;
  rec_qual.SetFilterEventMaskName(filter_mask_name);
  rec_qual.SetHitCriterionFlag();
  MWriteRootFile writer(file_out_path, "RECREATE", "Magic root-file", 2);
  writer.AddContainer("BMCArrayConfig", "ArrayConfig");
  writer.AddContainer("BGeomTel", "GeomTel");
  writer.AddContainer("BEvent", "Events");
  writer.AddContainer("BChannelMask", "Events");
  writer.AddContainer("MCEventMask", "Events");
  writer.AddContainer("MCEventSource", "Events");
  writer.AddContainer("BMCEvent", "Events");
  writer.AddContainer(filter_mask_name, "Events");
  writer.AddContainer("BFQualityEvent", "Events");
  writer.AddContainer("BRecParameters", "Events");

  plist.AddToList(geom);
  plist.AddToList(config);
  plist.AddToList(&tasks);
  tasks.AddToList(&reader);
  tasks.AddToList(&chmask);
  tasks.AddToList(&trigcut);
  tasks.AddToList(&muonrec);
  tasks.AddToList(&rec_qual);
  tasks.AddToList(&writer);
  magic.SetParList(&plist);
  magic.Eventloop();
}
