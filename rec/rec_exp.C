void rec_exp(){
  const char * path_rec = "/data/rec";
  const char * path_proc = "/data/processed";
  const int min_strings = 3;
  const int min_channels = 6;
  const char *config_path = gSystem->ExpandPathName("$BARSSYS/config");
  const char *channel_mask_path = Form("%s/2015_channelmask_vs_e0556",
                                       config_path);
  const char *file_in = "e0556_caus100_type3_AT3PE";
  const char *file_in_path = Form("%s/%s_off52.root", path_proc, file_in);
  const char *file_out_path = Form("%s/%s_off63_rec.root", path_rec, file_in);
  const char *filter_mask_name = "CausalityFilterMask";

  MEvtLoop magic;
  MParList plist;
  MTaskList tasks;
  BReadTree events_reader("Events", file_in_path, "EventReader", "EventReader");
  BReadTree config_reader("ConfigTel", file_in_path, "ConfigurationReader",
                          "ConfigurationReader");
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
  rec_qual.SetHitCriterionFlag();
  rec_qual.SetFilterEventMaskName(filter_mask_name);
  MWriteRootFile writer(file_out_path, "RECREATE", "Magic root-file", 2);
  writer.AddContainer("BConfigTel", "ConfigTel");
  writer.AddContainer("BGeomTel", "Events");
  writer.AddContainer("BEvent", "Events");
  writer.AddContainer("BChannelMask", "Events");
  writer.AddContainer(filter_mask_name, "Events");
  writer.AddContainer("BJointExtractedHeader", "Events");
  writer.AddContainer("BRecParameters", "Events");

  plist.AddToList(&tasks);
  tasks.AddToList(&events_reader);
  tasks.AddToList(&config_reader);
  tasks.AddToList(&chmask);
  tasks.AddToList(&trigcut);
  tasks.AddToList(&muonrec);
  tasks.AddToList(&rec_qual);
  tasks.AddToList(&writer);
  magic.SetParList(&plist);
  magic.Eventloop();
}
