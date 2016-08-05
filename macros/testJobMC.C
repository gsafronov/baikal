int testJobMC()
{
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
				      "*",
				      "TStreamerInfo",
				      "RIO",
				      "TStreamerInfo()");
  MParList plist;
  MTaskList tasks;
  plist.AddToList(&tasks);
  
  MReadTree  reader("Events", "OUT/mc_new_n8x10m_reco.root");
  reader.DisableAutoScheme();
  tasks.AddToList(&reader);
  
  //BRawMainInfo* maininfo = 0;
 
  //init configuration
  //  TFile f("OUT/2015/e0556.root", "READ");
  //  if (!f.IsOpen()) 
  //    return;
  
  //  TTree* t = (TTree*)f.Get("Configuration");
  // if (!t) 
  //  return 1;
  
  //  t->SetBranchAddress("BRawMainInfo.", &maininfo);
  //t->GetEntry(0);
  
  //  BTelConfigManager::Init(maininfo);
  
  // static geometry module
  //  BGeomApply geom;
  //  const char *data_type = "2015";
  //  geom.SetGeometry(Form("BGeomTel%s", data_type));
  //  tasks.AddToList(&geom);
  
  //dynamic geometry
  //  BEvoGeomApplyLinear geom;
  //  geom.SetGeometryData("aligned_beacon_data_Oct2015-Feb2016-guess_600s-align_60s.root");
  //  geom.SetStringTemplate("../config/stringtemplate_2sec");
  //  tasks.AddToList(&geom);
  
  BGeomApply geom;
  tasks.AddToList(&geom);
  geom.SetVerbose(1);
  geom.SetGeometry("BGeomTel2015MC");

  /*
  MWriteRootFile writer("foutMC_reco.root", "RECREATE", "Magic root-file", 2); 
  tasks.AddToList(&writer);
  writer.AddContainer("BRecParameters", "Events");  
  */

  BMyRecoReco myrecoTest("foutMC_reco.root", true);
  tasks.AddToList(&myrecoTest);

  MEvtLoop magic;
  magic.SetParList(&plist);
  magic.Eventloop(10000);


  return 0;
}


  //
  /*
    BJointMarkTelBuilder joint;
    char fnamediag[80];
    sprintf(fnamediag, "OUT/e%04d_diag.txt", nrun);
    ../libmars");
    
    int k;
    
    for(int i = nrun_min; i < nrun_max; i++) {
    k = ProcessRun(dirname, flag, data_type, i);
		switch(k) {
			case 777:
				cout << "Unknown data type!" << endl;
				return;
			break;
			case 888:
				cout << "Unknown process type!" << endl;
				return;
			break;
		}
	}
	
} 
*/

