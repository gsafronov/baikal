void compareMeth()
{
  TFile* fAll=new TFile("mcAnalysis_allHits.root");
  TFile* fDirect=new TFile("mcAnalysis_directHits.root");

  TH1F* hCluOffsets_all=(TH1F*)fAll->Get("hClusterOffsets");
  TH1F* hCluOffsets_direct=(TH1F*)fDirect->Get("hClusterOffsets");

  TH1F* hMethUncert=new TH1F("hMethUNcert","hMethUncert",1000,-10,10);
  
  for (int i=0; i<hCluOffsets_all->GetNbinsX(); i++){
    if (hCluOffsets_direct->GetBinContent(i+1)!=0) hMethUncert->Fill(hCluOffsets_direct->GetBinContent(i+1)-hCluOffsets_all->GetBinContent(i+1),1);
  }
  hMethUncert->Draw();
}
