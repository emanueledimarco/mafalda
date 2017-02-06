#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TString.h"
#include "TStyle.h"

#include <TracksAnalyzer.hh>

TracksAnalyzer::TracksAnalyzer() :
  _saveFigs(false) {
  _fileout = 0;
  _h1ds.clear();
}

void TracksAnalyzer::beginJob() {

  _fileout = TFile::Open("TrackAnalyzerOutput.root","recreate");

  TH1F *nTracks = new TH1F("nTracks","",20,0,20);
  TH1F *nHits = new TH1F("nHits","",30,0,30);
  TH1F *chi2 = new TH1F("chi2","",20,0,5);
  TH1F *direction = new TH1F("direction","",30,-2,2);

  _h1ds[std::string(nTracks->GetName())]=nTracks;
  _h1ds[std::string(nHits->GetName())]=nHits;
  _h1ds[std::string(chi2->GetName())]=chi2;
  _h1ds[std::string(direction->GetName())]=direction;

}

void TracksAnalyzer::analyze(SimpleTrackCollection tracks) {

  _h1ds["nTracks"]->Fill((int)tracks.size());

  for(SimpleTrackCollection::const_iterator tk=tracks.begin();tk!=tracks.end();++tk) {
    _h1ds["nHits"]->Fill(tk->hitsInPattern.size());
    _h1ds["chi2"]->Fill(tk->chi2);
    _h1ds["direction"]->Fill(tk->pars[1]);
  }

}

void TracksAnalyzer::endJob() {
  
  gStyle->SetOptStat(11111);

  TCanvas* c1 = new TCanvas("c1","",600,600);

  _fileout->cd();
  for(H1DCollection::const_iterator histo=_h1ds.begin(); histo!=_h1ds.end(); ++histo) {
    TH1F *h1d=histo->second;
    h1d->Write();
    if(_saveFigs) {
      h1d->SetLineColor(kRed+1);
      h1d->GetXaxis()->SetTitle(h1d->GetName());
      h1d->Draw();
      c1->SaveAs(Form("%s.pdf",h1d->GetName()));
    }
  }
  _fileout->Close();

}
