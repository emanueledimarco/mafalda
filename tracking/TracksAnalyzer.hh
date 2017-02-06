#ifndef TRACKS_ANALYZER_HH
#define TRACKS_ANALYZER_HH

#include <LinearTrackFinder.hh>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include <map>

class TracksAnalyzer {

  typedef std::map<std::string,TH1F*> H1DCollection;

public:
  TracksAnalyzer();
  
  void beginJob();
  void analyze(SimpleTrackCollection tracks);
  void endJob();
  
  void setSaveFigs(bool b) {_saveFigs = b;}

private:
  
  H1DCollection _h1ds;
  TFile* _fileout;
  TTree* _tree;
  bool _saveFigs;
  int _nTracks;
  double _nHits[1000], _chi2[1000], _direction[1000];

};

#endif
