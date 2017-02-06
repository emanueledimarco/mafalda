#ifndef TRACKS_ANALYZER_HH
#define TRACKS_ANALYZER_HH

#include <LinearTrackFinder.hh>

#include "TFile.h"
#include "TH1F.h"
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
  bool _saveFigs;

};

#endif
