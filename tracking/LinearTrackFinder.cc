#include <random>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>

#include <LinearTrackFinder.hh>

LinearTrackFinder::LinearTrackFinder(double x1, double y1, double x2, double y2) :
  _x1(x1),
  _y1(y1),
  _x2(x2),
  _y2(y2),
  _minDist(0),
  _xsize(1),
  _ysize(1),
  _nHitsMin(2),
  _maxTrackAttempts(5),
  _debugLevel(0) {
  _hits.clear();
  _usedSeeds.clear();
  }

LinearTrackFinder::LinearTrackFinder(){}

void LinearTrackFinder::loadHits(HitCollection hits) {
  _hits.clear();
  for (HitCollection::const_iterator hit=hits.begin(); hit!=hits.end(); ++hit) {
    if(hit->first > _x1 && hit->first < _x2 &&
       hit->second > _y1 && hit->second < _y2) _hits.push_back(*hit);
  }
}

HitCollection LinearTrackFinder::getInitialHits() {


  if(_debugLevel>0) std::cout << "Getting the seeding hits of the track" << std::endl;
  
  HitCollection out;
  int nHits = _hits.size();

  if(nHits<2) return out;

  bool newSeed=false;
  int nTrial=0;
  int random_one, random_two;
  double dist = 0;
  while (newSeed==false && nTrial<5) {
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(0,nHits-1); // guaranteed unbiased
    random_one = uni(rng);
    
    random_two = 0;
    while(dist<_minDist) {
      random_two = uni(rng);
      dist = distance(_hits[random_one],_hits[random_two]);
    }
    Seed thisSeed = std::make_pair(_hits[random_one],_hits[random_two]);
    SeedCollection::const_iterator sitr = std::find(_usedSeeds.begin(), _usedSeeds.end(), thisSeed);
    if (sitr == _usedSeeds.end()) newSeed = true;
    nTrial++;
  }

  out.push_back(_hits[random_one]);
  out.push_back(_hits[random_two]);
  if(_debugLevel>0) std::cout << "Got the following hits of the track: " << std::endl
                              << _hits[random_one].first << " , " << _hits[random_one].second << std::endl
                              << _hits[random_two].first << " , " << _hits[random_two].second << std::endl
                              << "with a distance of " << dist << std::endl;
  return out;
}

SimpleTrack LinearTrackFinder::getTrack(HitCollection hits) {
  if(_debugLevel>0) std::cout << "Make the linear fit" << std::endl;
  SimpleTrack track;
  TGraphErrors g(hits.size());
  for (int i=0; i<(int)hits.size(); ++i) {
    g.SetPoint(i,hits[i].first,hits[i].second);
    g.SetPointError(i,_hitUnc,_hitUnc);
    track.hitsInPattern.push_back(hits[i]);
  }
  
  // TCanvas *c1 = 0;
  // if (_debugLevel>2) {
  //   c1 = new TCanvas("c1","",600,600);
  //   g.Draw("APE");
  // }
  g.Fit("pol1");
  TF1 *func = g.GetFunction("pol1");
  for(int i=0; i<2; ++i) {
    track.pars[i] = func->GetParameter(i);
    track.errors[i] = func->GetParError(i);
  }
  track.chi2 = func->GetChisquare();
  if(_debugLevel>2) {
    TCanvas c1("c1","",600,600);
    g.Draw("APE");
    c1.SaveAs(Form("track_nHits%d_offs%f_coeff%f.pdf",(int)track.hitsInPattern.size(),track.pars[0],track.pars[1]));
  }
  return track;
}

void LinearTrackFinder::updateTrack(SimpleTrack &t, double x, double xsize, double ysize, double chi2max, int nhitsmin) {
  int initialHitsInPattern = t.hitsInPattern.size();
  double y = t.pars[0] + t.pars[1]*x;
  if(_debugLevel>3) {
    std::cout << "*** UPDATING TRACK WITH THE FOLLOWING HITS:" << std::endl;
    for (HitCollection::const_iterator thit = t.hitsInPattern.begin(); thit!=t.hitsInPattern.end(); ++thit) 
      std::cout << "\t\t" << thit->first << " , " << thit->second << std::endl;
  }

    
  
  for(HitCollection::const_iterator hit=_hits.begin(); hit<_hits.end(); ++hit) {
    if(_debugLevel>3) std::cout << "\tHit " << hit->first << " , " << hit->second << std::endl
                                << "\txwindow = " << x-xsize << " , " << x+xsize << std::endl
                                << "\tywindow = " << y-ysize << " , " << y+ysize << std::endl
                                << "\tHit belonging to track: " << hitInTrack(*hit,t) << std::endl;
    
    if(hitInTrack(*hit,t)==false &&
       hit->first > x-xsize && hit->first < x+xsize &&
       hit->second > y-ysize && hit->second < y+ysize) {
      t.hitsInPattern.push_back(*hit);
      if(_debugLevel>2) std::cout << "Added new hit to the hit pattern: " << hit->first << " , " << hit->second << std::endl;
    } else { if(_debugLevel>3) std::cout << "Hit NOT ADDED to pattern: " << hit->first << " , " << hit->second << std::endl; }
    
  }

  if (t.hitsInPattern.size() == initialHitsInPattern) return;

  if(_debugLevel>2) { 
    std::cout << "updating the track at point x = " << x << std::endl
              << "and linear parameters: " << t.pars[0] << ", " << t.pars[1] << std::endl
              << "starting with a track with " << initialHitsInPattern << " hits" << std::endl
              << "Found an additional hit along the candidate track. Refitting track. " << std::endl;
  }

  SimpleTrack tmpTrack = getTrack(t.hitsInPattern);

  std::cout << tmpTrack.chi2 << "   " << chi2max << "  " << (int)tmpTrack.hitsInPattern.size() << " " << nhitsmin << std::endl;
  t.good = (tmpTrack.chi2 < chi2max && (int)tmpTrack.hitsInPattern.size() >= nhitsmin);

  for(int i=0; i<2; ++i) {
    t.pars[i] = tmpTrack.pars[i];
    t.errors[i] = tmpTrack.pars[i];
  }
  t.chi2 = tmpTrack.chi2;
  return;
}

SimpleTrackCollection LinearTrackFinder::makeTracks() {

  SimpleTrackCollection out;

  if(_debugLevel > 0 ) std::cout << "Inital hits of the event = " << _hits.size() << std::endl;

  int failedTracks=0;
  while (_hits.size() > 3 && failedTracks < _maxTrackAttempts) {

    if(_debugLevel>0) std::cout << "Track collection size = " << out.size() << std::endl;
    HitCollection seedHits = getInitialHits();

    //    while (failedTracks < _maxTrackAttempts && _hits.size() > 3) {
      if(_debugLevel>2) std::cout << "LinearTrackFinder::makeTracks. Iterating over the hits. Residual hits = " << _hits.size() << std::endl
                                  << "max attempts for this track = " << failedTracks << std::endl;
  
      SimpleTrack track = getTrack(seedHits);
      for(double x=_x1; x<_x2; x+=_xsize) updateTrack(track,x,_xsize,_ysize,5,_nHitsMin);
      if (track.good) {
        std::cout << "\tGOOD TRACK!" << std::endl;
        out.push_back(track);
        failedTracks=0;
        
        std::cout << "%%% This track has hits: " << std::endl;
        for(HitCollection::const_iterator thit=track.hitsInPattern.begin(); thit!=track.hitsInPattern.end(); ++thit) {
          std::cout << "\t\t" << thit->first << " , " << thit->second << std::endl;
        }

        std::cout << "\tLOOKING FOR HITS TO BE ERASED (size = " << _hits.size() << "):" << std::endl;
        for(HitCollection::const_iterator ahit=_hits.begin(); ahit<_hits.end();) {
          std::cout << "\tCONSIDERING THIS HIT: " << ahit->first << " , " << ahit->second << std::endl; 
          if(hitInTrack(*ahit,track)) {
            std::cout << "\t\tERASING THIS  WHICH BELONGS TO A GOOD TRACK " << std::endl;
            ahit = _hits.erase(ahit);
          } else {
            std::cout << "\t\tLEAVING THIS HIT" << std::endl;
            ++ahit;
          }
        }
        // good track created. Re-seed from the residual hits
        seedHits = getInitialHits();
      } else {
        failedTracks++;
      } // if too-many attempts of track have been done for this seed, re-seed
  }
  return out;
}

bool LinearTrackFinder::hitInTrack(Hit hit, SimpleTrack t) {
  HitCollection::const_iterator posInPattern = std::find(t.hitsInPattern.begin(), t.hitsInPattern.end(), hit);
  return (posInPattern != t.hitsInPattern.end());
}
