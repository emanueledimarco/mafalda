#include <random>
#include <algorithm>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TAxis.h>
#include <TMath.h>

#include <LinearTrackFinder.hh>

LinearTrackFinder::LinearTrackFinder(double x1, double y1, double x2, double y2) :
  _x1(x1),
  _y1(y1),
  _x2(x2),
  _y2(y2),
  _minDist(0),
  _xsize(1),
  _ysize(1),
  _cmin(-1000),
  _cmax(1000),
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
  double seedSlope = 0;
  while (newSeed==false && nTrial<10) {
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(0,nHits-1); // guaranteed unbiased
    random_one = uni(rng);
    
    random_two = 0;
    int nComb=0;
    while(dist < _minDist && nComb < 5*nHits ) {
      nComb++;
      random_two = uni(rng);
      if(random_two==random_one) continue;
      dist = distance(_hits[random_one],_hits[random_two]);
      int one, two;
      if(_hits[random_one].first < _hits[random_two].first) {
        one = random_one; two = random_two;
      } else {
        one = random_two; two = random_one;
      }
      double dx = _hits[two].first - _hits[one].first;
      double dy = _hits[two].second - _hits[one].second;
      seedSlope = dy/dx;
      if(seedSlope < _cmin || seedSlope > _cmax) continue;
      if(_debugLevel>2) std::cout << "Seed slope = " << seedSlope << std::endl;
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
                              << "with a distance of " << dist << std::endl
                              << "and a seed slope of " << seedSlope << std::endl;
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
  
  TF1 *func = new TF1("line","pol1",_x1,_x2);
  func->SetParameter(1,0.5*(_cmin+_cmax));
  g.Fit("line","Q");
  for(int i=0; i<2; ++i) {
    track.pars[i] = func->GetParameter(i);
    track.errors[i] = func->GetParError(i);
  }
  track.chi2 = func->GetChisquare();
  track.good = false; // decided only at the update level
  if(_debugLevel>2) {
    TCanvas c1("c1","",600,600);
    g.SetTitle(""); g.SetName("");
    g.Draw("APE");
    g.GetXaxis()->SetLimits(_x1,_x2);
    g.GetYaxis()->SetLimits(_y1,_y2);
    g.GetXaxis()->SetRangeUser(_x1,_x2);
    g.GetYaxis()->SetRangeUser(_y1,_y2);
    g.GetXaxis()->SetTitle("x");
    g.GetYaxis()->SetTitle("y");
    g.SetMarkerStyle(kFullSquare);
    g.SetMarkerColor(kBlack);
    g.Draw("APE");
    c1.Update();

    TGraphErrors gothers(_hits.size());
    for(int i=0; i<(int)_hits.size(); ++i) {
      if(!hitInTrack(_hits[i],track)) {
        gothers.SetPoint(i,_hits[i].first,_hits[i].second);
        gothers.SetPointError(i,_hitUnc,_hitUnc);
      }
    }
    gothers.SetTitle(""); gothers.SetName("");
    gothers.SetMarkerStyle(kFullSquare);
    gothers.SetMarkerColor(kGray);
    gothers.SetLineColor(kGray);
    gothers.Draw("PE");
    std::cout << "Track parameters:" << track.pars[0] << "\t"  << track.pars[1] << std::endl;
    c1.SaveAs(Form("track_nTotHits%d_nTrackHits%d_offs%03.3f_coeff%03.3f.pdf",(int)_hits.size(),(int)track.hitsInPattern.size(),track.pars[0],track.pars[1]));
  }
  delete func;
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

  t.good = (tmpTrack.chi2 < chi2max && (int)tmpTrack.hitsInPattern.size() >= nhitsmin &&
            t.pars[1] > _cmin && t.pars[1] < _cmax);
  if(_debugLevel>2) {
    std::cout << "--- Gooddness of track:" << std::endl;
    std::cout << "chi2 = " << tmpTrack.chi2 << " (" << chi2max << ")" <<  std::endl
              << "hits in track = " << (int)tmpTrack.hitsInPattern.size() << " (" << nhitsmin << ")" << std::endl;
  }

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
  int maxAttempts = _maxTrackAttempts;
  while (_hits.size() > 3 && failedTracks < maxAttempts) {

    if(_debugLevel>0) std::cout << "Track collection size = " << out.size() << std::endl;
    HitCollection seedHits = getInitialHits();

    if(_debugLevel>2) std::cout << "LinearTrackFinder::makeTracks. Iterating over the hits. Residual hits = " << _hits.size() << std::endl
                                << "max attempts for this track = " << failedTracks << std::endl;

    maxAttempts = std::min(_maxTrackAttempts,(int)_hits.size());
  
    SimpleTrack track = getTrack(seedHits);
    
    for(double x=_x1; x<_x2; x+=_xsize) updateTrack(track,x,_xsize,_ysize,10,_nHitsMin);
    if (track.good) {
      out.push_back(track);
      failedTracks=0;
 
      if(_debugLevel>1) {
        std::cout << "\tGOOD TRACK!" << std::endl;       
        std::cout << "%%% This track has hits: " << std::endl;
        for(HitCollection::const_iterator thit=track.hitsInPattern.begin(); thit!=track.hitsInPattern.end(); ++thit) {
          std::cout << "\t\t" << thit->first << " , " << thit->second << std::endl;
        }
      }
      if(_debugLevel>1) std::cout << "\tLOOKING FOR HITS TO BE ERASED (size = " << _hits.size() << "):" << std::endl;
      for(HitCollection::const_iterator ahit=_hits.begin(); ahit<_hits.end();) {
        if(_debugLevel>1) std::cout << "\tCONSIDERING THIS HIT: " << ahit->first << " , " << ahit->second << std::endl; 
        if(hitInTrack(*ahit,track)) {
          if(_debugLevel>1) std::cout << "\t\tERASING THIS  WHICH BELONGS TO A GOOD TRACK " << std::endl;
          ahit = _hits.erase(ahit);
        } else {
          if(_debugLevel>1) std::cout << "\t\tLEAVING THIS HIT" << std::endl;
          ++ahit;
        }
      }
      // good track created. Re-seed from the residual hits
      _usedSeeds.clear();
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

int LinearTrackFinder::nPossibleSeeds(int nhits) {
  return (int) 0.5 * TMath::Factorial(nhits) / TMath::Factorial(nhits - 2);
}
