#include <random>
#include <TGraphErrors.h>
#include <TF1.h>

#include <LinearTrackFinder.hh>

LinearTrackFinder::LinearTrackFinder(double x1, double y1, double x2, double y2) :
  _x1(x1),
  _y1(y1),
  _x2(x2),
  _y2(y2) {
  _hits.clear();
  _hitsInWindow.clear();
  }

LinearTrackFinder::LinearTrackFinder(){}

HitCollection LinearTrackFinder::getInitialHits(double mindist) {
                                                                
  for(std::vector<Hit>::iterator aHit=_hits.begin(); aHit<_hits.end(); ++aHit){
    if(aHit->first > _x1 && aHit->second > _y1 && aHit->first < _x2 && aHit->second < _y2) _hitsInWindow.push_back(*aHit);
  }
  
  int nHits = _hitsInWindow.size();
  
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(0,nHits); // guaranteed unbiased
  auto random_one = uni(rng);

  double dist = 0;
  auto random_two = 0;
  while(dist<mindist) {
    random_two = uni(rng);
    dist = distance(_hitsInWindow[random_one],_hitsInWindow[random_two]);
  }
  HitCollection out;
  out.push_back(_hitsInWindow[random_one]);
  out.push_back(_hitsInWindow[random_two]);
  return out;
}

SimpleTrack LinearTrackFinder::getTrack(HitCollection hits) {
  SimpleTrack track;
  TGraphErrors g(hits.size());
  for (int i=0; i<(int)hits.size(); ++i) {
    g.SetPoint(i,hits[i].first,hits[i].second);
    g.SetPointError(i,_hitUnc,_hitUnc);
    track.hitsInPattern.push_back(hits[i]);
    // _hits.erase(_hits.begin()+i);
  }
  
  g.Fit("pol1");
  TF1 *func = g.GetFunction("pol1");
  for(int i=0; i<2; ++i) {
    track.pars[i] = func->GetParameter(i);
    track.errors[i] = func->GetParError(i);
  }
  track.chi2 = func->GetChisquare();
  return track;
}

void LinearTrackFinder::updateTrack(SimpleTrack &t, double x, double xsize, double ysize, double chi2max, double nhitsmin) {
  double y = t.pars[0] + t.pars[1]*x;
  for(HitCollection::const_iterator hit=_hits.begin(); hit<_hits.end(); ++hit) {
    if(hit->first > x-xsize && hit->first < x+xsize &&
       hit->second > y-ysize && hit->second < y+ysize) t.hitsInPattern.push_back(*hit);
  }

  SimpleTrack tmpTrack = getTrack(t.hitsInPattern);

  if(x>_x2 && tmpTrack.chi2 < chi2max && tmpTrack.hitsInPattern.size()>=nhitsmin) {
    _trackCollection.push_back(tmpTrack);
    for(int th=0; th<(int)t.hitsInPattern.size(); ++th) {
      bool found=false;
      for(HitCollection::const_iterator ahit=_hits.begin(); ahit<_hits.end() && (!found); ++ahit) {
        if(ahit->first == t.hitsInPattern[th].first && ahit->second == t.hitsInPattern[th].second) {
          _hits.erase(_hits.begin()+th);
          found=true;
        }
      }
    }
    return;
  }

  for(int i=0; i<2; ++i) {
    t.pars[i] = tmpTrack.pars[i];
    t.errors[i] = tmpTrack.pars[i];
  }
  t.chi2 = tmpTrack.chi2;
  t.hitsInPattern = tmpTrack.hitsInPattern;
  return;
}
