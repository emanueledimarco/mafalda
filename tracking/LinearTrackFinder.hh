#ifndef LINEAR_TRACK_FINDER_HH
#define LINEAR_TRACK_FINDER_HH

#include <vector>
#include <utility>

typedef std::pair<double,double> Hit;
typedef std::pair<Hit,Hit> Seed;
typedef std::vector<Seed> SeedCollection;
typedef std::vector<Hit> HitCollection;
struct SimpleTrack {
  double pars[2];
  double errors[2];
  double chi2;
  bool good;
  HitCollection hitsInPattern;
};
typedef std::vector<SimpleTrack> SimpleTrackCollection;

class LinearTrackFinder {
  
public:

  LinearTrackFinder(double x1, double y1, double x2, double y2);
  LinearTrackFinder();

  void loadHits(HitCollection hits);
  void setHitUncertainty(double a) { _hitUnc = a; }
  void setSeedingMinHitDistance(double d) { _minDist = d; }
  void setSearchWindowSize(double x, double y) { _xsize = x; _ysize = y; }
  void setnMinHits(int n) { _nHitsMin = n; }
  void setMaxTrackAttempts(int n) { _maxTrackAttempts = n; }
  void setTrackCoeffRange(double min, double max) {_cmin=min; _cmax=max; }
  void setDebugLevel(int n) { _debugLevel = n; }
  int getResidualHits() { return _hits.size(); }
  SimpleTrackCollection makeTracks();

private:

  double distance(Hit one, Hit two) { 
    double x1 = one.first;
    double y1 = one.second;
    double x2 = two.first;
    double y2 = two.second;
    return sqrt(pow((x2-x1),2)+pow((y2-y1),2)); 
  }

  HitCollection getInitialHits();
  SimpleTrack getTrack(HitCollection hits);
  void updateTrack(SimpleTrack &t, double x, double xsize, double ysize, double chi2max, int nhitsmin);
  bool hitInTrack(Hit hit, SimpleTrack t);
  int nPossibleSeeds(int nhits);

  double _x1,_y1,_x2,_y2;
  HitCollection _hits;
  SeedCollection _usedSeeds;
  double _hitUnc;
  double _minDist;
  double _xsize, _ysize;
  double _cmin, _cmax;
  int _nHitsMin;
  int _maxTrackAttempts;
  int _debugLevel;
  
};

#endif
