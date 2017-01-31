#include <vector>
#include <utility>

typedef std::pair<double,double> Hit;
typedef std::vector<Hit> HitCollection;
struct SimpleTrack {
  double pars[2];
  double errors[2];
  double chi2;
  HitCollection hitsInPattern;
};
typedef std::vector<SimpleTrack> SimpleTrackCollection;

class LinearTrackFinder {
  
public:

  LinearTrackFinder(double x1, double y1, double x2, double y2);
  LinearTrackFinder();

  void loadHits(std::vector<Hit> hits) { _hits = hits; }
  void setHitUncertainty(double a) { _hitUnc = a; }

private:

  double distance(Hit one, Hit two) { 
    double x1 = one.first;
    double y1 = one.second;
    double x2 = two.first;
    double y2 = two.second;
    return sqrt(pow((x2-x1),2)+pow((y2-y1),2)); 
  }

  HitCollection getInitialHits(double mindist);
  SimpleTrack getTrack(HitCollection hits);
  void updateTrack(SimpleTrack &t, double x, double xsize, double ysize, double chi2max, double nhitsmin);

  double _x1,_y1,_x2,_y2;
  HitCollection _hits;
  HitCollection _hitsInWindow;
  double _hitUnc;
  SimpleTrackCollection _trackCollection;
  
};
