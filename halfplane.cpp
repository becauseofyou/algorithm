point intersect(point P, Vector v, point Q, Vector w) {
  point p;
  Vector u = P - Q;
  double t = cross(w, u) / cross(v, w);
  p = P + v * t;
  return p;
}
int inhalfplane(point p,point s,point e) {
  return sgn(cross(e - s, p - s)) ;
}
std::vector<point> CUT(const std::vector<point> &p, point s, point e) {
  std::vector<point> q;
  int n = (int) p.size();
  for(int i = 0; i < n; i++) {
    int nowin = inhalfplane(p[i], s, e);
    int nextin = inhalfplane(p[(i + 1) % n], s, e);
    if(nowin >= 0) {
      q.push_back(p[i]);
    }
    if(nextin * nowin < 0) {
      q.push_back(intersect(p[i], p[(i + 1) % n] - p[i], s, e - s));
    }
  }
  return q;
}
