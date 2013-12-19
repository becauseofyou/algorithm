
struct point {
	double x,y;
	void in() {
		scanf("%lf%lf",&x,&y);
	}
	point (double x=0,double y=0):x(x),y(y){}
};
typedef point Vector;
inline Vector operator + (const Vector &a, const Vector &b) {
	return Vector(a.x + b.x, a.y + b.y);
}
inline Vector operator - (const Vector &a, const Vector &b) {
	return Vector(a.x - b.x, a.y - b.y);
}
inline Vector operator * (const Vector &a, double t) {
	return Vector(a.x * t, a.y * t);
}
inline Vector operator / (const Vector &a, double p) {
	return Vector(a.x / p, a.y / p);
}
inline int sgn(double x,double eps=1e-8) {
	return x < -eps ? -1 : x > eps;
}
bool operator < (const point &a, const point &b) {
	return sgn(a.x - b.x) < 0 || sgn(a.x - b.x) == 0 && sgn(a.y - b.y) < 0;
}
bool operator == (const point &a, const point &b) {
	return sgn(a.x - b.x) == 0 && sgn(a.y - b.y) == 0;
}
inline double cross(Vector a, Vector b) {
	return a.x * b.y - a.y * b.x;
}
inline double dot(Vector a, Vector b) {
	return a.x * b.x + a.y * b.y;
}
bool intersect(point P, Vector v, point Q, Vector w, point &p) {
	Vector u = P - Q;
	if(sgn(cross(v, w)) == 0) return false;
	double t = cross(w, u) / cross(v, w);
	p = P + v * t;
	return true;
}
bool dotOnseg(point p,seg L)  {
	return sgn(cross(L.s - p, L.e - p)) == 0 && sgn(dot(L.s - p, L.e - p)) <= 0;
}
