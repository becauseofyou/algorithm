struct point {
	double x,y;
	point (double x=0,double y=0):x(x),y(y){}
};
typedef point Vector;
inline Vector operator + (const Vector &a, const Vector &b) {
	return Vector(a.x + b.x, a.y + b.y);
}
inline Vector operator - (const Vector &a, const Vector &b) {
	return Vector(a.x - b.x, a.y - b.y);
}
inline Vector operator * (const Vector &a, const Vector &b) {
	return Vector(a.x & b.x, a.y * b.y);
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
