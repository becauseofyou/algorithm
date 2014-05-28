/* **********************************************
Author      : wuyiqi
Created Time: 2014年04月05日 星期六 11时11分57秒
File Name   : A.cpp
 *********************************************** */
#include <set>
#include <map>
#include <list>
#include <stack>
#include <queue>
#include <cmath>
#include <deque>
#include <bitset>
#include <cstdio>
#include <vector>
#include <string>
#include <complex>
#include <sstream>
#include <utility>
#include <climits>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>
const int MAX_CONVEX = 500;
const int MAX_ITEM = 500;
struct point {
    double x,y;
    void in() {
        scanf("%lf%lf",&x,&y);
    }
    point (double x=0,double y=0):x(x),y(y){}
};
struct seg {
    point s, e;
    seg(point s, point e) : s(s) , e(e) {}
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
inline double cross(point o, point a, point b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
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
//bool segcross(point a, point b, point c, point d)  {  
//    int d1,d2,d3,d4;  
//    d1=sgn(cross(a,b,c));  
//    d2=sgn(cross(a,b,d));  
//    d3=sgn(cross(c,d,a));  
//    d4=sgn(cross(c,d,b));  
//    if((d1^d2)==-2&&(d3^d4)==-2)  
//        return true;  
//    return false;  
//}  
bool segcross(point p1, point p2, point q1, point q2) {
    return (
            std::min(p1.x, p2.x) <= std::max(q1.x, q2.x) &&
            std::min(q1.x, q2.x) <= std::max(p1.x, p2.x) &&
            std::min(p1.y, p2.y) <= std::max(q1.y, q2.y) &&
            std::min(q1.y, q2.y) <= std::max(p1.y, p2.y) && /* 跨立实验 */
            cross(p1, q2, q1) * cross(p2, q2, q1) <= 0 && 
            cross(q1, p2, p1) * cross(q2, p2, p1) <= 0  /* 叉积相乘判方向 */
           );
}



#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
using namespace std;
const double eps = 1e-8;
const double pi = acos(-1.0);

struct Point {
    double x, y;
    Point operator - (const Point& t) const {
        Point tmp;
        tmp.x = x - t.x;
        tmp.y = y - t.y;
        return tmp;
    }
    Point operator + (const Point& t) const {
        Point tmp;
        tmp.x = x + t.x;
        tmp.y = y + t.y;
        return tmp;
    }
    bool operator == (const Point& t) const {
        return fabs(x-t.x) < eps && fabs(y-t.y) < eps;
    }
}GP; 
struct Cir {
	Point ct;
    double r;
};

inline double Cross(Point a, Point b, Point c) {					// 叉积
	return (b.x-a.x)*(c.y-a.y) - (c.x-a.x)*(b.y-a.y); 
}
inline double PPdis(Point a, Point b) {								// 点点距离
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
inline double PLdis(Point p,Point l1,Point l2){						// 点线距离
	return fabs(Cross(p,l1,l2))/PPdis(l1,l2);
}
inline bool same_dir(Point a, Point b) {							// 向量是否同向
    return fabs(a.x*b.y-b.x*a.y) < eps && a.x*b.x > -eps && a.y*b.y > -eps;
}
bool dotOnSeg(Point p, Point s, Point e) {							// 点是否在线段上
    if ( p == s || p == e )		// 看具体情况端点是否合法
        return true;
    return fabs(Cross(p,s,e)) < eps && 
        (p.x-s.x)*(p.x-e.x)<eps && (p.y-s.y)*(p.y-e.y)<eps;
}
bool Intersect(Point p1, Point p2, Point p3, Point p4, Point& p) {	// 直线相交
	double a1, b1, c1, a2, b2, c2, d;
	a1 = p1.y - p2.y; b1 = p2.x - p1.x; c1 = p1.x*p2.y - p2.x*p1.y;
	a2 = p3.y - p4.y; b2 = p4.x - p3.x; c2 = p3.x*p4.y - p4.x*p3.y;
	d = a1*b2 - a2*b1;
	if ( fabs(d) < eps )	return false;
	p.x = (-c1*b2 + c2*b1) / d;
	p.y = (-a1*c2 + a2*c1) / d;
	return true;
}


/***** ax+by+c=0 ***********************************************************************/
struct Line {
    double a, b, c;
};
double PL_dis(Point p, Line ln) {									// 点线距离
    return fabs(ln.a*p.x+ln.b*p.y+ln.c)/sqrt(ln.a*ln.a+ln.b*ln.b);
}
Line Turn(Point s, Point e) {										// 线段转直线表达式
    Line ln;
    ln.a = s.y - e.y;
    ln.b = e.x - s.x;
    ln.c = s.x*e.y - e.x*s.y;
    return ln;
}
bool Line_Inst(Line l1, Line l2, Point &p) {						// 直线相交 
    double d = l1.a*l2.b - l2.a*l1.b;
    if ( fabs(d) < eps )    return false;
    p.x = (-l1.c*l2.b + l2.c*l1.b) / d;
    p.y = (-l1.a*l2.c + l2.a*l1.c) / d;
    return true; 
}
bool Cir_Line(Point ct, double r, Line ln, Point& t1, Point& t2) {   // 直线与圆交点
    if ( PL_dis(ct, ln) > r + eps )
        return false;
    ln.c += ln.a*ct.x + ln.b*ct.y;
    if ( fabs(ln.b) < eps ) {
        t1.x = t2.x = -ln.c/ln.a;
        t1.y = sqrt(r*r - ln.c*ln.c/ln.a/ln.a);
        t2.y = -t1.y;
    } 
    else {
        double A, B, C;
        A = ln.a*ln.a + ln.b*ln.b;
        B = 2.0*ln.a*ln.c;
        C = ln.c*ln.c - ln.b*ln.b*r*r;
        t1.x = (-B - sqrt(B*B - 4.0*A*C))/2.0/A;
        t2.x = (-B + sqrt(B*B - 4.0*A*C))/2.0/A;
        t1.y = -t1.x*ln.a/ln.b - ln.c/ln.b;
        t2.y = -t2.x*ln.a/ln.b - ln.c/ln.b;
    }
    t1 = t1 + ct;
    t2 = t2 + ct;
    return true;
}
bool Cir_Cir(Point c1, double r1, Point c2, double r2, Point& t1, Point& t2) {  // 圆与圆交点
	double d = PPdis(c1, c2);
	if ( d > r1+r2+eps || d < fabs(r1-r2)-eps )
		return false;
    Line ln;
    ln.a = 2*(c1.x - c2.x);
    ln.b = 2*(c1.y - c2.y);
    ln.c = r1*r1 - r2*r2 - (c1.x*c1.x+c1.y*c1.y-c2.x*c2.x-c2.y*c2.y);
    Cir_Line(c1, r1, ln, t1, t2);
    return true;
}
/***************************************************************************************/


/***** 三角形 **************************************************************************/
bool Orthocenter(Point a, Point b, Point c, Point& p) {				// 垂心
	Point ta = a, tb = b;
	ta.x = a.x + b.y - c.y;
	ta.y = a.y + c.x - b.x;
	tb.x = b.x + a.y - c.y;
	tb.y = b.y + c.x - a.x;
	return Intersect(a, ta, b, tb, p);
}
bool Circumcenter(Point& u, Point& v, Point& w, Point& p) {			// 外接圆心
	Point a, b, c, d;
	a.x = (u.x + w.x) / 2;
	a.y = (u.y + w.y) / 2;
	b.x = a.x + u.y - w.y;
	b.y = a.y - u.x + w.x;

	c.x = (v.x + w.x) / 2;
	c.y = (v.y + w.y) / 2;
	d.x = c.x + v.y - w.y;
	d.y = c.y - v.x + w.x;
	return Intersect(a, b, c, d, p);
}
/***************************************************************************************/


/***** 圆 ******************************************************************************/
bool cir_line(Point ct, double r, Point l1, Point l2, Point& p1, Point& p2) {// 直线与圆
	if ( PLdis(ct, l1, l2) > r + eps )
		return false;
	double a1, a2, b1, b2, A, B, C, t1, t2;
	a1 = l2.x - l1.x; a2 = l2.y - l1.y;
	b1 = l1.x - ct.x; b2 = l1.y - ct.y;
	A = a1*a1 + a2*a2;
	B = (a1*b1 + a2*b2)*2;
	C = b1*b1 + b2*b2 - r*r;
	t1 = (-B - sqrt(B*B - 4.0*A*C))/2.0/A;
	t2 = (-B + sqrt(B*B - 4.0*A*C))/2.0/A;
	p1.x = l1.x + a1*t1; p1.y = l1.y + a2*t1;
	p2.x = l1.x + a1*t2; p2.y = l1.y + a2*t2;
	return true;
}
bool cir_cir(Point c1, double r1, Point c2, double r2, Point& p1, Point& p2) {// 圆与圆
	double d = PPdis(c1, c2);
	if ( d > r1+r2+eps || d < fabs(r1-r2)-eps )
		return false;
	Point u, v;
	double t=(1 + (r1*r1-r2*r2)/PPdis(c1,c2)/PPdis(c1,c2))/2;
	u.x = c1.x + (c2.x-c1.x)*t;
	u.y = c1.y + (c2.y-c1.y)*t;
	v.x = u.x + c1.y - c2.y;
	v.y = u.y + c2.x - c1.x;
	cir_line(c1, r1, u, v, p1, p2);
	return true;
}
double cir_area_inst(Point c1, double r1, Point c2, double r2) {			// 两圆面积交
	double a1, a2, d, ret;
    d = sqrt((c1.x-c2.x)*(c1.x-c2.x)+(c1.y-c2.y)*(c1.y-c2.y));
    if ( d > r1 + r2 - eps ) 
        return 0;
    if ( d < r2 - r1 + eps ) 
        return pi*r1*r1;
	if ( d < r1 - r2 + eps )
		return pi*r2*r2;
    a1 = acos((r1*r1+d*d-r2*r2)/2/r1/d);
    a2 = acos((r2*r2+d*d-r1*r1)/2/r2/d);
    ret = (a1-0.5*sin(2*a1))*r1*r1 + (a2-0.5*sin(2*a2))*r2*r2;
	return ret;
}
struct Ball { double x, y, z, r; };
double Dis(Ball p, Ball q) {												// 两球体积交
	return sqrt((p.x-q.x)*(p.x-q.x) + (p.y-q.y)*(p.y-q.y) + (p.z-q.z)*(p.z-q.z));
}
double cal(double a, double b, double r) {
	return pi*((r*r*b - b*b*b/3.0)-(r*r*a - a*a*a/3.0));
}
double Solve(Ball a, Ball b) {
	double Va, Vb, dis, ret;
	if ( a.r > b.r ) {
		Ball tmp = a; a = b; b = tmp;
	}
	dis = Dis(a, b);
	Va = 4.0/3*pi*(a.r*a.r*a.r);
	Vb = 4.0/3*pi*(b.r*b.r*b.r);
	if ( dis > a.r + b.r - eps ) 
		ret = Va + Vb;
	else if ( dis < b.r - a.r + eps ) 
		ret = Vb;
	else {
		double t = (a.r*a.r - b.r*b.r) / dis;
		double la = (dis + t) / 2.0;
		double lb = (dis - t) / 2.0;
		ret = Va + Vb - cal(la, a.r, a.r) - cal(lb, b.r, b.r);
	}
	return ret;
}
/***************************************************************************************/


/***** 多边形 **************************************************************************/
Point GravityCenter(Point* p, int n) {								// 多边形重心
	int i;
	double A=0, a;
	Point t;
	t.x = t.y = 0;
	p[n] = p[0];
	for (i=0; i < n; i++) {
		a = p[i].x*p[i+1].y - p[i+1].x*p[i].y;
		t.x += (p[i].x + p[i+1].x) * a;
		t.y += (p[i].y + p[i+1].y) * a;
		A += a;
	}
	t.x /= A*3;
	t.y /= A*3;
	return t;
}
bool point_in_polygon(Point o, Point* p, int n) {					// 点是否在多边形内
	int i, t;
	Point a, b;
	p[n] = p[0];
	for (i=0; i < n; i++) {
		if ( dotOnSeg(o, p[i], p[i+1]) )
			return true;
	}
	t = 0;
	for (i=0; i < n; i++) {
		a = p[i]; b = p[i+1];
		if ( a.y > b.y ) {
			Point tmp = a; a = b; b = tmp;
		}
		if ( Cross(o, a, b) < -eps && a.y < o.y-eps && o.y < b.y+eps )
			t++;
	}
	return t&1;
}
/***************************************************************************************/


/***** 凸包 ****************************************************************************/
bool cmpyx(Point a, Point b) {										
	if ( a.y != b.y )
		return a.y < b.y;
	return a.x < b.x;
}
void Grahamxy(Point *p, int &n) {									// 水平序(住:两倍空间)
	if ( n < 3 ) 
		return;
	int i, m=0, top=1;
	sort(p, p+n, cmpyx);	
	for (i=n; i < 2*n-1; i++)
		p[i] = p[2*n-2-i];
	for (i=2; i < 2*n-1; i++) {
		while ( top > m && Cross(p[top], p[i], p[top-1]) < eps ) 
			top--;
		p[++top] = p[i];
		if ( i == n-1 )	m = top;
	}
	n = top;
}

bool cmpag(Point a, Point b) {
    double t = (a-GP).x*(b-GP).y - (b-GP).x*(a-GP).y;
    return fabs(t) > eps ? t > 0 : PPdis(a, GP) < PPdis(b, GP);
}
void Grahamag(Point *p, int &n) {									// 极角序                    
    int i, top = 1;
    GP = p[0];
    for (i=1; i < n; i++) if(p[i].y<GP.y-eps || (fabs(p[i].y-GP.y)<eps && p[i].x<GP.x)) {
        GP = p[i];
    }
    sort(p, p+n, cmpag);
    for ( i=2; i < n; i++ ) {
        while ( top > 0 &&  Cross(p[top], p[i], p[top-1]) < eps )
            top--;
        p[++top] = p[i];
    }
    p[++top] = p[0];
    n = top;
}
/***************************************************************************************/
