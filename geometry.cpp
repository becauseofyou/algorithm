#include <iostream>
#include <fstream>
#include <cstring>
#include <climits>
#include <deque>
#include <cmath>
#include <queue>
#include <stack>
#include <list>
#include <map>
#include <set>
#include <utility>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <bitset>
#include <functional>
#include <algorithm>
typedef long long LL;
#define for_each(it, container)  for(__typeof(container.begin()) it = container.begin(); it != container.end(); ++it)
#define mp std::make_pair
#define pii std::pair<int, int>
#define sqr(x) ((x) * (x))
const double eps = 1e-8;
int sgn(double x) 
{
        return x < -eps ? -1 : x > eps;
}
struct Point;
typedef Point Vector;
struct Point 
{
	double x, y;
	void in() {
		scanf("%lf%lf", &x, &y);
	}
	void print() {
		printf("%.2lf %.2lf\n", x, y);
	}
	Point(double x = 0, double y = 0) : x(x), y(y) {
	}
        inline Vector rotate(double ang) {
                return Vector(x * cos(ang) - y * sin(ang), x * sin(ang) + y * cos(ang));
        }
        inline double dot(const Vector &a) {
                return x * a.x + y * a.y;
        }
        inline bool operator == (const Point &a) const {
                return sgn(x - a.x) == 0 && sgn(y - a.y) == 0;
        }
        inline bool operator < (const Point &a) const {
                return sgn(x - a.x) < 0 ¦¦ sgn(x - a.x) == 0 && sgn(y - a.y) < 0;
        }
        inline Vector operator + (const Vector &a) const {
                return Vector(x + a.x, y + a.y);
        }
        inline Vector operator - (const Vector &a) const {
                return Vector(x - a.x, y - a.y);
        }
        inline double operator * (const Vector &a) const {
                return x * a.y - y * a.x;
        }
        inline Vector operator * (double t) const {
                return Vector(x * t, y * t);
        }
        inline Vector operator / (double t) {
                return Vector(x / t, y / t);
        }
        inline double vlen() {
                return sqrt(x * x + y * y);
        }
};

struct Cir
{
        Point ct;
        double r;
        void in() {
                ct.in();
                scanf("%lf", &r);
        }
};
struct Seg 
{
        Point s, e;
        Seg() {
        }
        Seg(Point s, Point e): s(s), e(e) {
        }
        void in() {
                s.in();
                e.in();
        }
};
struct Line
{
        int a, b, c;
};

inline bool cmpyx(const Point &a, const Point &b) 
{
        if(a.y != b.y) {
                return a.y < b.y;
        }
        return a.x < b.x;
}

double cross(Point a, Point b, Point c) 
{
	return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
}
bool same_dir(Vector a, Vector b) 
{						
	return sgn(a.x * b.y - b.x * a.y) == 0 && sgn(a.x * b.x) >= 0 && sgn(a.y * b.y) >= 0;
}
bool dot_on_seg(Point p, Seg L) 
{
        return sgn((L.s - p) * (L.e - p)) == 0 && sgn((L.s - p).dot(L.e - p)) <= 0;
}
double ppdis(Point a, Point b) 
{							
	return sqrt((a - b).dot(a - b));
}
double pldis(Point p,Point l1,Point l2)
{					
	return fabs(cross(p,l1,l2))/ppdis(l1,l2);
}
double pldis(Point p, Line ln) 
{								
	return fabs(ln.a * p.x + ln.b * p.y + ln.c) / sqrt(ln.a * ln.a + ln.b * ln.b);
}
bool point_in_circle(Point &a, Cir cr) 
{
	return sgn(ppdis(a, cr.ct) - cr.r) <= 0;
}
bool intersect(Point P, Vector v, Point Q, Vector w, Point &ret) 
{
        Vector u = P - Q;
        if(sgn(v * w) == 0) return false;
        double t = w * u / (v * w);
        ret = P + v * t;
        return true;
}
Point intersect(Point P, Vector v, Point Q, Vector w)
{
        Point ret;
        Vector u = P - Q;
        if(sgn(v * w) == 0) return false;
        double t = w * u / (v * w);
        ret = P + v * t;
        return ret;
}
//点到直线距离
Point disptoline(Point p, Seg l)
{
        return fabs(cross(p, l.s, l.e)) / ppdis(l.s, l.e);
}
//点到直线的垂足(最近点)
Point ptoline(Point p, Seg l)
{
        Point vec = l.s - l.e;
        return intersect(p, vec.norm(), l.s, vec);
}
//点到线段的最近点
Point ptoseg(Point p, Seg l)
{
        Point norm = (l.s - l.e).norm();
        if(sgn(norm * (p - l.s)) * sgn(norm * (p - l.e)) > 0) {
                double sa = ppdis(p, l.s);
                double sb = ppdis(p, l.e);
                return sgn(sa - sb) < 0 ? l.s : l.e;
        }
        return intersect(p, norm, l.s, l.e - l.s);
}

bool segcross(Point p1, Point p2, Point q1, Point q2) 
{
	return (
			std::min(p1.x, p2.x) <= std::max(q1.x, q2.x) &&
			std::min(q1.x, q2.x) <= std::max(p1.x, p2.x) &&
			std::min(p1.y, p2.y) <= std::max(q1.y, q2.y) &&
			std::min(q1.y, q2.y) <= std::max(p1.y, p2.y) && /* 跨立实验 */
			cross(p1, q2, q1) * cross(p2, q2, q1) <= 0 && 
			cross(q1, p2, p1) * cross(q2, p2, p1) <= 0  /* 叉积相乘判方向 */
	       );
}

// 水平序, 注意两倍空间
struct Convex_Hull 
{
        static const int N = 100010;
        Point p[2 * N];
        int n;
        void init() {
                n = 0;
        }
        void in() {
                p[n].in();
                n++;
        }
        inline void push_back(const Point &np) {
                p[n++] = np;
        }
        void gao() {
                if(n < 3) {
                        return ;
                }
                std::sort(p, p + n);
                std::copy(p, p + n - 1, p + n);
                std::reverse(p + n, p + 2 * n - 1);
                int m = 0, top = 0;
                for(int i = 0; i < 2 * n - 1; i++) {
                        while(top >= m + 2 && sgn((p[top - 1] - p[top - 2]) * (p[i] - p[top - 2])) <= 0) {
                                top --;
                        }
                        p[top++] = p[i];
                        if(i == n - 1) {
                                m = top - 1;
                        }
                }
                n = top - 1;
        }
        void print() {
                for(int i = 0; i < n; i++) {
                        p[i].print();
                }
        }
        double get_area() {
                double ret = 0;
                Point ori(0, 0);
                for(int i = 0; i < n; i++) {
                        ret += (p[i] - ori) * (p[i + 1] - ori);
                }
                return fabs(ret) / 2;
        }
        double rotate() {
                if(n == 2) {
                        return (p[0] - p[1]).dot(p[0] - p[1]);
                }
                int i = 0, j = 0;
                for(int k = 0; k < n; k++) {
                        if(!(p[k] < p[i])) i = k;
                        if(!(p[j] < p[k])) j = k;
                }
                double ret = 0;
                int si = i, sj = j;
                while(i != sj || j != si) {
                        ret = std::max(ret, (p[i]-p[j]).dot(p[i]-p[j]));
                        if(sgn ( (p[(i + 1) % n] - p[i]) * (p[(j + 1) % n] - p[j]) ) < 0) {
                                i = (i + 1) % n;
                        } else {
                                j = (j + 1) % n;
                        }
                }
                return ret;
        }
}convex;

Line turn(Point s, Point e) 
{
	Line ln;
	ln.a = s.y - e.y;
	ln.b = e.x - s.x;
	ln.c = s.x * e.y - e.x * s.y;
	return ln;
}
//圆的折射，输入圆心，p->inter是射线，inter是圆上一点， ref是折射率，返回折射向量
Vector reflect_vector(Point center, Point p, Point inter, double ref) 
{
        Vector p1 = inter - p, p2 = center - inter;
        double sinang = p1 * p2 / (p1.vlen() * p2.vlen()) / ref;
        double ang = asin(fabs(sinang));
        return sinang > eps ? p2.rotate(-ang) : p2.rotate(ang);
}

bool cir_line(Point ct, double r, Point l1, Point l2, Point& p1, Point& p2) 
{
	if ( sgn (pldis(ct, l1, l2) - r ) > 0)
		return false;
	double a1, a2, b1, b2, A, B, C, t1, t2;
	a1 = l2.x - l1.x; a2 = l2.y - l1.y;
	b1 = l1.x - ct.x; b2 = l1.y - ct.y;
	A = a1 * a1 + a2 * a2;
	B = (a1 * b1 + a2 * b2) * 2;
	C = b1 * b1 + b2 * b2 - r * r;
	t1 = (-B - sqrt(B * B - 4.0 * A * C)) / 2.0 / A;
	t2 = (-B + sqrt(B * B - 4.0 * A * C)) / 2.0 / A;
	p1.x = l1.x + a1 * t1; p1.y = l1.y + a2 * t1;
	p2.x = l1.x + a1 * t2; p2.y = l1.y + a2 * t2;
	return true;
}
bool cir_cir(Point c1, double r1, Point c2, double r2, Point& p1, Point& p2) 
{
	double d = ppdis(c1, c2);
	if ( sgn(d - r1 - r2) > 0|| sgn (d - fabs(r1 - r2) ) < 0 )
		return false;
	Point u, v;
	double t = (1 + (r1 * r1 - r2 * r2) / ppdis(c1, c2) / ppdis(c1, c2)) / 2;
	u.x = c1.x + (c2.x - c1.x) * t;
	u.y = c1.y + (c2.y - c1.y) * t;
	v.x = u.x + c1.y - c2.y;
	v.y = u.y + c2.x - c1.x;
	cir_line(c1, r1, u, v, p1, p2);
	return true;
}
/***************************************************************************************/


/***** 三角形 **************************************************************************/
#include <math.h>
struct point{double x,y;};
struct line{point a,b;};

double distance(point p1,point p2){
	return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}

point intersection(line u,line v){
	point ret=u.a;
	double t=((u.a.x-v.a.x)*(v.a.y-v.b.y)-(u.a.y-v.a.y)*(v.a.x-v.b.x))
			/((u.a.x-u.b.x)*(v.a.y-v.b.y)-(u.a.y-u.b.y)*(v.a.x-v.b.x));
	ret.x+=(u.b.x-u.a.x)*t;
	ret.y+=(u.b.y-u.a.y)*t;
	return ret;
}

//外心
point circumcenter(point a,point b,point c){
	line u,v;
	u.a.x=(a.x+b.x)/2;
	u.a.y=(a.y+b.y)/2;
	u.b.x=u.a.x-a.y+b.y;
	u.b.y=u.a.y+a.x-b.x;
	v.a.x=(a.x+c.x)/2;
	v.a.y=(a.y+c.y)/2;
	v.b.x=v.a.x-a.y+c.y;
	v.b.y=v.a.y+a.x-c.x;
	return intersection(u,v);
}

//内心
point incenter(point a,point b,point c){
	line u,v;
	double m,n;
	u.a=a;
	m=atan2(b.y-a.y,b.x-a.x);
	n=atan2(c.y-a.y,c.x-a.x);
	u.b.x=u.a.x+cos((m+n)/2);
	u.b.y=u.a.y+sin((m+n)/2);
	v.a=b;
	m=atan2(a.y-b.y,a.x-b.x);
	n=atan2(c.y-b.y,c.x-b.x);
	v.b.x=v.a.x+cos((m+n)/2);
	v.b.y=v.a.y+sin((m+n)/2);
	return intersection(u,v);
}

//垂心
point perpencenter(point a,point b,point c){
	line u,v;
	u.a=c;
	u.b.x=u.a.x-a.y+b.y;
	u.b.y=u.a.y+a.x-b.x;
	v.a=b;
	v.b.x=v.a.x-a.y+c.y;
	v.b.y=v.a.y+a.x-c.x;
	return intersection(u,v);
}

//重心
//到三角形三顶点距离的平方和最小的点
//三角形内到三边距离之积最大的点
point barycenter(point a,point b,point c){
	line u,v;
	u.a.x=(a.x+b.x)/2;
	u.a.y=(a.y+b.y)/2;
	u.b=c;
	v.a.x=(a.x+c.x)/2;
	v.a.y=(a.y+c.y)/2;
	v.b=b;
	return intersection(u,v);
}

//费马点
//到三角形三顶点距离之和最小的点
point fermentpoint(point a,point b,point c){
	point u,v;
	double step=fabs(a.x)+fabs(a.y)+fabs(b.x)+fabs(b.y)+fabs(c.x)+fabs(c.y);
	int i,j,k;
	u.x=(a.x+b.x+c.x)/3;
	u.y=(a.y+b.y+c.y)/3;
	while (step>1e-10)
		for (k=0;k<10;step/=2,k++)
			for (i=-1;i<=1;i++)
				for (j=-1;j<=1;j++){
					v.x=u.x+step*i;
					v.y=u.y+step*j;
					if (distance(u,a)+distance(u,b)+distance(u,c)>distance(v,a)+distance(v,b)+distance(v,c))
						u=v;
				}
	return u;
}

/***************************************************************************************/


/***** 圆 ******************************************************************************/

double cir_area_inst(Point c1, double r1, Point c2, double r2)
{
        double a1, a2, d, ret;
        d = sqrt((c1 - c2).dot(c1 - c2) );
        if(sgn(d - r1 - r2) >= 0) 
                return 0;
        if(sgn(d - r2 + r1) <= 0)
                return pi * r1 * r1;
        if(sgn(d - r1 + r2) <= 0)
                return pi * r2 * r2;
        a1 = acos((r1 * r1 + d * d - r2 * r2) / 2 / r1 / d);
        a2 = acos((r2 * r2 + d * d - r1 * r1) / 2 / r2 / d);
        ret = (a1 - 0.5 * sin(2 * a1)) * r1 * r1 + (a2 - 0.5 * sin(2 * a2)) * r2 * r2;
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
Point gravity_center(Point* p, int n) {								// 多边形重心
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

bool point_in_polygon(Point o, Point *p, int n)
{
        int t;
        Point a, b;
        p[n] = p[0];
        for(int i = 0; i < n; i++) {
                if(dot_on_seg(o, Seg(p[i], p[i + 1]))) {
                        return true;
                }
        }
        t = 0;
        for(int i = 0; i < n; i++) {
                a = p[i]; b = p[i + 1];
                if(a.y > b.y) {
                        std::swap(a, b);
                }
                if(sgn((a - o) * (b - o)) < 0 && sgn(a.y - o.y) < 0 && sgn(o.y - b.y) <= 0) {
                        t++;
                }
        }
        return t & 1;
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
