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
        inline Vector norm() {
                return Point(-y, x);
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
        bool in_convex(Point pt) {
                if(sgn((p[1]-p[0])*(pt-p[0])) <= 0 || sgn((p[n-1]-p[0])*(pt-p[0])) >= 0) {
                        return false;
                }
                int l = 1, r = n - 2, best = -1;
                while(l <= r) {
                        int mid = l + r >> 1;
                        int f = sgn((p[mid]-p[0])*(pt-p[0]));
                        if(f >= 0) {
                                best = mid;
                                l = mid + 1;
                        } else {
                                r = mid - 1;
                        }
                }
                return     sgn((p[best+1]-p[best])*(pt-p[best])) > 0;
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
struct Point {
    double x, y;
	Point(){}
	Point (double tx,double ty)
	{
		this->x = tx;
		this->y = ty;
	}
	bool operator == (const Point &t) const {
		return t.x == x && t.y == y;
	}
	Point operator - (const Point &t) const {
		Point res;
		res.x = x - t.x;
		res.y = y - t.y;
		return res;
	}
}; 
double multi(Point &o, Point &a, Point &b) {							// 点积
	return (a.x-o.x)*(b.x-o.x) + (a.y-o.y)*(b.y-o.y);
}
double cross(Point &o, Point &a, Point &b) {							// 叉积
	return (a.x-o.x)*(b.y-o.y) - (b.x-o.x)*(a.y-o.y); 
}
double cp(Point &a, Point &b) {									
	return a.x*b.y - b.x*a.y;
}
double angle(Point &a, Point &b) {												// 两向量夹角
	double ans = fabs((atan2(a.y, a.x) - atan2(b.y, b.x)));
	return ans > pi+eps ? 2*pi-ans : ans;
}
double cir_polygon(Point ct, double R, Point *p, int n) {					// 圆与简单多边形
	Point o, a, b, t1, t2;
	double sum=0, res, d1, d2, d3, sign;
	o.x = o.y = 0;
	p[n] = p[0];
	for (int i=0; i < n; i++) {
		a = p[i]-ct;
		b = p[i+1]-ct;
		sign = cp(a,b) > 0 ? 1 : -1;
		d1 = a.x*a.x + a.y*a.y;
		d2 = b.x*b.x + b.y*b.y;
		d3 = sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
		if (d1 < R*R+eps && d2 < R*R+eps) { //两个点都在圆内
			res = fabs(cp(a, b));
		}
		else if (d1 < R*R-eps || d2 < R*R-eps) { //一个点在圆内
			cir_line(o, R, a, b, t1, t2);
			if ((a.x-t2.x)*(b.x-t2.x) < eps && (a.y-t2.y)*(b.y-t2.y) < eps) {
				t1 = t2;
			}
			if (d1 < d2) 
				res = fabs(cp(a, t1)) + R*R*angle(b, t1);
			else
				res = fabs(cp(b, t1)) + R*R*angle(a, t1);
		}
		else if (fabs(cp(a, b))/d3 > R-eps) { // 两个点都在园外，且线段与圆之多只有一个交点
			res = R*R*angle(a, b);
		}
		else {  // 线段与圆有两个交点
			cir_line(o, R, a, b, t1, t2); 
			if (multi(t1, a, b) > eps || multi(t2, a, b) > eps) { // a b 在圆的同一侧
				res = R*R*angle(a, b);
			}
			else {
				res = fabs(cp(t1, t2));
				if (cross(t1, t2, a) < eps) 
					res += R*R*(angle(a, t1) + angle(b, t2));
				else 
					res += R*R*(angle(a, t2) + angle(b, t1));
			}
		}			
		sum += res * sign;
	}
	return fabs(sum)/2.0;
}
Point p[55];
int main()
{
    int t,ca=1;
	double x1,x2,x3,x4,y1,y2,y3,y4,R;
	while(scanf("%lf",&R)!=EOF)
	{
		Point cen = Point(0,0);
		int n;
		scanf("%d",&n);
		for(int i = 0; i < n; i++)
		{
			scanf("%lf%lf",&p[i].x,&p[i].y);
		}
		printf("%.2lf\n",cir_polygon(cen,R,p,n));
	}
	return 0;
}
//判断n+1个圆是否有交,R[n]=mid
//nlogn判断n个圆是否有交是这样的，我们先锁定x的范围，也就是所有圆的右边界的最小值right，
//以及所有圆的左边界的最大值left  那么n个圆的公共部分肯定在left right之间，
//而且有一个很重要的性质，如果n个圆的公共部分的左右区间是L,R,假设我们当前枚举的答案是mid，
//如果x=mid这条直线与所有的圆没有公共部分，
//那么我们找两个圆与直线的公共部分不相交（其实就是上边界最小以及下边界最大的两个圆），
//判断这两个圆的公共部分在x=mid的哪一侧即可，其实就是满足二分性质
bool judge(double mid)  
{  
    double Left,Right;  
    for(int i = 0; i <= n; i++) {  
        if(i == 0) {  
            Left = p[i].x - R[i];  
            Right = p[i].x + R[i];  
        } else{  
            if(p[i].x-R[i] > Left) Left = p[i].x-R[i];  
            if(p[i].x+R[i] < Right) Right = p[i].x+R[i];  
        }  
    }  
  
    if(Left - Right > eps) return false;  
    int step = 50;  
    while(step--) {  
        double mid = (Left + Right)*0.5;  
        double low,high,uy,dy;  
        int low_id,high_id;  
        for(int i = 0; i <= n; i++) {  
            double d = sqrt(R[i]*R[i]-(p[i].x-mid)*(p[i].x-mid));  
            uy = p[i].y + d;  
            dy = p[i].y - d;  
            if(i == 0) {  
                low_id = high_id = 0;  
                low = dy; high = uy;  
            } else {  
                if(uy < high) high = uy,high_id = i;  
                if(dy > low) low = dy,low_id = i;  
            }  
        }  
          
        if(high - low > -eps) {  
            return 1;  
        }  
        Point a,b;  
        if(Cir_Cir(p[high_id],R[high_id],p[low_id],R[low_id],a,b)) {  
            if((a.x+b.x)*0.5 < mid) {  
                Right = mid;  
            } else Left = mid;  
        } else return false;  
    }  
    return false;  
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

//半平面交
Point intersect(Seg s1, Seg s2)
{
        return intersect(s1.s, s1.e-s1.s, s2.s, s2.e-s2.s);
}

int ord[N], dq[N];
Seg sg[N];  
double at2[N];

int cmp(int a, int b)
{
        if(sgn(at2[a]-at2[b]) != 0) {
                return sgn(at2[a] - at2[b]) < 0;
        }
        return sgn((sg[a].e-sg[a].s)*(sg[b].e-sg[a].s)) < 0;
}
bool is_right(int a, int b, int c) 
{
        Point t;
        t = intersect(sg[a], sg[b]);
        return sgn((sg[c].e-sg[c].s)*(t-sg[c].s)) < 0;
}
int HPI(int n, std::vector<Point>&p) {
        int i, j, l=1, r=2;
        for(i = 0; i < n; i++) {
                at2[i] = atan2(sg[i].e.y-sg[i].s.y, sg[i].e.x-sg[i].s.x);
                ord[i] = i;
        }
        std::sort(ord, ord + n, cmp);
        for(i=j=1; i < n; i++) if(sgn(at2[ord[i]]-at2[ord[i-1]]) > 0) {
                ord[j++] = ord[i];
        }
        n = j;
        p.clear();
        dq[l] = ord[0];dq[r] = ord[1];
        for(i=2; i < n; i++) {
                for(; l < r && is_right(dq[r-1],dq[r],ord[i]); r--) {
                        if(sgn(at2[ord[i]] - at2[dq[r-1]] - pi) >= 0) 
                                return -1;
                }
                while(l < r && is_right(dq[l], dq[l+1], ord[i])) l++;
                dq[++r] = ord[i];
        }
        while(l < r && is_right(dq[r-1], dq[r], dq[l])) r--;
        while(l < r && is_right(dq[l], dq[l+1], dq[r])) l++;
        dq[l-1] = dq[r]; p.resize(r-l+1);
        for(int i = l; i <= r; i++) {
                p[i-l]=intersect(sg[dq[i]], sg[dq[i-1]]);
        }
        return r - l + 1;
}


//三维几何////////////////////////////////////////////////////////////
struct Point;
typedef double DB;
typedef Point PT;
#define op operator
inline int sgn(DB x) {
         return x < -eps ? -1 : x > eps;
}
struct Point {
        DB x, y, z;
        Point(DB x=0,DB y=0, DB z=0): x(x), y(y), z(z) {}
        void in() {
                scanf("%lf%lf%lf", &x, &y, &z);
        }
        void print() {
                printf("%lf %lf %lf\n", x, y, z);
        }
        inline PT op * (const PT &t) const {
                return PT(y * t.z - t.y * z, z * t.x - t.z * x, x * t.y - t.x * y);
        }
        inline PT op - (const PT &t) const {
                return PT(x - t.x, y - t.y, z - t.z);
        }
        inline PT op + (const PT &t) const {
                return PT(x + t.x, y + t.y, z + t.z);
        }
        inline PT op / (DB d) {
                return PT(x / d, y / d, z / d);
        }
        inline PT op * (DB d) {
                return PT(x * d, y * d, z * d);
        }
        inline bool op == (const PT &t) const {
                return sgn(x - t.x) == 0 && sgn(y - t.y) == 0 && sgn(z - t.z) == 0;
        }
        inline bool op < (const PT &t) const {
                int f1=sgn(x-t.x), f2=sgn(y-t.y), f3=sgn(z-t.z);
                return f1<0||(f1==0&&f2<0)||(f1==0&&f2==0&&f3<0);
        }
        inline DB vlen() {
                return sqrt(x * x + y * y + z * z);
        }
        inline DB dot(const PT& t) const {
                return x * t.x + y * t.y + z * t.z;
        }
};
struct line{PT a, b;};
struct plane{ PT a,b,c; };
PT norm(PT a, PT b, PT c) {
        return (b - a) * (c - a);
}
PT norm(plane s) {
        return (s.b - s.a) * (s.c - s.a);
}
PT corplane(PT a, PT b, PT c, PT d) {
        return sgn(norm(a, b, c).dot(d - a)) == 0;
}
bool same_side(PT p1, PT p2, PT a, PT b, PT c) {
        PT v = norm(a, b, c);
        return sgn(v.dot(p1-a)) * sgn(v.dot(p2-a)) >= 0;
}
double ptoplane(PT p, plane s) {
        return fabs(norm(s).dot(p - s.a)) / norm(s).vlen();
}
double ptoplane(PT p, PT s1, PT s2, PT s3) {
        return fabs(norm(s1, s2, s3).dot(p - s1)) / norm(s1, s2, s3).vlen();
}
double area(PT a, PT b, PT c) {
        return ((b-a)*(c-a)).vlen();
}
//四面体面积*6
double volume(PT a, PT b, PT c, PT d) {
        return ((b-a)*(c-a)).dot(d-a);
}
inline Point get_point(Point st,Point ed,Point tp){
        double t1=(tp-st).dot(ed-st);
        double t2=(ed-st).dot(ed-st);
        double t=t1/t2;
        Point ans=st + ((ed-st)*t);
        return ans;
}
inline double dist(Point st,Point ed) {
        return sqrt((ed-st).dot(ed-st));
}
//沿着直线st-ed旋转角度A，从ed往st看是逆时针
Point rotate(Point st,Point ed,Point tp,double A){
        Point root=get_point(st,ed,tp);
        Point e=(ed-st)/dist(ed,st);
        Point r=tp-root;
        Point vec=e*r;
        Point ans=r*cos(A)+vec*sin(A)+root;
        return ans;
}

struct face { 
        int a,b,c;
        bool ok;
};
struct CH3D {
        static const int MAXN = 55;
        int n;//初始顶点数
        PT P[MAXN];//初始顶点
        int num;//凸包表面的三角形个数
        face F[8*MAXN];//凸包表面的三角形
        int g[MAXN][MAXN];
        double vol(PT &p, face &f) {
                Point m=P[f.b]-P[f.a];
                Point n=P[f.c]-P[f.a];
                Point t=p-P[f.a];
                return (m*n).dot(t);
        }
        void deal(int p,int a,int b) {
                int f=g[a][b];
                face add;
                if(F[f].ok) {
                        if(sgn(vol(P[p],F[f])) > 0) {
                                dfs(p,f);
                        } else {
                                add.a=b; add.b=a; add.c=p; add.ok=true;
                                g[p][b]=g[a][p]=g[b][a]=num;
                                F[num++]=add;
                        }
                }
        }
        void dfs(int p,int now) {
                F[now].ok=false;
                deal(p,F[now].b,F[now].a);
                deal(p,F[now].c,F[now].b);
                deal(p,F[now].a,F[now].c);
        }
        bool same(int s,int t) {
                Point a=P[F[s].a];
                Point b=P[F[s].b];
                Point c=P[F[s].c];
                return sgn(volume(a,b,c,P[F[t].a])) == 0 &&
                        sgn(volume(a,b,c,P[F[t].b])) == 0 &&
                        sgn(volume(a,b,c,P[F[t].c])) == 0;
        }
        void create() {
                int i,j,tmp;
                face add;
                num=0;
                if(n<4)return;
                bool flag=true;
                for(i=1;i<n;i++) {
                        if(sgn((P[0]-P[i]).vlen()) > 0) {
                                std::swap(P[1],P[i]);
                                flag=false;
                                break;
                        }
                }
                if(flag)return;
                flag=true;
                for(i=2;i<n;i++) {
                        if(sgn(area(P[0], P[1], P[i])) > 0) {
                                std::swap(P[2],P[i]);
                                flag=false;
                                break;
                        }
                }
                if(flag)return;
                flag=true;
                for(i=3;i<n;i++) {
                        if(sgn(fabs(volume(P[0], P[1], P[2], P[i]))) > 0) {
                                std::swap(P[3],P[i]);
                                flag=false;
                                break;
                        }
                }
                if(flag)return;
                for(i=0;i<4;i++) {
                        add.a=(i+1)%4;
                        add.b=(i+2)%4;
                        add.c=(i+3)%4;
                        add.ok=true;
                        if(sgn(vol(P[i], add)) > 0) {
                                std::swap(add.b,add.c);
                        }
                        g[add.a][add.b]=g[add.b][add.c]=g[add.c][add.a]=num;
                        F[num++]=add;
                }
                for(i=4;i<n;i++) {
                        for(j=0;j<num;j++) {
                                if(F[j].ok && sgn(vol(P[i],F[j])) > 0) {
                                        dfs(i,j);
                                        break;
                                }
                        }
                }
                tmp=num;
                for(i=num=0;i<tmp;i++)
                        if(F[i].ok)
                                F[num++]=F[i];
        }
        double ptoface(PT p,int i) {
                return fabs(volume(P[F[i].a],P[F[i].b],P[F[i].c],p)/((P[F[i].b]-P[F[i].a])*(P[F[i].c]-P[F[i].a])).vlen());
        }
}hull;
//旋转点集，使法向量为v的平面水平
void rotate_to_horizontal(int n, PT *p, PT v) {
        if(sgn(v.x)==0 && sgn(v.y)==0) {
                return ;
        }
        double a, c, s;
        a = atan2(v.y, v.x);
        c = cos(a), s = sin(a);
        for(int i = 0; i < n; i++) {
                PT t = p[i];
                p[i].x = t.x * c + t.y * s;
                p[i].y = t.y * c - t.x * s;
        }
        a = atan2(sqrt(v.x*v.x+v.y*v.y), v.z);
        c = cos(a), s = sin(a);
        for(int i = 0; i < n; i++) {
                PT t = p[i];
                p[i].z = t.z * c + t.x * s;
                p[i].x = t.x * c - t.z * s;
        }
}
struct Convex_Hull {
        static const int N = 110;
        Point p[2 * N];
        int n;
        void init() {
                n = 0;
        }
        inline void push_back(const Point &np) {
                p[n++] = np;
        }
        DB cross(PT a, PT b) {
                return a.x * b.y - a.y * b.x;
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
                        while(top >= m + 2 && sgn(cross(p[top-1]-p[top-2],p[i]-p[top-2])) <= 0) {
                                top --;
                        }
                        p[top++] = p[i];
                        if(i == n - 1) {
                                m = top - 1;
                        }
                }
                n = top - 1;
        }
        double get_area() {
                double ret = 0;
                for(int i = 0; i < n; i++) {
                        ret += p[i].x * (p[i + 1].y) - p[i].y*p[i+1].x;
                }
                return fabs(ret) / 2;
        }
}convex;
PT p[55]; int n;
bool input()
{
        scanf("%d", &n);
        if(n == 0) return false;
        hull.n = n;
        for(int i = 0; i < n; i++) {
                hull.P[i].in();
        }
        double ans_h = 0, ans_area = 1e30;
        hull.create();
        for(int i = 0; i < hull.num; i++) {
                for(int j = 0; j < n; j++) {
                        p[j] = hull.P[j];
                }
                PT v = norm(p[hull.F[i].b], p[hull.F[i].a], p[hull.F[i].c]);
                rotate_to_horizontal(n, p, v);
                double z = p[hull.F[i].a].z;
                for(int j = 0; j < n; j++) {
                        p[j].z -= z;
                }
                double H = 0;
                for(int j = 0; j < n; j++) {
                        H = std::max(H, hull.ptoface(hull.P[j], i));
                }
                convex.init();
                for(int j = 0; j < n; j++) {
                        convex.push_back(p[j]);
                }
                convex.gao();
                double S = convex.get_area();
                if(sgn(H - ans_h) > 0 || sgn(H - ans_h)==0 && sgn(ans_area-S) > 0) {
                        ans_h = H;
                        ans_area = S;
                }
        }
        printf("%.3f %.3f\n", ans_h, ans_area);
        return true;
}
