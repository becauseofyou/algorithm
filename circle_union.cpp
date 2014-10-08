#include<iostream>
#include<stdio.h>
#include<string.h>
#include<cmath>
#include<ctime>
#include<algorithm>
#include<map>
using namespace std;

typedef double db;
const int maxn = 1005 ;
const db EPS = 1e-8;
const db PI = acos(- 1.0);
int sign(db x){return x < - EPS ? - 1 : x > EPS;}
db sqr(db x){return x * x;}

struct TPoint{
    db x,y;
    TPoint(){}
    TPoint(db xx,db yy):x(xx),y(yy){}
    TPoint operator+(const TPoint P){return TPoint(x + P.x, y + P.y) ;}
    TPoint operator-(const TPoint P){return TPoint(x - P.x, y - P.y) ;}
	TPoint operator/(const db k){return TPoint(x / k, y / k);}
	TPoint operator*(const db k){return TPoint(x * k, y * k);}
    bool operator<(const TPoint P)const{
        return x - P.x < - EPS || (sign(x - P.x) == 0 && y - P.y < - EPS);
    }
    db len(){return sqrt(sqr(x) + sqr(y));}
    db X(TPoint P){return x * P.y - y * P.x;}
    void get(){scanf("%lf%lf",&x,&y);}
	void out(){cout <<'('<<x<<','<<y<<')'<<endl;}
};    

struct TCircle{
    TPoint O; db r;
	db area(){return PI * sqr(r) ;}
    db get_area(db k){ return sqr(r) * (k - sin(k)) ;}
    bool operator<(const TCircle C)const{return r - C.r > EPS;}
    void get(){O.get(); scanf("%lf",&r); }
	int Rel(TCircle C){
		double d = (O-C.O).len();int d1 = sign(d - fabs(r - C.r));
		if( d1 <= 0 || !sign(d) )return 0;
		return d < r + C.r - EPS ? 1 : 2;
	}
	void get_jiao(TCircle C, TPoint &p1, TPoint &p2){
		TPoint c1 = O, c2 = C.O; db r1 = r, r2 = C.r;
		db d2 = (c1.x - c2.x) * (c1.x - c2.x) + (c1.y - c2.y) * (c1.y - c2.y);
		db co = (r1 * r1 + d2 - r2 * r2) / (2 * r1 * sqrt(d2));
		TPoint v1 = (c2 - c1) / sqrt(d2), v2 = TPoint(-v1.y, v1.x) * (r1 * sqrt(1 - co * co));
		TPoint X = c1 + v1 * (r1 * co);
		p1 = X + v2;p2 = X - v2;
	}
};    

struct Event{
	db x; 
	int typ ;
	TPoint p;
	bool operator<(const Event e)const{
		return x - e.x < - EPS || (sign(x - e.x) == 0 && typ > e.typ);
	}
	void add(db xx, int t, TPoint pp){
		x = xx; p = pp; typ = t;
	}
}s[ maxn << 1 ];

TCircle C[maxn],tC[maxn];
db area[ maxn ];

void Cir_union(TCircle C[], int &n){
	memset( area, 0, sizeof area);
    int i,j,k = 0;
    sort(C, C + n) ;
	for( i = 0; i < n; ++ i){
		int len = 0,nc = 0;
		TPoint a,b,c =TPoint(- C[i].r, 0) + C[i]. O;
		s[len++].add(- PI,1, c);s[len++].add( PI, - 1, c);
		for(j = 0; j < n; ++ j) if(i ^ j){
			int rl = C[i].Rel( C[j]);
			if(rl == 2 ) continue;
			if(! rl) {
				if( j < i) s[len++].add(- PI,1, c),s[len++].add( PI, - 1, c);
				continue;
			}
			C[i].get_jiao( C[j],a,b); 
			TPoint aa = a - C[i].O, bb = b - C[i].O;
			db jL = atan2(aa.y, aa.x), jR = atan2(bb.y, bb.x) ;
			if(jR - jL> EPS) ++ nc;
			s[len ++].add( jL, - 1, a),s[len ++].add(jR,  1, b);
		}
		sort(s, s + len);
        db pj =- PI; TPoint pp = c;
		for(j = 0; j < len; ++ j){
			db ts = C[i].get_area(s[j].x - pj) + pp.X(s[j].p);
			area[ nc ] += ts;
			if(nc > 1) area[ nc - 1] -= ts;
			nc += s[j].typ;
			pj = s[j].x,pp = s[j].p;
		}
	}
}
int n ;
bool get(){
    if(EOF == scanf("%d",&n)) return 0;
    int i ;
    for(i = 0;  i < n; ++ i) C[i].O.get(),scanf("%lf",&C[i].r);
    return 1;
}
void work(){
	Cir_union(C, n) ;
	int i;
	for(i = 1; i <= n; ++ i)  printf("[%d] = %.3f\n", i, area[i] * 0.5);
}
int main(){
    while(get()) work();
    return 0;
}


#include<cstdio>
#include<cstring>
#include<cmath>
#include<algorithm>
#include<vector>
#include<queue>
using namespace std;

struct Point
{
    double x,y;
    Point(double a=0.0,double b=0.0){x=a;y=b;}
    Point operator+(const Point&a)const{return Point(x+a.x,y+a.y);}
    Point operator-(const Point&a)const{return Point(x-a.x,y-a.y);}
    Point operator*(const double&a)const{return Point(x*a,y*a);}
    Point operator/(const double&a)const{return Point(x/a,y/a);}
    double operator*(const Point&a)const{return x*a.y-y*a.x;}
    double operator/(const Point&a)const{return sqrt((a.x-x)*(a.x-x)+(a.y-y)*(a.y-y));}
    double operator%(const Point&a)const{return x*a.x+y*a.y;}
}po[1005];
double r[1005];
const double eps = 1e-7;
const double pi=acos(-1.0);
inline int sgn(double x)
{return fabs(x)<eps?0:(x>0.0?1:-1);}
pair<double,bool>arg[2005];
double cir_union(Point c[],double r[],int n)
{
    double sum=0.0,sum1=0.0,d,p1,p2,p3;
    for(int i=0;i<n;i++)
    {
        bool f=1;
        for(int j=0;f&&j<n;j++)
            if(i!=j&&sgn(r[j]-r[i]-c[i]/c[j])!=-1)f=0;
        if(!f)swap(r[i],r[--n]),swap(c[i--],c[n]);
    }
    for(int i=0;i<n;i++)
    {
        int k=0,cnt=0;
        for(int j=0;j<n;j++)
            if(i!=j&&sgn((d=c[i]/c[j])-r[i]-r[j])<=0)
            {
                p3=acos((r[i]*r[i]+d*d-r[j]*r[j])/(2.0*r[i]*d));
                p2=atan2(c[j].y-c[i].y,c[j].x-c[i].x);
                p1=p2-p3;p2=p2+p3;
                if(sgn(p1+pi)==-1)p1+=2*pi,cnt++;
                if(sgn(p2-pi)==1)p2-=2*pi,cnt++;
                arg[k++]=make_pair(p1,0);arg[k++]=make_pair(p2,1);
            }
        if(k)
        {
            sort(arg,arg+k);
            p1=arg[k-1].first-2*pi;
            p3=r[i]*r[i];
            for(int j=0;j<k;j++)
            {
                p2=arg[j].first;
                if(cnt==0)
                {
                    sum+=(p2-p1-sin(p2-p1))*p3;
                    sum1+=(c[i]+Point(cos(p1),sin(p1))*r[i])*(c[i]+Point(cos(p2),sin(p2))*r[i]);
                }
                p1=p2;
                arg[j].second?cnt--:cnt++;
            }
        }
        else sum+=2*pi*r[i]*r[i];
    }
    return (sum+fabs(sum1))*0.5;
}
int main()
{
    int n;scanf("%d",&n);
    for(int i=0;i<n;i++)
        scanf("%lf%lf%lf",&po[i].x,&po[i].y,r+i);
    printf("%.3f\n",cir_union(po,r,n));
    return 0;
}
