// vimrc
set guifont=Monaco:h08:b
set cin nu rnu ts=4 sw=4 sts=4 et noswapfile nobackup
set so=100
set backspace=eol,start,indent
"colorscheme last256
syntax on
map <F4> :!g++ -std=c++11 %<.cpp -m32 -Wall -o %<<cr>
map <F5> :!%< < %<.in<cr>
map <F6> :vsplit %<.in<cr>

typedef double DB;
#define op operator 
const DB eps = 1e-8;
inline int sgn(DB x) {  return x < -eps ? -1 : x > eps; }
struct PT {
    DB x, y;
    PT (DB x=0, DB y = 0) : x(x), y(y){ }
    PT norm() { return PT(-y, x); }
    PT rotate(DB ang) { return PT(x * cos(ang) - y * sin(ang), x * sin(ang) + y * cos(ang)); }
};
//矢量V以P为顶点逆时针旋转angle并放大scale倍
point rotate(point v,point p,double angle,double scale){
	point ret=p;
	v.x-=p.x,v.y-=p.y;
	p.x=scale*cos(angle);
	p.y=scale*sin(angle);
	ret.x+=v.x*p.x-v.y*p.y;
	ret.y+=v.x*p.y+v.y*p.x;
	return ret;
}
struct Seg { PT s, e; }
struct Line{ int a, b, c; };
struct Cir{PT ct; DB r;}
bool dot_on_seg(PT p, Seg L) { return sgn((L.s - p) * (L.e - p)) == 0 && sgn((L.s - p).dot(L.e - p)) <= 0; }
DB ppdis(PT a, PT b) { return sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) ); }
bool intersect(PT P, PT v, PT Q, PT w, PT &p) {
    PT u = P - Q;
    if(sgn(v * w) == 0) return false;
    double t = w * u / (v * w);
    p = P + v * t;
    return true;
}
double cross(PT a, PT b, PT c){ return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);}
PT disptoline(PT p, Seg l) { return fabs(cross(p, l.s, l.e)) / ppdis(l.s, l.e); }
PT ptoline(PT p, Seg l) { PT vec = l.s - l.e; return intersect(p, vec.norm(), l.s, vec); }
PT ptoseg(PT p, Seg l) { 
    PT norm = (l.s - l.e).norm();
    if(sgn(norm * (p - l.s)) * sgn(norm * (p - l.e)) > 0) {
        double sa = ppdis(p, l.s);
        double sb = ppdis(p, l.e);
        return sgn(sa - sb) < 0 ? l.s : l.e;
    }
    return intersect(p, norm, l.s, l.e - l.s);
}
bool point_in_polygon(PT o, PT *p, int n, bool flag) {  //传入flag表示在边界上算不算在里面
    int t = 0; PT a, b;
    for(int i = 0; i < n; i++) if(dot_on_seg(o, Seg(p[i], p[(i + 1) % n]))) return flag;
    for(int i = 0; i < n; i++) {
        a = p[i]; b = p[(i + 1) % n];
        if(a.y > b.y) std::swap(a, b);
        if(sgn((a - o) * (b - o)) < 0 && sgn(a.y - o.y) < 0 && sgn(o.y - b.y) <= 0) 
            t++;
    }
    return t & 1;
}
bool segcross(PT p1, PT p2, PT q1, PT q2)  { // 快速判断线段相交
    return (
            std::min(p1.x, p2.x) <= std::max(q1.x, q2.x) &&
            std::min(q1.x, q2.x) <= std::max(p1.x, p2.x) &&
            std::min(p1.y, p2.y) <= std::max(q1.y, q2.y) &&
            std::min(q1.y, q2.y) <= std::max(p1.y, p2.y) && 
            cross(p1, q2, q1) * cross(p2, q2, q1) <= 0 && 
            cross(q1, p2, p1) * cross(q2, p2, p1) <= 0  
           ); 
}
void convex(PT *p, int &n) {
    if(n < 3) { return ; }
    std::sort(p, p + n);
    std::copy(p, p + n - 1, p + n);
    std::reverse(p + n, p + 2 * n - 1);
    int m = 0, top = 0;
    for(int i = 0; i < 2 * n - 1; i++) {
        while(top >= m + 2 && sgn((p[top - 1] - p[top - 2]) * (p[i] - p[top - 2])) <= 0) { top --; }
        p[top++] = p[i];
        if(i == n - 1) { m = top - 1; }
    }
    n = top - 1;
}
int inhalfplane(point p,point s,point e) { return sgn(cross(e - s, p - s)) ; }
std::vector<point> CUT(const std::vector<point> &p, point s, point e) {
    std::vector<point> q;
    int n = (int) p.size();
    for(int i = 0; i < n; i++) {
        int nowin = inhalfplane(p[i], s, e);
        int nextin = inhalfplane(p[(i + 1) % n], s, e);
        if(nowin >= 0) 
            q.push_back(p[i]);
        if(nextin * nowin < 0) 
            q.push_back(intersect(p[i], p[(i + 1) % n] - p[i], s, e - s));
    }
    return q;
}
bool in_convex(PT *p, PT pt, int n) {
    if(sgn((p[1]-p[0])*(pt-p[0])) <= 0 || sgn((p[n-1]-p[0])*(pt-p[0])) >= 0) { return false; }
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
//圆的折射，输入圆心，p->inter是射线，inter是圆上一点， ref是折射率，返回折射向量
Vector reflect_vector(PT center, PT p, PT inter, double ref) {
    Vector p1 = inter - p, p2 = center - inter;
    double sinang = p1 * p2 / (p1.vlen() * p2.vlen()) / ref;
    double ang = asin(fabs(sinang));
    return sinang > eps ? p2.rotate(-ang) : p2.rotate(ang);
}
bool cir_line(PT ct, double r, PT l1, PT l2, PT& p1, PT& p2) {
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
bool cir_cir(PT c1, double r1, PT c2, double r2, PT& p1, PT& p2) {
    double d = ppdis(c1, c2);
    if ( sgn(d - r1 - r2) > 0|| sgn (d - fabs(r1 - r2) ) < 0 )
        return false;
    PT u, v;
    double t = (1 + (r1 * r1 - r2 * r2) / ppdis(c1, c2) / ppdis(c1, c2)) / 2;
    u.x = c1.x + (c2.x - c1.x) * t;
    u.y = c1.y + (c2.y - c1.y) * t;
    v.x = u.x + c1.y - c2.y;
    v.y = u.y + c2.x - c1.x;
    cir_line(c1, r1, u, v, p1, p2);
    return true;
}
struct Point // 求圆并
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
double cir_union(Point c[],double r[],int n) {
    double sum=0.0,sum1=0.0,d,p1,p2,p3;
    for(int i=0;i<n;i++) {
        bool f=1;
        for(int j=0;f&&j<n;j++)
            if(i!=j&&sgn(r[j]-r[i]-c[i]/c[j])!=-1)f=0;
        if(!f)swap(r[i],r[--n]),swap(c[i--],c[n]);
    }
    for(int i=0;i<n;i++) {
        int k=0,cnt=0;
        for(int j=0;j<n;j++)
            if(i!=j&&sgn((d=c[i]/c[j])-r[i]-r[j])<=0) {
                p3=acos((r[i]*r[i]+d*d-r[j]*r[j])/(2.0*r[i]*d));
                p2=atan2(c[j].y-c[i].y,c[j].x-c[i].x);
                p1=p2-p3;p2=p2+p3;
                if(sgn(p1+pi)==-1)p1+=2*pi,cnt++;
                if(sgn(p2-pi)==1)p2-=2*pi,cnt++;
                arg[k++]=make_pair(p1,0);arg[k++]=make_pair(p2,1);
            }
        if(k) {
            sort(arg,arg+k);
            p1=arg[k-1].first-2*pi;
            p3=r[i]*r[i];
            for(int j=0;j<k;j++) {
                p2=arg[j].first;
                if(cnt==0) {
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
double multi(PT &o, PT &a, PT &b) {	return (a.x-o.x)*(b.x-o.x) + (a.y-o.y)*(b.y-o.y); }
double cross(PT &o, PT &a, PT &b) {return (a.x-o.x)*(b.y-o.y) - (b.x-o.x)*(a.y-o.y); }
double cp(PT &a, PT &b) {return a.x*b.y - b.x*a.y; }
double angle(PT &a, PT &b) {												// 两向量夹角
    double ans = fabs((atan2(a.y, a.x) - atan2(b.y, b.x)));
    return ans > pi+eps ? 2*pi-ans : ans;
}
double cir_polygon(PT ct, double R, PT *p, int n) {					// 圆与简单多边形
    PT o, a, b, t1, t2;
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
double polyUnion( int n ) { //多边形面积并
    double sum = 0;
    for( int i = 0; i < n; ++i ) for( int ii = 0; ii < g[i].sz; ++ii ) {
        int tot = 0;
        c[tot++] = MP(0, 0);
        c[tot++] = MP(1, 0);
        for( int j = 0; j < n; ++j ) if( i != j ) for( int jj = 0; jj < g[j].sz; ++jj ) {
            int d1 = dcmp(cross(g[i].p[ii+1] - g[i].p[ii], g[j].p[jj] - g[i].p[ii]));
            int d2 = dcmp(cross(g[i].p[ii+1] - g[i].p[ii], g[j].p[jj+1] - g[i].p[ii]));
            if( !d1 && !d2 ) {
                point t1 = g[j].p[jj+1] - g[j].p[jj];
                point t2 = g[i].p[ii+1] - g[i].p[ii];
                if( dcmp( dot(t1, t2) ) > 0 && j < i ) {
                    c[tot++] = MP(segP(g[j].p[jj], g[i].p[ii], g[i].p[ii+1]), 1);
                    c[tot++] = MP(segP(g[j].p[jj+1], g[i].p[ii], g[i].p[ii+1]), -1);
                }
            }
            else if( d1 >= 0 && d2 < 0 || d1 < 0 && d2 >= 0 ) {
                double tc = cross(g[j].p[jj+1] - g[j].p[jj], g[i].p[ii] - g[j].p[jj]);
                double td = cross(g[j].p[jj+1] - g[j].p[jj], g[i].p[ii+1] - g[j].p[jj]);
                if( d2 < 0 )
                    c[tot++] = MP(tc / (tc - td), 1);
                else c[tot++] = MP(tc / (tc - td), -1);
            }
        }
        sort(c, c + tot);
        double cur = min(max(c[0].first, 0.0), 1.0);
        int sgn = c[0].second;
        double s = 0;
        for( int j = 1; j < tot; ++j ) {
            double nxt = min(max(c[j].first, 0.0), 1.0);
            if( !sgn ) s += nxt - cur;
            sgn += c[j].second;
            cur = nxt;
        }
        sum += cross(g[i].p[ii], g[i].p[ii+1]) * s;
    }
    return sum / 2;
}
bool judge(double mid)  {  //判断n+1个圆是否有交,R[n]=mid
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

        if(high - low > -eps) {  return 1;  }  
        PT a,b;  
        if(Cir_Cir(p[high_id],R[high_id],p[low_id],R[low_id],a,b)) {  
            if((a.x+b.x)*0.5 < mid) {  
                Right = mid;  
            } else Left = mid;  
        } else return false;  
    }  
    return false;  
}  
struct point{double x,y;};
struct line{point a,b;};
double distance(point p1,point p2){ return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)); }
point intersection(line u,line v){
    point ret=u.a;
    double t=((u.a.x-v.a.x)*(v.a.y-v.b.y)-(u.a.y-v.a.y)*(v.a.x-v.b.x))
        /((u.a.x-u.b.x)*(v.a.y-v.b.y)-(u.a.y-u.b.y)*(v.a.x-v.b.x));
    ret.x+=(u.b.x-u.a.x)*t;
    ret.y+=(u.b.y-u.a.y)*t;
    return ret;
}
point circumcenter(point a,point b,point c){//外心
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
point incenter(point a,point b,point c){//内心
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
point perpencenter(point a,point b,point c){//垂心
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
double cir_area_inst(PT c1, double r1, PT c2, double r2) {
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
PT gravity_center(PT* p, int n) {								// 多边形重心
    int i;
    double A=0, a;
    PT t;
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
bool cmpag(PT a, PT b) {
    double t = (a-GP).x*(b-GP).y - (b-GP).x*(a-GP).y;
    return fabs(t) > eps ? t > 0 : PPdis(a, GP) < PPdis(b, GP);
}
void Grahamag(PT *p, int &n) {									// 极角序                    
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
/**半平面交***************************/
int ord[N], dq[N];
Seg sg[N];  
double at2[N];
int cmp(int a, int b) {
    if(sgn(at2[a]-at2[b]) != 0) {
        return sgn(at2[a] - at2[b]) < 0;
    }
    return sgn((sg[a].e-sg[a].s)*(sg[b].e-sg[a].s)) < 0;
}
bool is_right(int a, int b, int c) {
    PT t;
    t = intersect(sg[a], sg[b]);
    return sgn((sg[c].e-sg[c].s)*(t-sg[c].s)) < 0;
}
PT intersect(Seg s1, Seg s2) { return intersect(s1.s, s1.e-s1.s, s2.s, s2.e-s2.s); }
int HPI(int n, std::vector<PT>&p) {
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
//返回两点所在大圆劣弧对应圆心角,0<=angle<=pi
double angle(double lng1,double lat1,double lng2,double lat2){//计算圆心角lat表示纬度,-90<=w<=90,lng表示经度
    double dlng=fabs(lng1-lng2)*pi/180;
    while (dlng>=pi+pi)
        dlng-=pi+pi;
    if (dlng>pi)
        dlng=pi+pi-dlng;
    lat1*=pi/180,lat2*=pi/180;
    return acos(cos(lat1)*cos(lat2)*cos(dlng)+sin(lat1)*sin(lat2));
}
double line_dist(double r,double lng1,double lat1,double lng2,double lat2){//计算距离,r为球半径
    double dlng=fabs(lng1-lng2)*pi/180;
    while (dlng>=pi+pi)
        dlng-=pi+pi;
    if (dlng>pi)
        dlng=pi+pi-dlng;
    lat1*=pi/180,lat2*=pi/180;
    return r*sqrt(2-2*(cos(lat1)*cos(lat2)*cos(dlng)+sin(lat1)*sin(lat2)));
}
inline double sphere_dist(double r,double lng1,double lat1,double lng2,double lat2){//计算球面距离,r为球半径
    return r*angle(lng1,lat1,lng2,lat2);
}

#define zero(x) (((x)>0?(x):-(x))<eps)
struct point3{double x,y,z;};
struct line3{point3 a,b;};
struct plane3{point3 a,b,c;};
point3 xmult(point3 u,point3 v){//计算cross product U x V
    point3 ret;
    ret.x=u.y*v.z-v.y*u.z;
    ret.y=u.z*v.x-u.x*v.z;
    ret.z=u.x*v.y-u.y*v.x;
    return ret;
}
double dmult(point3 u,point3 v){//计算dot product U . V
    return u.x*v.x+u.y*v.y+u.z*v.z;
}
point3 subt(point3 u,point3 v){//矢量差 U - V
    point3 ret;
    ret.x=u.x-v.x;
    ret.y=u.y-v.y;
    ret.z=u.z-v.z;
    return ret;
}
double volume(point3 a, point3 b, point3 c, point3 d) {//四面体体积
    return fabs(dmult( (b - a) * (c - a) , (d - a) ) ) / 6.0;
}
point3 shade_ptoplane(point3 p, point3 a, point3 b, point3 c) {//点到平面的投影
    point3 nor = (b - a) * (c - a);
    point3 nor0 = ( nor * dmult(nor, p - a) ) / vlen(nor) / vlen(nor);
    return (p - nor0);
}
point3 pvec(point3 s1,point3 s2,point3 s3){//取平面法向量
    return xmult(subt(s1,s2),subt(s2,s3));
}
double distance(point3 p1,point3 p2){//两点距离,单参数取向量大小
    return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z));
}
double vlen(point3 p){//向量大小
    return sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
}
int dots_inline(point3 p1,point3 p2,point3 p3){//判三点共线
    return vlen(xmult(subt(p1,p2),subt(p2,p3)))<eps;
}
int dots_onplane(point3 a,point3 b,point3 c,point3 d){//判四点共面
    return zero(dmult(pvec(a,b,c),subt(d,a)));
}
int dot_online_in(point3 p,point3 l1,point3 l2){//判点是否在线段上,包括端点和共线
    return zero(vlen(xmult(subt(p,l1),subt(p,l2))))&&(l1.x-p.x)*(l2.x-p.x)<eps&&
        (l1.y-p.y)*(l2.y-p.y)<eps&&(l1.z-p.z)*(l2.z-p.z)<eps;
}
int dot_online_ex(point3 p,point3 l1,point3 l2){//判点是否在线段上,不包括端点
    return dot_online_in(p,l1,l2)&&(!zero(p.x-l1.x)||!zero(p.y-l1.y)||!zero(p.z-l1.z))&&
        (!zero(p.x-l2.x)||!zero(p.y-l2.y)||!zero(p.z-l2.z));
}
int dot_inplane_in(point3 p,point3 s1,point3 s2,point3 s3){//判点是否在空间三角形上,包括边界,三点共线无意义
    return zero(vlen(xmult(subt(s1,s2),subt(s1,s3)))-vlen(xmult(subt(p,s1),subt(p,s2)))-
            vlen(xmult(subt(p,s2),subt(p,s3)))-vlen(xmult(subt(p,s3),subt(p,s1))));
}
int dot_inplane_ex(point3 p,point3 s1,point3 s2,point3 s3){//判点是否在空间三角形上,不包括边界,三点共线无意义
    return dot_inplane_in(p,s1,s2,s3)&&vlen(xmult(subt(p,s1),subt(p,s2)))>eps&&
        vlen(xmult(subt(p,s2),subt(p,s3)))>eps&&vlen(xmult(subt(p,s3),subt(p,s1)))>eps;
}
int same_side(point3 p1,point3 p2,point3 l1,point3 l2){//判两点在线段同侧,点在线段上返回0,不共面无意义
    return dmult(xmult(subt(l1,l2),subt(p1,l2)),xmult(subt(l1,l2),subt(p2,l2)))>eps;
}
int opposite_side(point3 p1,point3 p2,point3 l1,point3 l2){//判两点在线段异侧,点在线段上返回0,不共面无意义
    return dmult(xmult(subt(l1,l2),subt(p1,l2)),xmult(subt(l1,l2),subt(p2,l2)))<-eps;
}
int same_side(point3 p1,point3 p2,point3 s1,point3 s2,point3 s3){//判两点在平面同侧,点在平面上返回0
    return dmult(pvec(s1,s2,s3),subt(p1,s1))*dmult(pvec(s1,s2,s3),subt(p2,s1))>eps;
}
int opposite_side(point3 p1,point3 p2,point3 s1,point3 s2,point3 s3){//判两点在平面异侧,点在平面上返回0
    return dmult(pvec(s1,s2,s3),subt(p1,s1))*dmult(pvec(s1,s2,s3),subt(p2,s1))<-eps;
}
int parallel(point3 u1,point3 u2,point3 v1,point3 v2){//判两直线平行
    return vlen(xmult(subt(u1,u2),subt(v1,v2)))<eps;
}
int parallel(point3 u1,point3 u2,point3 u3,point3 v1,point3 v2,point3 v3){//判两平面平行
    return vlen(xmult(pvec(u1,u2,u3),pvec(v1,v2,v3)))<eps;
}
int parallel(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){//判直线与平面平行
    return zero(dmult(subt(l1,l2),pvec(s1,s2,s3)));
}
int perpendicular(point3 u1,point3 u2,point3 v1,point3 v2){//判两直线垂直
    return zero(dmult(subt(u1,u2),subt(v1,v2)));
}
int perpendicular(point3 u1,point3 u2,point3 u3,point3 v1,point3 v2,point3 v3){//判两平面垂直
    return zero(dmult(pvec(u1,u2,u3),pvec(v1,v2,v3)));
}
int perpendicular(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){//判直线与平面平行
    return vlen(xmult(subt(l1,l2),pvec(s1,s2,s3)))<eps;
}
int intersect_in(point3 u1,point3 u2,point3 v1,point3 v2){//判两线段相交,包括端点和部分重合
    if (!dots_onplane(u1,u2,v1,v2))
        return 0;
    if (!dots_inline(u1,u2,v1)||!dots_inline(u1,u2,v2))
        return !same_side(u1,u2,v1,v2)&&!same_side(v1,v2,u1,u2);
    return dot_online_in(u1,v1,v2)||dot_online_in(u2,v1,v2)||dot_online_in(v1,u1,u2)||dot_online_in(v2,u1,u2);
}
int intersect_ex(point3 u1,point3 u2,point3 v1,point3 v2){//判两线段相交,不包括端点和部分重合
    return dots_onplane(u1,u2,v1,v2)&&opposite_side(u1,u2,v1,v2)&&opposite_side(v1,v2,u1,u2);
}
int intersect_in(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){//判线段与空间三角形相交,包括交于边界和(部分)包含
    return !same_side(l1,l2,s1,s2,s3)&&!same_side(s1,s2,l1,l2,s3)&&
        !same_side(s2,s3,l1,l2,s1)&&!same_side(s3,s1,l1,l2,s2);
}
int intersect_ex(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){//判线段与空间三角形相交,不包括交于边界和(部分)包含
    return opposite_side(l1,l2,s1,s2,s3)&&opposite_side(s1,s2,l1,l2,s3)&&
        opposite_side(s2,s3,l1,l2,s1)&&opposite_side(s3,s1,l1,l2,s2);
}
point3 intersection(point3 u1,point3 u2,point3 v1,point3 v2){//计算两直线交点,注意事先判断直线是否共面和平行!
    point3 ret=u1;
    double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
        /((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
    ret.x+=(u2.x-u1.x)*t;
    ret.y+=(u2.y-u1.y)*t;
    ret.z+=(u2.z-u1.z)*t;
    return ret;
}
point3 intersection(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){//计算直线与平面交点,注意事先判断是否平行,并保证三点不共线!
    point3 ret=pvec(s1,s2,s3);
    double t=(ret.x*(s1.x-l1.x)+ret.y*(s1.y-l1.y)+ret.z*(s1.z-l1.z))/
        (ret.x*(l2.x-l1.x)+ret.y*(l2.y-l1.y)+ret.z*(l2.z-l1.z));
    ret.x=l1.x+(l2.x-l1.x)*t;
    ret.y=l1.y+(l2.y-l1.y)*t;
    ret.z=l1.z+(l2.z-l1.z)*t;
    return ret;
}
line3 intersection(point3 u1,point3 u2,point3 u3,point3 v1,point3 v2,point3 v3){//计算两平面交线,注意事先判断是否平行,并保证三点不共线!
    line3 ret;
    ret.a=parallel(v1,v2,u1,u2,u3)?intersection(v2,v3,u1,u2,u3):intersection(v1,v2,u1,u2,u3);
    ret.b=parallel(v3,v1,u1,u2,u3)?intersection(v2,v3,u1,u2,u3):intersection(v3,v1,u1,u2,u3);
    return ret;
}
double ptoline(point3 p,point3 l1,point3 l2){//点到直线距离
    return vlen(xmult(subt(p,l1),subt(l2,l1)))/distance(l1,l2);
}
double ptoplane(point3 p,point3 s1,point3 s2,point3 s3){//点到平面距离
    return fabs(dmult(pvec(s1,s2,s3),subt(p,s1)))/vlen(pvec(s1,s2,s3));
}
double linetoline(point3 u1,point3 u2,point3 v1,point3 v2){//直线到直线距离
    point3 n=xmult(subt(u1,u2),subt(v1,v2));
    return fabs(dmult(subt(u1,v1),n))/vlen(n);
}
double angle_cos(point3 u1,point3 u2,point3 v1,point3 v2){//两直线夹角cos值
    return dmult(subt(u1,u2),subt(v1,v2))/vlen(subt(u1,u2))/vlen(subt(v1,v2));
}
double angle_cos(point3 u1,point3 u2,point3 u3,point3 v1,point3 v2,point3 v3){//两平面夹角cos值
    return dmult(pvec(u1,u2,u3),pvec(v1,v2,v3))/vlen(pvec(u1,u2,u3))/vlen(pvec(v1,v2,v3));
}
double angle_sin(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){//直线平面夹角sin值
    return dmult(subt(l1,l2),pvec(s1,s2,s3))/vlen(subt(l1,l2))/vlen(pvec(s1,s2,s3));
}
int grid_onedge(int n,point* p){//多边形上的网格点个数
    int i,ret=0;
    for (i=0;i<n;i++)
        ret+=gcd(abs(p[i].x-p[(i+1)%n].x),abs(p[i].y-p[(i+1)%n].y));
    return ret;
}
int grid_inside(int n,point* p){//多边形内的网格点个数
    int i,ret=0;
    for (i=0;i<n;i++)
        ret+=p[(i+1)%n].y*(p[i].x-p[(i+2)%n].x);
    return (abs(ret)-grid_onedge(n,p))/2+1;
}
struct face { int a,b,c; bool ok; };
struct CH3D {
    static const int MAXN = 55;
    int n;//初始顶点数
    PT P[MAXN];//初始顶点
    int num;//凸包表面的三角形个数
    face F[8*MAXN];//凸包表面的三角形
    int g[MAXN][MAXN];
    double vol(PT &p, face &f) {
        PT m=P[f.b]-P[f.a];
        PT n=P[f.c]-P[f.a];
        PT t=p-P[f.a];
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
        PT a=P[F[s].a];
        PT b=P[F[s].b];
        PT c=P[F[s].c];
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
void rotate_to_horizontal(int n, PT *p, PT v) {//旋转点集，使法向量为v的平面水平
    if(sgn(v.x)==0 && sgn(v.y)==0) { return ; }
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
bool input(){
    scanf("%d", &n);
    if(n == 0) return false;
    hull.n = n;
    for(int i = 0; i < n; i++) { hull.P[i].in(); }
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
/*******directed_mst***/
struct Edge {
        int u, v;
        int cost;
}edge[M];
int pre[N], id[M], f[N], e;
item in[N];
void add_edge(int a, int b, int cost) { edge[e].u = a; edge[e].v = b; edge[e].cost = cost; e++; }
item directed_mst(int root, int n, int m) {
    item ret = 0;
    int u, v;
    while(1) {
        std::fill(in, in + n, inf);
        for(int i = 0; i < m; i++) {
            u = edge[i].u, v = edge[i].v;
            if(edge[i].cost < in[v] && u != v) {
                pre[v] = u;
                in[v] = edge[i].cost;
            }
        }
        for(int i = 0; i < n; i++) {
            if(i != root && in[i] == inf) {
                return -1;
            }
        }
        int cirs = 0;
        std::fill(id, id + n, -1);
        std::fill(f, f + n, -1);
        in[root] = 0;
        for(int i = 0; i < n; i++) {
            ret += in[i];
            int v = i;
            while(f[v] != i && id[v] == -1 && v != root) {
                f[v] = i;
                v = pre[v];
            }
            if(v != root && id[v] == -1) {
                for(int u = pre[v]; u != v; u = pre[u]) {
                    id[u] = cirs;
                }
                id[v] = cirs++;
            }
        }
        if(cirs == 0) { break; }
        for(int i = 0; i < n; i++) {
            if(id[i] == -1) {
                id[i] = cirs++;
            }
        }
        for(int i = 0; i < m; i++) {
            v = edge[i].v;
            edge[i].u = id[edge[i].u];
            edge[i].v = id[edge[i].v];
            if(edge[i].u != edge[i].v) {
                edge[i].cost -= in[v];
            }
        }
        n = cirs;
        root = id[root];
    }
    return ret;
}
/**********一般图匹配***************/
#define SET(a,b) memset(a,b,sizeof(a))
deque<int> Q;
//g[i][j]存放关系图：i,j是否有边,match[i]存放i所匹配的点
bool g[MAXN][MAXN],inque[MAXN],inblossom[MAXN];
int match[MAXN],pre[MAXN],base[MAXN];
int findancestor(int u,int v) {//找公共祖先
        bool inpath[MAXN]={false};
        while(1) {
                u=base[u];
                inpath[u]=true;
                if(match[u]==-1)break;
                u=pre[match[u]];
        }
        while(1) {
                v=base[v];
                if(inpath[v])return v;
                v=pre[match[v]];
        }
}
void reset(int u,int anc) {//压缩花
        while(u!=anc) {
                int v=match[u];
                inblossom[base[u]]=1;
                inblossom[base[v]]=1;
                v=pre[v];
                if(base[v]!=anc)pre[v]=match[u];
                u=v;
        }
}
void contract(int u,int v,int n) {
        int anc=findancestor(u,v);
        SET(inblossom,0);
        reset(u,anc);reset(v,anc);
        if(base[u]!=anc)pre[u]=v;
        if(base[v]!=anc)pre[v]=u;
        for(int i=1;i<=n;i++)
                if(inblossom[base[i]]) {
                        base[i]=anc;
                        if(!inque[i]) {
                                Q.push_back(i);
                                inque[i]=1;
                        }
                }
}
bool dfs(int S,int n) {
    for(int i=0;i<=n;i++)pre[i]=-1,inque[i]=0,base[i]=i;
    Q.clear();Q.push_back(S);inque[S]=1;
    while(!Q.empty()) {
        int u=Q.front();Q.pop_front();
        for(int v=1;v<=n;v++) {
            if(g[u][v]&&base[v]!=base[u]&&match[u]!=v) {
                if(v==S||(match[v]!=-1&&pre[match[v]]!=-1))contract(u,v,n);
                else if(pre[v]==-1) {
                    pre[v]=u;
                    if(match[v]!=-1)Q.push_back(match[v]),inque[match[v]]=1;
                    else {
                        u=v;
                        while(u!=-1) {
                            v=pre[u];
                            int w=match[v];
                            match[u]=v;
                            match[v]=u;
                            u=w;
                        }
                        return true;
                    } } } } } return false;
}
int solve(int n) {
    SET(match,-1);
    int ans=0;
    for(int i=1;i<=n;i++)
        if(match[i]==-1&&dfs(i,n))
            ans++;
    return ans;
}
/*********max_flow********************/
template<class T>
struct Max_Flow {
    int s, t, n;
    int Q[N], sign;
    int head[N], level[N], cur[N], pre[N];
    int nxt[M], pnt[M], E;
    T cap[M];
    void Init(int n, int s, int t) {
        this->n = n;
        this->s = s;
        this->t = t;
        E = 0;
        std::fill(head, head + n, -1);
    }
    void Add_edge(int from, int to, T c) {
        pnt[E] = to;
        cap[E] = c;
        nxt[E] = head[from];
        head[from] = E++;
        pnt[E] = from;
        cap[E] = 0;
        nxt[E] = head[to];
        head[to] = E++;
    }
    bool Bfs(int s, int t) {
        sign = t;
        std::fill(level, level + n, -1);
        int *front = Q, *tail = Q;
        *tail++ = t; level[t] = 0;
        while(front < tail && level[s] == -1) {
            int u = *front++;
            for(int e = head[u]; e != -1; e = nxt[e]) {
                if(cap[e ^ 1] > 0 && level[pnt[e]] < 0) {
                    level[pnt[e]] = level[u] + 1;
                    *tail ++ = pnt[e];
                }
            }
        }
        return level[s] != -1;
    }
    void Push(int t, T &flow) {
        T mi = INF;
        int p = pre[t];
        for(int p = pre[t]; p != -1; p = pre[pnt[p ^ 1]]) {
            mi = std::min(mi, cap[p]);
        }
        for(int p = pre[t]; p != -1; p = pre[pnt[p ^ 1]]) {
            cap[p] -= mi;
            if(!cap[p]) { sign = pnt[p ^ 1]; }
            cap[p ^ 1] += mi;
        }
        flow += mi;
    }
    void Dfs(int u, int t, T &flow) {
        if(u == t) { Push(t, flow); return ; }
        for(int &e = cur[u]; e != -1; e = nxt[e]) {
            if(cap[e] > 0 && level[u] - 1 == level[pnt[e]]) {
                pre[pnt[e]] = e;
                Dfs(pnt[e], t, flow);
                if(level[sign] > level[u]) {
                    return ;
                }
                sign = t;
            }
        }
    }
    T Dinic() {
        pre[s] = -1;
        T flow = 0;
        while(Bfs(s, t)) {
            std::copy(head, head + n, cur);
            Dfs(s, t, flow);
        }
        return flow;
    }
};
/***********费用流*******************/
template<int N,typename T> 
struct CostFlow {
    int s,t,head[N],etot,prevv[N],preve[N],inq[N],que[N],qf,qe;
    T dis[N];
    struct Edge {int v,next; T cap,cost;} g[501000];
    void init() {
        memset(head,-1,sizeof(head)); etot = 0;
    }
    void add_edge(int u,int v,T cap,T cost) {
        g[etot].v = v; g[etot].cap = cap; g[etot].cost = cost; g[etot].next = head[u]; head[u] = etot ++;
        g[etot].v = u; g[etot].cap = 0; g[etot].cost = -cost; g[etot].next = head[v]; head[v] = etot ++;
    }
    void mcmf(int _s,int _t,T &cost,T &flow) {
        s = _s; t = _t; cost = flow = 0;
        while (true) {
            for (int i = 0; i < N; i ++) dis[i] = (T)1e30;
            dis[s] = 0;
            qf = qe = 0;
            que[qe++] = s;
            while (qf!=qe) {
                int u = que[qf++]; inq[u] = 0; if (qf==N) qf = 0;
                for (int i = head[u]; i != -1; i = g[i].next) {
                    Edge &e = g[i];
                    if (e.cap && dis[e.v]>dis[u]+e.cost) {
                        dis[e.v] = dis[u]+e.cost;
                        prevv[e.v] = u; preve[e.v] = i;
                        if (!inq[e.v]) {
                            que[qe++] = e.v;
                            if (qe==N) qe = 0;
                            inq[e.v] = 1;
                        }
                    }
                }
            }
            if (dis[t]==T(1e30)) break;
            T f = (T)1e30;
            for (int u = t; u != s; u = prevv[u])
                f = min(f,g[preve[u]].cap);
            cost += f*dis[t];
            flow += f;
            for (int u = t; u != s; u = prevv[u])
                g[preve[u]].cap -= f,g[preve[u]^1].cap += f;
        }
    }
};
/************2sat***************/
template<int N> 
struct Sat {
    int head[N],etot,stack[N],top;
    bool mark[N];
    struct Edge {int v,next;} g[501000];
    void init() {
        memset(head,-1,sizeof(head)); etot = 0;
    }
    void addEdge(int u,int v) {
        g[etot].v = v; g[etot].next = head[u]; head[u] = etot ++;
    }
    bool dfs(int u) {
        if (mark[u^1]) return false;
        if (mark[u]) return true;
        mark[u] = true;
        stack[top++] = u;
        for (int i = head[u]; i != -1; i = g[i].next) {
            int v = g[i].v;
            if (!dfs(v)) return false;
        }
        return true;
    }
    bool work(int n) {
        memset(mark,0,sizeof(mark));
        for (int i = 0; i < n; i += 2) {
            if (!mark[i] && !mark[i+1]) {
                top = 0;
                if (!dfs(i)) {
                    while (top) mark[stack[--top]] = false;
                    if (!dfs(i+1)) return false;
                }
            }
        }
        return true;
    }
};
/*************KM*********************/
const int INF = 0x3f3f3f3f;
template<int N> 
struct KM {
    int lv[N],rv[N],left[N][N],n,m; // 左部n个点，右部m个点，有左部的完美匹配
    int la[N],ra[N],G[N][N];
    bool expath(int u) {
        lv[u] = 1;
        for (int i = 0; i < m; i ++) 
            if (!rv[i] && la[u]+ra[i]==G[u][i]) {
                rv[i] = 1;
                if (left[i]==-1 || expath(left[i]))
                    return left[i] = u,true;
            }
        return false;
    }
    void init() {
        for (int i = 0; i < N; i ++)
            for (int j = 0; j < N; j ++) G[i][j] = -INF;
    }
    int km(int _n,int _m) {
        n = _n; m = _m;
        memset(left,-1,sizeof(left));
        for (int i = 0; i < n; i ++) {
            la[i] = -INF;
            for (int j = 0; j < m; j ++) la[i] = std::max(la[i],G[i][j]);
        }
        for (int i = 0; i < m; i ++) ra[i] = 0;
        for (int u = 0; u < n; u ++) {
            for (int i = 0; i < n; i ++) lv[i] = 0;
            for (int i = 0; i < m; i ++) rv[i] = 0;
            while (!expath(u)) {
                int d = INF;
                for (int i = 0; i < n; i ++) if (lv[i])
                    for (int j = 0; j < m; j ++) if (!rv[j])
                        d = std::min(d,la[i]+ra[j]-G[i][j]);
                for (int i = 0; i < n; i ++) 
                    if (lv[i]) la[i] -= d,lv[i] = 0;
                for (int i = 0; i < m; i ++)
                    if (rv[i]) ra[i] += d, rv[i] = 0;
            }
        }
        int ret = 0;
        for (int i = 0; i < m; i ++) if (left[i]!=-1) ret += G[left[i]][i];
        return ret;
    }
};
/*****上下界网络流*************/
1：无源汇的可行流 ： 新建源点，汇点，M[i]为每个点进来的下界流减去出去的下界流，
如果M[i]为正，由源点向改点建M[i]的边，反之，由该点向汇点建M[i]的边，
原图中的边为每条边的上界建去下界。跑一遍最大流，每条边的流量加上下界流就是答案。
2：有源汇的最大流： 从汇点向源点建一条容量为INF的边，用上面的方法判断是否有解，
有解就再跑一遍从原图中源点到汇点的最大流
add_edge(T,S,INF,0);  
SAP(ss,tt,tt+1);  
for(i = head[ss];i!=-1;i=edge[i].next)  {  
    if(edge[i].c) return false;  
}  
return SAP(S,T,tt+1) == max_flow;  
3：有源汇的最小流：先跑一遍最大流，然后连上从汇点到源点的边，再跑一次最大流
SAP(ss,tt,tt+1);
add_edge(t,s,INF,0);
SAP(ss,tt,tt+1);
bool flag = false;
for(int i = head[ss];i!=-1;i=edge[i].next) {
    if(edge[i].c) {
        flag = true;
        break;
    }
}
if(flag) {
    puts("impossible");
}
else {
    int ans = 0;
    for(int i = head[t]; i != -1; i = edge[i].next) {
        if(edge[i].v == s ) {
            ans = edge[i^1].c;
            break;
        }
    }
    printf("%d\n",ans);
}
/*******matrix_tree*************/
#define zero(x) ((x>0? x:-x)<1e-15)
double a[N][N];
double b[N][N];
int g[110][110];
double det(double a[][N],int n) {
    int i, j, k, sign = 0;
    double ret = 1, t;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            b[i][j] = a[i][j];
    for (i = 0; i < n; i++) {
        if (zero(b[i][i])) {
            for (j = i + 1; j < n; j++)
                if (!zero(b[j][i]))
                    break;
            if (j == n)
                return 0;
            for (k = i; k < n; k++)
                swap(b[i][k], b[j][k]);
            sign++;
        }
        ret *= b[i][i];
        for (k = i + 1; k < n; k++)
            b[i][k] /= b[i][i];
        for (j = i + 1; j < n; j++)
            for (k = i + 1; k < n; k++)
                b[j][k] -= b[j][i] * b[i][k];
    }
    if (sign & 1)
        ret = -ret;
    return ret;
}
int main() {
    /* 0 based
     * a[i][i] = degree[i];
     * if(g[i][j]) {
     *     a[i][j]--;
     * }
     * ans = det(a, n - 1);
     */
}

/** 无向图最小环输出解，pre[i][j]表示i到j最短路径上的第一个点***/
void solve() {
    ret = INF;
    for(int k = 0; k < n; k++) {
        for(int i = 0; i < k; i++) {
            for(int j = i + 1; j < k; j++) {
                if(update(ret, map[i][k] + dis[i][j] + map[k][j])) {
                    tot = 0;
                    ans[tot++] = k;
                    int p = i;
                    while(p != -1) {
                        ans[tot++] = p;
                        p = pre[p][j];
                    } } } }
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if(dis[i][k] + dis[k][j] < dis[i][j]) {
                    pre[i][j] = pre[i][k];
                    dis[i][j] = dis[i][k] + dis[k][j];
                } } }
    }
}
/*********表达式计算*************/
struct Exp {
    bool error;
    int tot,top;
    int num[N];
    char op[N];
    void ini() {
        error=false;
        tot=0;
        top=1;
        op[1]='(';
    }
    bool prior(char a,char b) {
        if(b=='+'||b=='-')
            return a!='(';
        return a=='*'||a=='/';
    }
    int cal(char c,int a,int b) {
        if(c=='+') return a+b;
        if(c=='-') return a-b;
        if(c=='*') return a*b;
        if(b!=0) return a/b;
        error=true;
        return 0;
    }
    bool digit(char ch) {
        return ch>='0'&&ch<='9';
    }
    int solve(char *s,int len) {
        s[len++]=')';
        for(int i=0;i<len;i++) {
            if(s[i]=='(') op[++top]=s[i];
            else if(s[i]==')') {
                while(top>0&&op[top]!='(') {
                    num[tot-1]=cal(op[top],num[tot-1],num[tot]);
                    tot--;
                    top--;
                }
                top--;
            }
            else if(s[i]=='-'&&(i==0||s[i-1]=='(')) {
                num[++tot]=0;
                op[++top]='-';
            }
            else if(digit(s[i])) {
                int t=s[i]-'0';
                for(i++;digit(s[i]);i++) {
                    t=t*10+s[i]-'0';
                }
                num[++tot]=t;
                i--;
            }
            else {
                while(top>0&&prior(op[top],s[i])) {
                    num[tot-1]=cal(op[top],num[tot-1],num[tot]);
                    tot--; top--;
                }
                op[++top]=s[i];
            }
        }
        return num[1];
    }
}E;
/********点双连通缩点******/
struct Edge{
    int id;
    int s,t,next;
}list[M];
int head[N],E,dfn[N],low[N], stk[M], tot, n,nn, ok[N], belong[M],child, Btype,tdfn,ve[M];
void init(){
    memset(ve,false,sizeof(ve));
    memset(ok,false,sizeof(ok));
    tdfn = 0; Btype = 0;
    memset(dfn,0,sizeof(dfn));
}
void add(int u,int v,int _id){
}
void dfs(int u,int f){
    low[u] = dfn[u] = ++tdfn;
    for(int i = head[u]; i != -1; i = list[i].next){
        if(ve[list[i].id])  continue;
        int v = list[i].t;
        ve[list[i].id] = true;
        stk[++tot] = i;
        if(!dfn[v]) {    
            if(u == f) child ++;//f is root
            dfs(v,f);
            if(low[v] < low[u]) low[u] = low[v];
            if(low[v] >= dfn[u]){
                if(u!=f) ok[u] = true; // u is gedian 
                ++Btype;
                int id;
                do {
                    id = stk[tot--];
                    belong[id/2+1] = Btype;
                }while(list[id].s!=u);
            }
        }
        else if(dfn[v] < low[u]) low[u] = dfn[v];
    }
}
void SCC(){
    init();
    for(int i = 1; i <= n; i++) {
        child = 0; tot = 0;
        if(!dfn[i]) {
            dfs(i,i);
            ok[i] = false;
            if(child > 1) ok[i] = true;
        }
    }
    nn = 0;
    for(int i = 1; i <= n; i++) if(ok[i]) nn++;
    nn += Btype;
    for(int i = 0; i < N; i++) edge[i].clear();
    for(int i = 1,now=0; i <= n; i++)if(ok[i]){
        ++now;
        for(int j = head[i]; j != -1; j = list[j].next){
            edge[now+Btype].push_back(belong[j/2+1]);
            edge[belong[j/2+1]].push_back(now+Btype);
        }
    }
}
/***强连通缩点*****************/
struct Edge {
    int head[N], nxt[M], pnt[M], E;
    void init() {
        E = 0;
        memset(head, -1, sizeof(head));
    } void add_edge(int a, int b){
        pnt[E] = b;
        nxt[E] = head[a];
        head[a] = E++;
    }
}ori, ne, rne;
int in[N], bel[N], dfn[N], low[N], stk[N], n, m, Time, top, btype;
std::vector<int> con[N];
void dfs(int u) {
    low[u] = dfn[u] = ++Time; in[u] = 1; stk[++top] = u;
    for(int i = ori.head[u]; ~i; i = ori.nxt[i]) {
        int v = ori.pnt[i];
        if(!dfn[v]) {
            dfs(v);
            low[u]=std::min(low[u], low[v]);
        } else if(in[v]) {
            low[u]=std::min(low[u], low[v]);
        }
    }
    if(low[u] == dfn[u]) {
        btype++; int v;
        do {
            v = stk[top--], in[v] = 0, bel[v] = btype;
            con[btype].push_back(v);
        }while(v != u);
    }
}
void scc() {
    std::fill(in, in + n + 1, 0);
    std::fill(dfn, dfn + n + 1, 0);
    btype = Time = top = 0;
    for(int i = 1; i <= n; i++) if(!dfn[i]) {
        dfs(i);
    }
}
/************边双连通****************/
int Btype,Time,N,M,dfn[maxn],low[maxn],Belong[maxn],st[maxn],Top;
void dfs(int s) {
    int i,t,id;
    st[++Top]=s;
    dfn[s]=low[s]=++Time;
    for (i=head[s];i!=-1;i=edge[i].next) {	
        if(edge[i].vis)continue;	
        edge[i].vis=edge[i^1].vis=1;
        t=edge[i].t;
        if (!dfn[t]) {
            dfs(t);
            low[s]=min(low[s],low[t]);
        }
        else low[s]=min(low[s],dfn[t]);
    }
    if(dfn[s]==low[s]) {
        Btype++;
        do{
            t=st[Top--];
            Belong[t]=Btype;
        }while(t!=s);
    }
}
void SCC(int n) {
    int i;
    Time=0;Btype=0;Top=0;
    memset(dfn,0,sizeof(int)*(n+1));
    for(i=1;i<=n;i++)if(!dfn[i])
        dfs(i);
}
/*************tree_divide_conquer********/
int bfs(int root, int fa) {
    int head = 0, tail = 0;
    queue[tail++] = root;
    parent[root] = fa;
    while(head < tail) {
        int u = queue[head++];
        for(int i = first[u]; ~i; i = nxt[i]) {
            int v = pnt[i];
            if(v != parent[u] && !solved[v]) {
                queue[tail++] = v;
                parent[v] = u;
                depth[v] = depth[u] + weight[i];
            }
        }
    }
    return tail;
}
int get_root(int root, int &sz) {
    int tail = bfs(root, -1);
    for(int i = 0; i < tail; i++) {
        int u = queue[i];
        size[u] = mx[u] = 1;
    }
    for(int i = tail - 1; i >= 1; i--) {
        int u = queue[i];
        size[parent[u]] += size[u];
        mx[parent[u]] = std::max(mx[parent[u]], size[u]);
    }
    for(int i = 0; i < tail; i++) {
        int u = queue[i];
        mx[u] = std::max(mx[u], tail - size[u]);
        if(mx[u] < mx[root]) {
            root = u;
        }
    }
    sz = tail;
    return root;
}
/**********suffix array**********************/
int sa[N],X[N],Y[N],b[N],a[N],h[N],r[N];
bool comp(int *r,int a,int b,int le) {
    return r[a] == r[b] && r[a+le] == r[b+le];
}
void sort(int *Rank,int *Id,int n,int m) {
    std::fill(b,b+m,0);
    for(int i = n-1; i >= 0; i--) b[Rank[i]]++;
    for(int i = 1; i < m; i++) b[i] += b[i-1];
    for(int i = n-1; i >= 0; i--) sa[--b[Rank[Id[i]]]] = Id[i];
}
void calh(int n) {
    for(int i = 1; i <= n; i++) r[sa[i]] = i;
    int height = 0;
    for(int i = 0; i < n; i++){
        if(height) height--;
        int j = sa[r[i]-1];
        while(a[j+height]==a[i+height]) height++;
        h[r[i]] = height;
    }
}
void suffix(int n,int m=500) {
    int *Rank = X, *Id = Y, p = 1;
    for(int i = 0; i < n; i++) Rank[i] = a[i], Id[i] = i;
    sort(Rank,Id,n,m);
    for(int j = 1; p < n; j <<= 1){
        p = 0;
        for(int i = n-j; i < n; i++) Id[p++] = i;
        for(int i = 0; i < n; i++) if(sa[i] >= j) Id[p++] = sa[i] - j;
        sort(Rank,Id,n,p);
        std::swap(Rank,Id);
        Rank[sa[0]] = 0, p = 1;
        for(int i = 1; i < n; i++)
            Rank[sa[i]] = comp(Id,sa[i-1],sa[i],j) ? p-1 : p++;
        m = p;
    }
    calh(n-1);
}
struct Mancher {  
    char str[maxn];//start from index 1  
    int p[maxn];  
    char s[maxn];  
    int n;  
    void checkmax(int &ans,int b){  
        if(b>ans) ans=b;  
    }  
    inline int min(int a,int b){  
        return a<b?a:b;  
    }  
    void kp(){  
        int i;  
        int mx = 0;  
        int id;  
        for(i=1; i<n; i++){  
            if( mx > i )  
                p[i] = min( p[2*id-i], p[id]+id-i );  
            else  
                p[i] = 1;  
            for(; str[i+p[i]] == str[i-p[i]]; p[i]++) ;  
            if( p[i] + i > mx ) {  
                mx = p[i] + i;  
                id = i;  
            }  
        }  
    }  
    void pre()  {  
        int i,j,k;  
        n = strlen(s);  
        str[0] = '$';  
        str[1] = '#';  
        for(i=0;i<n;i++)  {  
            str[i*2 + 2] = s[i];  
            str[i*2 + 3] = '#';  
        }  
        n = n*2 + 2;  
        str[n] = 0;  
    }  
    void solve() {// 求出所有的最长回文子串所在的区间  
        int & tot = M::n;  
        tot = 0;  
        for(int i = 2; i < n; i++)  {  
            if(i%2&&p[i]==1) continue;  
            if(i%2)  {  M::in[tot++] = M::node(i/2-p[i]/2+1,i/2+p[i]/2);  }  
            else   {  M::in[tot++] = M::node(i/2-(p[i]/2-1),i/2+(p[i]/2-1));  }  
        }  
    }  
}task1;  
/************Aho-Corasick*************/
int ch[M][CD], val[M], fail[M], Q[M], used[M], sz, ID[256];
void init(){
    fail[0]=0;
    for(int i=0;i<CD;i++) ID[i+'a']=i;
}
void Reset(){
    memset(ch[0],0,sizeof(ch));
    sz=1;
}
void Insert(char *a){
    int p=0;
    for(;*a;a++)  {
        int c=ID[*a];
        if(!ch[p][c])  {
            memset(ch[sz],0,sizeof(ch[sz]));
            val[sz]=0;
            used[sz]=false;
            ch[p][c]=sz++;
        }
        p=ch[p][c];
    }
    val[p]++;
}
void Construct(){
    int *s=Q,*e=Q;
    for(int i=0;i<CD;i++) {
        if(ch[0][i]){
            fail[ ch[0][i] ] = 0;
            *e ++ = ch[0][i];
        }
    }while(s!=e){
        int u = *s++;
        for(int i = 0; i < CD ;i++){
            int &v = ch[u][i];
            if(v){
                *e++ = v;
                fail[v] = ch[ fail[u] ][i];
            } else  {
                v=ch[ fail[u] ][i];
            }
        }
    }
}
/**********KMP**********************/
inline void make(char *buf,int *fal) {
    static int i,j;
    fal[0]=-1;
    for(i=1,j=-1;buf[i];++i) {
        while(j>=0 && buf[j+1]!=buf[i])
            j=fal[j];
        if(buf[j+1]==buf[i])
            ++j;
        fal[i]=j;
    }
}
inline int match(char *p,char *t,int* fal) {
    static int i,j,re;
    for(re=i=0,j=-1;t[i];++i) {
        while(j>=0 && p[j+1]!=t[i])
            j=fal[j];
        if(p[j+1]==t[i])
            ++j;
        if(!p[j+1])
        {
            ++re;
            j=fal[j];
        }
    }
    return re;
}
//******Heavy-light-decomposition*********************/
int parent[N], success[N], depth[N], size[N];
void dfs(int u, int fa=-1) {
    parent[u] = fa, size[u] = 1, success[u] = -1;
    depth[u] = ~fa ? depth[fa] + 1 : 0;
    int v, mx = 0;
    for(int i = first[u]; ~i; i = nxt[i]) {
        if(fa == (v = pnt[i])) {
            continue;
        }
        dfs(v, u);
        if(size[v] > mx) {
            success[u] = v, mx = size[v];
        }
        size[u] += size[v];
    }
}
int head[N], id[N], cnt, mp[N];
void make_chain(int u, bool in_chain=false) {
    head[u] = in_chain ? head[parent[u]] : u;
    id[u] = ++cnt;
    mp[cnt] = u;
    if(~success[u]) {
        make_chain(success[u], true);
    }
    for(int i = first[u]; ~i; i = nxt[i]) {
        if(success[u] == pnt[i] || parent[u] == pnt[i]) {
            continue;
        }
        make_chain(pnt[i]);
    }
}
std::vector<int> operations[N];
void color(int a, int b, int c) {
    while(head[a] != head[b]) {
        if(depth[head[a]] < depth[head[b]]) {
            std::swap(a, b);
        }
        operations[id[head[a]]].push_back(c);
        operations[id[a] + 1].push_back(-c);
        a = parent[head[a]];
    }
    if(depth[a] > depth[b]) {
        std::swap(a, b);
    }
    operations[id[a]].push_back(c);
    operations[id[b] + 1].push_back(-c);
}
/************** treap*****************/
struct Node *nill;
struct Node {
    Node *ch[2];
    int val;
    int sz;

    void up() {
        if (this == nill) return ;
        sz = ch[0]->sz + ch[1]->sz + 1;
    }
};
void split(Node *a,Node *&b,Node *&c,int val) {
    if (a == nill) {
        b = c = nill;
    } else if (a->val <= val) {
        b = a;
        split(a->ch[1],b->ch[1],c,val);
        b->up();
    } else {
        c = a;
        split(a->ch[0],b,c->ch[0],val);
        c->up();
    }
}
unsigned ran() {
    static unsigned ranx = 233233233;
    return ranx += ranx << 2 | 1;
}
bool roll(int a,int b) {
    return ran() % (a+b) < a;
}
void merge(Node *&a,Node *b,Node *c) {
    if (b == nill) {
        a = c;
    } else if (c == nill) {
        a = b;
    } else if (roll(b->sz,c->sz)) {
        a = b;
        merge(a->ch[1],b->ch[1],c);
        a->up();
    } else {
        a = c;
        merge(a->ch[0],b,c->ch[0]);
        a->up();
    }
}
/*************splay************************/
struct Node *nill;
struct Node {
    int sz;
    bool rev;
    Node *ch[2],*fa;
    void up() {
        if (this==nill) return ;
        sz = ch[0]->sz + ch[1]->sz + 1;
    }
    void down() {
        if (this==nill) return ;
        if (!rev) return ;
        rev = false;
        ch[0]->reverse();
        ch[1]->reverse();
    }
    bool d() {
        return fa->ch[1]==this;
    }
    void reverse() {
        rev ^= 1;
        std::swap(ch[0],ch[1]);
    }
    void setc(Node *o,int c) {
        ch[c] = o;
        o->fa = this;
        up();
    }
    void rot() {
        int c = d(),cc = fa->d();
        Node *z = fa->fa;
        Node *tmp = fa;
        fa->setc(ch[c^1],c);
        setc(tmp,c^1);
        z->setc(this,cc);
    }
    void D() {
        if (this==nill) return ;
        fa->D();
        down();
    }
    void splay(Node *aim = nill) {
        D();
        while (fa!=aim) {
            if (fa->fa!=aim) {
                d()==fa->d() ? fa->rot() : rot();
            }
            rot();
        }
    }
};
/***link_cut_tree***********************/
const int N = 100000 + 5;
struct Node *nill;
struct Node {
    Node *fa,*ch[2];
    bool rev;
    void down() {
        if (this==nill) return ;
        if (rev==false) return ;
        rev = false;
        ch[0]->reverse();
        ch[1]->reverse();
    }
    void up() { if (this == nill) return ; }
    void reverse() { std::swap(ch[0],ch[1]); rev ^= 1; }
    bool d() { return fa->ch[1]==this; }
    bool isroot() { return fa==nill || fa->ch[0]!=this && fa->ch[1]!=this; }
    void setc(Node *p,int c) { ch[c] = p; p->fa = this; up(); }
    void rot() {
        Node *f = fa,*ff = fa->fa;
        int c = d(),cc = fa->d();
        f->setc(ch[c^1],c);
        this->setc(f,c^1);
        if (ff->ch[cc]==f) ff->setc(this,cc);
        else this->fa = ff;
    }
    void D() { if (this==nill) return ; fa->D(); down(); }
    Node *splay() {
        D();
        while (!isroot()) {
            if (!fa->isroot())
                d()==fa->d() ? fa->rot() : rot();
            rot();
        }
        return this;
    }
    Node *access() {
        for (Node *p = this,*q = nill; p!=nill; q = p,p = p->fa)
            p->splay()->setc(q,1);
        return splay();
    }
    void link(Node *p) { splay()->fa = p; }
    void rootable() { access()->reverse(); }
};
/************Math*******************/
1:∫0dx=c 不定积分的定义2:∫x^udx=(x^(u+1))/(u+1)+c 3:∫1/xdx=ln|x|+c 4)∫a^xdx=(a^x)/lna+c 5:∫e^xdx=e^x+c 
6:∫sinxdx=-cosx+c 7:∫cosxdx=sinx+c 8:∫1/(cosx)^2dx=tanx+c 9:∫1/(sinx)^2dx=-cotx+c 
10:∫1/√（1-x^2) dx=arcsinx+c 11:∫1/(1+x^2)dx=arctanx+c 12:∫1/(a^2-x^2)dx=(1/2a)ln|(a+x)/(a-x)|+c 
13:∫secxdx=ln|secx+tanx|+c 14:∫1/(a^2+x^2)dx=1/a*arctan(x/a)+c 15:∫1/√(a^2-x^2) dx=(1/a)*arcsin(x/a)+c 
16: ∫sec^2 x dx=tanx+c; 17: ∫shx dx=chx+c; 18: ∫chx dx=shx+c; 19: ∫thx dx=ln(chx)+c;
1. sum( k ) = n(n+1)/2
2. sum( 2k-1 ) = n^2
3. sum( k^2 ) = n(n+1)(2n+1)/6
4. sum( (2k-1)^2 ) = n(4n^2-1)/3
5. sum( k^3 ) = (n(n+1)/2)^2
6. sum( (2k-1)^3 ) = n^2(2n^2-1)
7. sum( k^4 ) = n(n+1)(2n+1)(3n^2+3n-1)/30
8. sum( k^5 ) = n^2(n+1)^2(2n^2+2n-1)/12
9. sum( k(k+1) ) = n(n+1)(n+2)/3
10. sum( k(k+1)(k+2) ) = n(n+1)(n+2)(n+3)/4
12. sum( k(k+1)(k+2)(k+3) ) = n(n+1)(n+2)(n+3)(n+4)/5
void FFT(Complex *a, int n, int rev) {
    for (int i = 1,j = 0; i < n; ++ i) {
        for (int k = n>>1; k > (j^=k); k >>= 1);
        if (i<j) swap(a[i],a[j]);
    }
    for (int m = 2; m <= n; m <<= 1) {
        Complex wm(cos(2*pi*rev/m),sin(2*pi*rev/m));
        for (int i = 0; i < n; i += m) {
            Complex w(1.0,0.0);
            for (int j = i; j < i+m/2; ++ j) {
                Complex t = w*a[j+m/2];
                a[j+m/2] = a[j] - t;
                a[j] = a[j] + t;
                w = w * wm;
            }
        }
    }
    if (rev==-1) {
        for (int i = 0; i < n; ++ i) a[i].x = (a[i].x/n+eps);
    }
}
void FFT(int *a, int n, int rev) {
    if(n == 1)  return ;
    memcpy(b, a, n<<2);
    for(int i = 0;i < n; i++) {
        a[(i%3)*(n/3)+i/3] = b[i];
    }
    FFT(a, n/3, rev);
    FFT(a+n/3, n/3, rev);
    FFT(a+n/3*2, n/3, rev);
    for(int i = 0;i < n; i++)   m2[i] = m1[i*(Top/n)];
    w[0] = 0; w[1] = (rev*n/3 + n)%n; w[2] = (rev*n/3*2 + n)%n;
    for(int i = 0;i < n/3; i++) {
        for(int j = 0;j < 3; j++) {
            b[i+j*(n/3)] = (a[i]+1ll*a[i+n/3]*m2[w[j]]+1ll*a[i+n/3*2]*m2[w[j]]%Mod*m2[w[j]])%Mod;
            w[j] += rev;
            if(w[j] >= n)   w[j] -= n;
            if(w[j] < 0)    w[j] += n;
        }
    }
    memcpy(a, b, n<<2);
}
void FFT(int *a, int n, int rev) {
    if(n == 1)  return ;
    for(int j = 0;j < n; j++) b[j] = a[j];
    for(int j = 0;j < n; j++) {
        a[(j%2)*(n/2)+j/2] = b[j];
    }
    FFT(a, n/2, rev);
    FFT(a+n/2, n/2, rev);
    for(int i = 0;i < n; i++)   m2[i] = m1[i*(Top/n)];
    int wm = rev==1 ? 1 : -1;
    int w[2];
    w[0] = 0; w[1] = n/2;
    for(int i = 0;i < n/2; i++) {
        for(int j = 0;j < 2;j ++) {
            b[i+j*(n/2)] = (a[i]+1ll*a[i+n/2]*m2[w[j]])%MOD;
            w[j] += wm;
            if(w[j] >= n)   w[j] -= n;
            if(w[j] < 0)    w[j] += n;
        }
    }
    for(int i = 0;i < n; i++)   a[i] = b[i];
    if(rev == -1 && n == realn) {
        int inv = pow_mod(n, MOD-2);
        for(int i = 0;i < n; i++)   a[i] = 1ll*a[i]*inv%MOD;
    }
}
pii linear(int A[], int B[], int M[], int n) {
    int x = 0, m = 1;
    for(int i = 0;i < n; i++) {
        int a = A[i]*m, b = B[i] - A[i]*x, d = gcd(M[i], a);
        if(b % d != 0)  return mkp(0, -1);
        int t = b/d * inv(a/d, M[i]/d) % (M[i]/d);
        x = x  +m*t;
        m *= M[i]/d;
    }
    x = (x%m + m)%m;
    return mkp(x, m);
}
double simpson(double a, double b) {
    double c = a + (b - a)/2;
    return (F(a) + 4*F(c) + F(b))*(b-a)/6;
}
double asr(double a, double b, double eps, double A) {
    double c = a + (b - a)/2;
    double L = simpson(a, c), R = simpson(c, b);
    if(fabs(L + R - A) <= 15*eps)   return L+R+(L+R-A)/15.0;
    return asr(a, c, eps/2, L) + asr(c, b, eps/2, R);
}
double asr(double a, double b, double eps) {
    return asr(a, b, eps, simpson(a, b));
}
// 精确覆盖
const int maxc = (1<<10)+5; // 各列结点数
const int maxn = (1<<15)+5; // 总的结点数
const int maxr = (1<<12)+5; // 各行结点数
struct DLX {
    int n, sz;  // 列数，结点总数
    int S[maxc];    // 各列当前结点数
    int row[maxn], col[maxn];   // 各结点行列编号
    int L[maxn], R[maxn], U[maxn], D[maxn]; // 十字链表
    int ansd, ans[maxr];    // 解
    void init(int n) {      // n是列数
        this->n = n; // 虚拟结点
        for(int i = 0;i <= n; i++) {
            U[i] = D[i] = i; L[i] = i-1; R[i] = i+1;
        }
        R[n] = 0; L[0] = n;
        sz = n+1;
        memset(S, 0, sizeof(S));
    }
    void addRow(int r, vector<int>  columns) {
        int first = sz;
        for(int i = 0;i < columns.size(); i++) {
            int c = columns[i];
            L[sz] = sz-1; R[sz] = sz+1; D[sz] = c; U[sz] = U[c];
            D[U[c]] = sz; U[c] = sz;
            row[sz] = r; col[sz] = c;
            S[c]++; sz++;
        }
        R[sz-1] = first; L[first] = sz-1;
    }
    // 顺着链表A， 遍历除s外的其他元素
#define FOR(i, A, s) for(int i = A[s];i != s; i = A[i])
    void remove(int c) {
        L[R[c]] = L[c];
        R[L[c]] = R[c];
        FOR(i, D, c) FOR(j, R, i) {
            U[D[j]] = U[j]; D[U[j]] = D[j]; --S[col[j]];
        }
    }
    void restore(int c) {
        FOR(i, D, c) FOR(j, R, i) {
            ++S[col[j]]; U[D[j]] = j; D[U[j]] = j;
        }
        L[R[c]] = c;
        R[L[c]] = c;
    }
    // d为递归深度
    bool dfs(int d) {
        if(R[0] == 0) {
            ansd = d;                        // 找到解
            return true;                   // 记录解的长度
        }
        // 找s最小的列c
        int c = R[0];                       // 第一个未删除的列
        FOR(i, R, 0) if(S[i] < S[c])    c = i;
        remove(c);                          // 删除第c列
        FOR(i, D, c) {
            ans[d] = row[i];
            FOR(j, R, i) remove(col[j]);    // 用结点i所在行覆盖第c列
            if(dfs(d+1))    return true;
            FOR(j, L, i)    restore(col[j]); // 恢复结点i所在行能覆盖的所有其他列
        }
        restore(c);                          // 恢复第c列
        return false;
    }
    bool solve(vector<int>  &v) {
        v.clear();
        if(!dfs(0)) return false;
        for(int i = 0;i < ansd; i++)
            v.push_back(ans[i]);
        return true;
    }
}dlx;
// 重复覆盖
const int maxc = 77;
const int maxn = 5000;
const int maxr = 77;
struct  DLX {
    int n, sz;
    int S[maxc];
    int row[maxn], col[maxn];
    int L[maxn], R[maxn], U[maxn], D[maxn];
    void init(int n) {
        this->n = n;

        for(int i = 0;i <= n; i++) {
            U[i] = D[i] = i; L[i] = i-1; R[i] = i+1;
        }
        R[n] = 0; L[0] = n;
        sz = n+1;
        memset(S, 0, sizeof(S));
    }
    void addRow(int r, vector<int>  columns) {
        int first = sz;
        for(int i = 0;i < columns.size(); i++) {
            int c = columns[i];
            L[sz] = sz-1; R[sz] = sz+1; D[sz] = c; U[sz] = U[c];
            D[U[c]] = sz; U[c] = sz;
            row[sz] = r; col[sz] = c;
            S[c]++; sz++;
        }
        R[sz-1] = first; L[first] = sz-1;
    }
#define FOR(i, A, s) for(int i = A[s];i != s;i = A[i])
    void remove(int x) {
        FOR(i, U, x) {
            L[R[i]] = L[i];
            R[L[i]] = R[i];
        }
    }
    void restore(int x) {
        FOR(i, D, x) {
            L[R[i]] = i;
            R[L[i]] = i;
        }
    }
    bool vis[maxc];
    int get_h() {
        for(int i = 0;i <= n; i++)  vis[i] = false;
        int ret = 0;
        FOR(i, R, 0) if(!vis[i]) {
            ret++; vis[i] = true;
            FOR(j, D, i) FOR(k, R, j)
                vis[col[k]] = true;
        }
        return ret;
    }
    int ans;
    void dfs(int d) {
        if(d+get_h() >= ans)   return ;
        if(R[0] == 0) {
            if(d < ans) ans = d;
            return ;
        }
        int c = R[0];
        FOR(i, R, 0) if(S[i] < S[c])    c = i;
        FOR(i, D, c) {
            remove(i);
            FOR(j, R, i)    remove(j);
            dfs(d+1);
            FOR(j, L, i)    restore(j);
            restore(i);
        }
    }
    int solve() {
        ans = 1111;
        dfs(0);
        return ans;
    }
}dlx;
/* 用于求整数解得方程组. */
const int maxn = 105;
// 有equ个方程，var个变元。增广阵行数为equ, 分别为0到equ - 1，列数为var + 1，分别为0到var.
int equ, var; 
int a[maxn][maxn];
int x[maxn]; // 解集.
bool free_x[maxn]; // 判断是否是不确定的变元.
//(-2表示有浮点数解，但无整数解，-1表示无解，0表示唯一解，大于0表示无穷解，并返回自由变元的个数)
int Gauss(void)
{
    int row = 0, col = 0;
    for ( ; row < equ && col < var; row++, col++) {
        // 枚举当前处理的行.
        int maxr = row;
        for (int i = row + 1; i < equ; i++)
            if (abs(a[i][col]) > abs(a[maxr][col])) maxr = i;
        if (maxr != row) // 与第k行交换.
            for(int i = row;i < var + 1; i++)   swap(a[row][i], a[maxr][i]);
        if (a[row][col] == 0){
            // 说明该col列第k行以下全是0了，则处理当前行的下一列.
            row--; continue;
        }
        for (int i = row + 1; i < equ; i++) {
            // 枚举要删去的行.
            if (a[i][col] != 0) {
                int LCM = lcm(abs(a[i][col]), abs(a[row][col]));
                int ta = LCM / abs(a[i][col]), tb = LCM / abs(a[row][col]);
                if (a[i][col] * a[row][col] < 0) tb = -tb; // 异号的情况是两个数相加.
                for (int j = col; j < var + 1; j++) {
                    a[i][j] = a[i][j] * ta - a[row][j] * tb;
                }
            }
        }
    }
    Debug();
    // 1. 无解的情况: 化简的增广阵中存在(0, 0, ..., a)这样的行(a != 0).
    for (int i = row; i < equ; i++)
        if (a[i][col] != 0) return -1;
    // 2. 无穷解的情况: 在var * (var + 1)的增广阵中出现(0, 0, ..., 0)这样的行
    if (row < var) {
        // 首先，自由变元有var - k个，即不确定的变元至少有var - k个.
        for (int i = row - 1; i >= 0; i--) {
            int free_x_num = 0, free_index;
            for (int j = 0; j < var; j++) {
                if (a[i][j] != 0 && free_x[j]) free_x_num++, free_index = j;
            }
            if (free_x_num > 1) continue; // 无法求解出确定的变元.
            int temp = a[i][var];
            for (int j = 0; j < var; j++)
                if (a[i][j] != 0 && j != free_index) temp -= a[i][j] * x[j];
            x[free_index] = temp / a[i][free_index]; // 求出该变元.
            free_x[free_index] = 0; // 该变元是确定的.
        }
        return var - row; // 自由变元有var - k个.
    }
    // 3. 唯一解的情况: 在var * (var + 1)的增广阵中形成严格的上三角阵.
    for (int i = var - 1; i >= 0; i--) {
        int temp = a[i][var];
        for (int j = i + 1; j < var; j++)
            if (a[i][j] != 0) temp -= a[i][j] * x[j];
        if (temp % a[i][i] != 0) return -2; // 说明有浮点数解，但无整数解.
        x[i] = temp / a[i][i];
    }
    return 0;
}
// 高斯消元解异或方程组   poj 1753
int a[22][22], x[22], var, equ, n, pos[22];
char s[22][22];
int min_sum;
void dfs(int row, int col) {
    if(col == -1 && row == -1) {
        int sum = 0;
        for(int i = 0;i < var; i++)   sum += x[i];
        min_sum = min(min_sum, sum);
        return ;
    }
    if(pos[row] == col) {
        x[col] = a[row][var];
        for(int i = var-1;i > col; i--) x[col] ^= x[i]&a[row][i];
        dfs(row-1, col-1);
    }
    else {
        x[col] = 1; dfs(row, col-1);
        x[col] = 0; dfs(row, col-1);
    }
}

int gauss() {
    int row = 0, col = 0;
    for( ; row < equ && col < var; row++, col++) {
        int maxr = row;
        for(int i = row+1;i < equ; i++) if(a[i][col])
            maxr = i;
        if(a[maxr][col] == 0) {
            row--; continue;
        }
        if(row != maxr) {
            for(int i = col;i < var+1; i++)
                swap(a[maxr][i], a[row][i]);
        }
        for(int i = row+1;i < equ; i++) {
            if(a[i][col]==0)    continue;
            for(int j = col;j < var+1; j++)
                a[i][j] ^= a[row][j];
        }
    }
    // 无解
    for(int i = row;i < equ; i++) if(a[i][var])
        return -1;
    if(row < var) {
        //  多解情况，求最小解
        min_sum = 1<<30;
        for(int i = 0;i < row; i++) {
            for(int j = 0;j < var; j++) if(a[i][j]) {
                pos[i] = j; break;
            }
        }
        dfs(row-1, var-1);
        return min_sum;
    }
    // 唯一解
    int ans = 0;
    for(int i = var-1;i >= 0; i--) {
        x[i] = a[i][var];
        for(int j = i+1;j < var; j++)
            x[i] ^= x[j]&a[i][j];
        ans += x[i];
    }
    return ans;
}
// 实数高斯消元
void gauss() {
    int n = var, r;
    for(int i = 0;i < n ;i++) {
        r = i;
        for(int j = i+1;j < n; j++) if(fabs(a[j][i]) > fabs(a[r][i]))
            r = j;
        if(r != i) for(int j = 0;j <= n; j++)   swap(a[r][j], a[i][j]);
        for(int j = n;j >= i; j--) {
            for(int k = i+1;k < n; k++) {
                a[k][j] -= a[k][i]/a[i][i]*a[i][j];
            }
        }
    }
    for(int i = n-1;i >= 0; i--) {
        x[i] = a[i][n];
        for(int j = i+1;j < n; j++)
            x[i] -= x[j]*a[i][j];
        x[i] /= a[i][i];
    }
}
//取模高斯消元 p为质数
void gauss()  {  
    //debug();  
    int i,j,k;  
    int row = 0 , col = 0;  
    for( ; row < equ && col < var; row++ , col++)  {  
        int maxr = row;  
        for(i = row+1;i < equ; i++)  
            if(abs(a[i][col]) >  abs(a[maxr][col]))  
                maxr = i;  
        if(a[maxr][col] == 0)  {  
            row--;  
            continue;  
        }  
        if(maxr != row)  
            for(i = col;i < equ+1;i ++)  
                swap(a[row][i] , a[maxr][i]);  
        for(i = row+1;i < equ; i++)  {  
            if(a[i][col] == 0)  
                continue;  
            LL LCM = lcm(abs(a[i][col]) , abs(a[row][col])) ;  
            LL aa = LCM / abs(a[row][col]) , bb = LCM / abs(a[i][col]);  
            if(a[row][col] * a[i][col] < 0)  
                aa = -aa;  
            for(j = col; j < var +1; j++)  
                a[i][j] = (a[i][j]*bb - aa*a[row][j])%p;  
        }  
    }  
    for(i = var-1;i >= 0; i--)  {  
        x[i] = a[i][var];  
        for(j = i + 1;j < var; j++)  
            x[i] -= a[i][j]*x[j];  
        if(x[i] == 0)  
            continue;  
        LL xx , yy;  
        LL ans = exgcd(a[i][i] , p , xx , yy);  
        xx = xx/ans*x[i];  
        x[i] = (xx%p+p)%p;  
    }  
    for(i = 0;i < var; i ++)  
        printf("%d ", x[i]);  
    puts("");  
}  
