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
