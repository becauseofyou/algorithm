#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <iomanip>
using namespace std;

#define MP make_pair

const double eps = 1e-8;

int dcmp( double x ){
    return (x > eps) - (x < -eps) ;
}

struct point {
    double x, y;
    point(){}
    point(double x, double y):x(x), y(y){}
    point operator - (const point &b) const {
        return point(x - b.x, y - b.y);
    }
    void input() {
        scanf( "%lf%lf", &x, &y );
    }
};

double cross( point a, point b ) {
    return a.x * b.y - a.y * b.x;
}


double dot( point a, point b ) {
    return a.x * b.x + a.y * b.y;
}

struct polygon {
    point p[5];
    int sz;
    void init() {
        p[sz] = p[0];
    }
}g[505];

pair<double, int> c[100000];

double segP( point a, point b, point c ) {
    if( dcmp(b.x - c.x) )
        return (a.x - b.x) / (c.x - b.x);
    return (a.y - b.y) / (c.y - b.y);
}

double polyUnion( int n )
{
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

int main()
{
    int n;
    while( scanf("%d", &n) != EOF ) {
        double sum = 0;
        for( int i = 0; i < n; ++i ) {
            //scanf("%d", &g[i].sz);
            g[i].sz = 4;
            for( int j = 0; j < g[i].sz; ++j )
                g[i].p[j].input();
            g[i].init();
            double tmp = 0;
            for( int j = 1; j <= 4; ++j )
                tmp += cross(g[i].p[j-1], g[i].p[j]);
            sum += fabs(tmp);
            if( tmp < 0 ) swap(g[i].p[1], g[i].p[3]);
        }
        printf( "%.10lf\n", sum * 0.5 / polyUnion(n) );
    }
    return 0;
}
