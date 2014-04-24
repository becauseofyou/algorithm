#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cmath>
typedef long long type;
struct comp{
    double x, y;
    comp(double _x=0, double _y=0) : x(_x), y(_y) {}
};
namespace FFT{
    const int N = 1<<18;
    const double pi2 = 3.1415926535897932 * 2;
    comp a[N], b[N], tmp[N];
    int n, bn;
    type res[N];
    inline comp W(int n, bool inv) {
        double ang = inv ? -pi2 / n : pi2 / n;
        return comp(cos(ang), sin(ang));
    }
    int bitrev(int x) {
        int ans = 0;
        for (int i=1; i<=bn; ++i)
            ans <<= 1, ans |= x & 1, x >>= 1;
        return ans;
    }
    void dft(comp *a,bool inv) {
        int step, to; comp w, wi, A, B;
        for (int i=0; i<n; ++i) {
            to = bitrev(i);
            if (to > i) std::swap(a[to], a[i]);
        }
        for (int i=1; i<=bn; ++i) {
            wi = W(1<<i, inv); w = comp(1, 0);
            step = 1 << (i-1);
            for (int k=0; k<step; ++k) {
                for (int j=0; j<n; j+=1<<i) {
                    int t = j | k, d = j|k|step;
                    A = a[t];
                    B.x  = w.x * a[d].x - w.y * a[d].y;
                    B.y  = w.x * a[d].y + w.y * a[d].x;
                    a[t].x = A.x + B.x, a[t].y = A.y + B.y;
                    a[d].x = A.x - B.x, a[d].y = A.y - B.y;
                }
                comp tmp;
                tmp.x = w.x * wi.x - w.y * wi.y;
                tmp.y = w.x * wi.y + w.y * wi.x;
                w = tmp;
            }
        }
    }
    int mul(int n1, int *x1, int n2, int *x2) {
        n = std::max(n1, n2);
        for (bn = 0; (1<<bn) < n; ++bn); ++bn;
        n = 1 << bn;
        for (int i=0; i<n; ++i) a[i] = b[i] = comp(0, 0);
        for (int i=0; i<n1; ++i) a[i] = comp(x1[i], 0);
        for (int i=0; i<n2; ++i) b[i] = comp(x2[i], 0);
        dft(a, false); dft(b, false);
        for (int i=0; i<n; ++i) {
            tmp[i].x = a[i].x * b[i].x - a[i].y * b[i].y;
            tmp[i].y = a[i].x * b[i].y + a[i].y * b[i].x;
        }
        dft(tmp, true);
        for (int i=0; i<n; ++i) res[i] = (type)(tmp[i].x/n + 0.1);
        for (--n; n && !res[n]; --n);
        return n+1;
    }
}
