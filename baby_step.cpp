#include <cstdio>
#include <algorithm>
#include <cmath>
#include <map>
#include <cstring>
const int MOD = 100000007;
int k, ret;
int pow_mod(int a, int b)
{
        long long ret = 1;
        while(b) {
                if(b & 1) {
                        ret = ret * a % MOD;
                }
                b >>= 1; a = 1LL * a * a % MOD;
        }
        return (int)ret;
}
int exgcd(int a, int b, int &x, int &y) 
{
        int gcd = a;
        if(b == 0) {
                x = 1;
                y = 0;
        } else {
                gcd = exgcd(b, a % b, x, y);
                x -= (a / b) * y;
                std::swap(x, y);
        }
        return gcd;
}
int main()
{
        int T, ca = 1;
        scanf("%d", &T);
        while(T--) {
                scanf("%d%d", &k, &ret);
                int m = (int)sqrt(1.0 * MOD) + 10;
                std::map<int, int> mp;
                for(int i = 0; i < m; i++) {
                        int tmp = pow_mod(k, i);
                        if(mp.find(tmp) == mp.end()) {
                                mp[tmp] = i;
                        }
                }
                int ans = MOD;
                for(int i = 0; i < m; i++) {
                        int tmp = pow_mod(pow_mod(k, i), m);
                        int x, y;
                        exgcd(tmp, MOD, x, y);
                        while(x < 0) x += MOD;
                        x = 1LL * x * ret % MOD;
                        y = 1LL * y * ret % MOD;
                        if(mp.find(x) != mp.end()) {
                                tmp = i * m + mp[x];
                                if(tmp < ans) {
                                        ans = tmp;
                                }
                        }
                }
                printf("Case %d: %d\n", ca++, ans);
        }
        return 0;
}

