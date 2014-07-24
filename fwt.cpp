//hdu 4867
#include <cstdio>
#include <cstring>
#include <algorithm>
const int N = 20010;
const int V = 1024;
const int MOD = 1000000007;
int temp[V], temp1[V];

void fwt(int *a, int *b, int l, int r) {
        if(l == r) {
                b[l] = a[l];
                return ;
        }
        int mid = (l + r) >> 1, half = (r - l + 1) >> 1;
        fwt(a, b, l, mid), fwt(a, b, mid + 1, r);
        std::copy(b + l, b + r + 1, temp + l);
        for(int i = l; i <= mid; i++) {
                b[i] = temp[i] - temp[i + half];
                if(b[i] < 0) 
                        b[i] += MOD;
        }
        for(int i = mid + 1; i <= r; i++) {
                b[i] = temp[i] + temp[i - half];
                if(b[i] >= MOD)
                        b[i] -= MOD;
        }
}
int inv2;
void rfwt(int *a, int *b, int l, int r) {
        if(l == r) {
                b[l] = a[l];
                return ;
        }
        int mid = (l + r) >> 1, half = (r - l + 1) >> 1;
        int temp[V];
        for(int i = l; i <= mid; i++) {
                temp[i] = (1LL * (a[i] + a[i + half]) * inv2) % MOD;
        }
        for(int i = mid + 1; i <= r; i++) {
                temp[i] = (1LL * (a[i] - a[i - half]) * inv2) % MOD;
                if(temp[i] < 0) 
                        temp[i] += MOD;
        }
        rfwt(temp, b, l, mid);
        rfwt(temp, b, mid + 1, r);
}
int tot;
struct Seg {
        struct Node {
                int cnt[V];
        }node[V << 2];
        int ans[V];
        void init() {
                for(int i = 0; i < tot; i++) {
                        temp[i] = 0;
                }
                temp[0] =  1;
                fwt(temp, temp1, 0, tot - 1);
                for(int i = 0; i < (tot << 2); i++) {
                        for(int j = 0; j < tot; j++) {
                                node[i].cnt[j] = temp1[j];
                        }
                }
                std::fill(ans, ans + tot, 1);
                ans[0] = 1;
        }
        void insert(int p, int *a, int l, int r, int x) {
                if(l == r) {
                        std::copy(a, a + tot, node[x].cnt);
                        return ;
                }
                int m = (l + r) >> 1;
                if(p <= m) insert(p, a, l, m, x + x);
                else insert(p, a, m + 1, r, x + x + 1);
                for(int i = 0; i < tot; i++) {
                        node[x].cnt[i] = 1LL * node[x + x].cnt[i] * node[x + x + 1].cnt[i] % MOD;
                }
        }
        void update() {
                rfwt(node[1].cnt, ans, 0, tot - 1);
        }
}tree;
int a[N], cnt[V], tmp[V];
void gan(int w) {
        if(cnt[w] == 0) {
                for(int i = 0; i < tot; i++) {
                        temp[i] = 0;
                }
                temp[0] = 1;
                fwt(temp, temp1, 0, tot - 1);
                tree.insert(w, temp1, 0, tot - 1, 1);
                return ;
        }
        for(int i = 0; i < tot; i++) {
                temp[i] = (i <= w);
                temp1[i] = 1;
        }
        fwt(temp, tmp, 0, tot - 1);
        int c = cnt[w];
        while(c) {
                if(c & 1) 
                        for(int i = 0; i < tot; i++) temp1[i] = 1LL * temp1[i] * tmp[i] % MOD;
                c >>= 1;
                for(int i = 0; i < tot; i++) tmp[i] = 1LL * tmp[i] * tmp[i] % MOD;
        }
        tree.insert(w, temp1, 0, tot - 1, 1);
        tree.update();
}
int main() {
        inv2 = 1;
        int tm = 2, p = MOD - 2;
        while(p) {
                if(p & 1) {
                        inv2 = 1LL * inv2 * tm % MOD;
                }
                p >>= 1;
                tm = 1LL * tm * tm % MOD;
        }
        int t, n, m, x, y;
        scanf("%d", &t);
        while(t--) {
                std::fill(cnt, cnt + V, 0);
                scanf("%d%d", &n, &m);
                int mx = 0;
                for(int i = 0; i < n; i++) {
                        scanf("%d", &a[i]);
                        cnt[a[i]]++;
                        if(a[i] > mx) {
                                mx = a[i];
                        }
                }
                for(tot = 1; (1<<tot) < mx; tot++);
                tot = 1 << tot;
              //  printf("tot=%d\n", tot);
                tree.init();
                for(int i = 0; i < V; i++) {
                        gan(i);
                }
                char op[10];
                for(int i = 0; i < m; i++) {
                        scanf("%s", op);
                        if(op[0] == 'C') {
                                scanf("%d%d", &x, &y);
                                cnt[a[x]]--;
                                gan(a[x]);
                                a[x] = y;
                                cnt[y]++;
                                gan(y);
                        } else {
                                scanf("%d", &x);
                                printf("%d\n", tree.ans[x]);
                        }
                }
        }
        return 0;
}
