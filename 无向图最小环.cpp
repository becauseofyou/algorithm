#include <cstdio>
#include <cstring>
const int INF = 1000000;
int map[110][110];
int dis[110][110];
int pre[110][110];
int n;
bool update(int &x, int y) {
        if(y < x) {
                x = y;
                return true;
        }
        return false;
}
int tot;
int ans[110];
int ret;
// pre[i][j]表示i到j最短路径上的第一个点
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
                                        }
                                }
                        }
                }
                for(int i = 0; i < n; i++) {
                        for(int j = 0; j < n; j++) {
                                if(dis[i][k] + dis[k][j] < dis[i][j]) {
                                        pre[i][j] = pre[i][k];
                                        dis[i][j] = dis[i][k] + dis[k][j];
                                }
                        }
                }
        }
}
void print() {
        if(ret == INF) {
                printf("No solution.\n");
        } else {
         //       printf("ret=%d\n",ret);
                for(int i = 0; i < tot; i++) {
                        if(i) printf(" ");
                        printf("%d", ans[i] + 1);
                }
                puts("");
        }
}
int main() {
        int t, m, a, b, c;
        while(scanf("%d", &n) != EOF) {
                if(n == -1) break;
                scanf("%d", &m);
                memset(pre, -1, sizeof(pre));
                for(int i = 0; i < n; i++) {
                        map[i][i] = dis[i][i] = 0;
                        for(int j = 0; j < n; j++) if(i != j) {
                                dis[i][j] = map[i][j] = INF;
                        }
                }
                for(int i = 0; i < m; i++) {
                        scanf("%d%d%d", &a, &b, &c);
                        a--; b--;
                        if(c < map[a][b]) {
                                dis[b][a] = map[b][a] = c;
                                dis[a][b] = map[a][b] = c;
                                pre[a][b] = b;
                                pre[b][a] = a;
                        }
                }
                solve();
                print();
        }
        return 0;
}
