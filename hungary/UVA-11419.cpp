/*
 * https://vjudge.net/problem/UVA-11419
 * 最小点覆盖 输出答案
 * 左边为行号 右边为列号
 * 每一个格点相当于一条边，连接一个行 和 一个列 号
 * 最终答案就是在二分图中选择最少的点覆盖所有的边
 *
 * 二分图中一共有这么些边
 * 1：匹配边
 * 2：左盖右未盖
 * 3：左不盖右盖
 * 4：左盖右盖 但非匹配边
 *
 * 2，3两种情况是类似的
 *
 * 文字的局限性
 */
#include <cstdio>
#include <cstring>
#include <vector>
using namespace std;
const int N = 1111;

vector<int> edge[N];
bool sx[N], sy[N];
int mx[N], my[N];

bool dfs(int u) {
    sx[u] = true;
    for(int v: edge[u]) {
        if(!sy[v]) {
            sy[v] = true;
            if(my[v] == -1 || dfs(my[v])) {
                mx[u] = v;
                my[v] = u;
                return true;
            }
        }
    }
    return false;
}

int main () {
    int n, m, k, x, y;
    while(scanf("%d%d%d", &n, &m, &k) == 3) {
        if(n + m + k == 0) break;
        for(int i = 1; i <= n; i++) edge[i].clear();
        for(int i = 0; i < k; i++) {
            scanf("%d%d", &x, &y);
            edge[x].push_back(y);
        }
        memset(my, -1, sizeof(my));
        memset(mx, -1, sizeof(mx));
        int ret = 0;
        for(int i = 1; i <= n; i++) {
            memset(sy, false, sizeof(sy));
            if(dfs(i)) 
                ret++;
        }
        vector<int> row, col;
        memset(sy, false, sizeof(sy));
        memset(sx, false, sizeof(sx));
        for(int i = 1; i <= n; i++) if(mx[i] == -1) {
            dfs(i);
        }
        for(int i = 1; i <= n; i++) if(!sx[i]) {
            row.push_back(i);
        }
        for(int i = 1; i <= m; i++) if(sy[i]) {
            col.push_back(i);
        }
        printf("%d", ret);
        for(int r : row) printf(" r%d", r); 
        for(int c : col) printf(" c%d", c);
        puts("");
    }
    return 0;
}
