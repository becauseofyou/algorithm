#include <cstdio>
#include <cstring>
#include <vector>
using namespace std;
vector<int> edge[1111];
bool sy[1111];
int match[1111];
bool dfs(int u) {
    for(vector<int>::iterator it = edge[u].begin(); it != edge[u].end(); ++it) {
        if(!sy[*it]) {
            sy[*it] = true;
            if(match[*it] == -1 || dfs(match[*it])) {
                match[*it] = u;
                return true;
            }
        }
    }
    return false;
}
int main() {
    int n;
    while(scanf("%d", &n) == 1) {
        for(int i = 0, k; i < n; i++) {
            edge[i].clear();
            scanf("%d: (%d)", &i, &k);
            for(int j = 0, x; j < k; j++) {
                scanf("%d", &x);
                edge[i].push_back(x);
            }
        }
        int ret = 0;
        memset(match, -1, sizeof(match));
        for(int i = 0; i < n; i++) {
            memset(sy, false, sizeof(sy));
            if(dfs(i)) {
                ret++;
            }
        }
        //ret实际上是最大匹配数的两倍
        printf("%d\n", n - ret/2);
    }
    return 0;
}
