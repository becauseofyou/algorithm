//最小路径覆盖,每增加一个匹配相当于结合了两条路径
//可以少去一条路径的代价
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

int main () {
    int t, n, m, a, b;
    scanf("%d", &t);
    while(t--) {
        scanf("%d%d", &n, &m);
        for(int i = 0; i < n; i++) edge[i].clear();
        for(int i = 0; i < m; i++) {
            scanf("%d%d", &a, &b);
            a--; b--;
            edge[a].push_back(b);
        }
        memset(match, -1, sizeof(match));
        int ret = 0;
        for(int i = 0; i < n; i++) {
            memset(sy, false, sizeof(sy));
            if(dfs(i)) {
                ret++;
            }
        }
        printf("%d\n", n - ret);
    }
    return 0;
}
