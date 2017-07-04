//×îÐ¡µã¸²¸Ç
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
    int n, m, k, mode1, mode2, job;
    while(scanf("%d", &n) == 1) {
        if(n == 0) break;
        scanf("%d%d", &m, &k);
        for(int i = 0; i < n; i++) edge[i].clear();
        for(int i = 0; i < k; i++) {
            scanf("%d%d%d", &job, &mode1, &mode2);
            if(mode1 && mode2)
                edge[mode1].push_back(mode2);
        }
        memset(match, -1, sizeof(match));
        int ret = 0;
        for(int i = 0; i < n; i++) {
            memset(sy, false, sizeof(sy));
            if(dfs(i)) {
                ret++;
            }
        }
        printf("%d\n", ret);
    }
    return 0;
}
