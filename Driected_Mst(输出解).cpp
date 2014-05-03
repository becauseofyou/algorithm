#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>

#define Clr(a,b) memset(a,b,sizeof(a))

using namespace std;

const int MAXN = 100010;

struct COST;
vector<COST*> cov;

struct COST {
    int c, id, used;
    COST * a, * b; // c = a - b;

    bool operator < (const COST & t) const {
        return c < t.c;
    }
    COST(int c, int id):c(c),id(id),a(0),b(0),used(0) {
        cov.push_back(this);
    }
    COST(COST * a, COST * b):a(a),b(b),c(a->c-b->c),id(-1),used(0) {
        cov.push_back(this);
    }
    void push() {
        if(id == -1) {
            a -> used += used;
            b -> used -= used;
        }
    }
    void use() { ++used; }
};

struct EDGE {
    int u, v;
    COST * c;

    int from, to, val;

    EDGE(){}
    EDGE(int u, int v, int cost, int id):u(u),v(v) {
        c = new COST(cost, id);
        from = u;
        to = v;
        val = cost;
    }
}   edge[MAXN];

int edgel;

inline bool cmp(COST * a, COST * b) { // a < b
    if(a == 0 || b == 0) return b == 0;
    else return a->c < b->c;
}

bool better(COST *a, COST *b) { //a better than b?
    if (a == 0 || b == 0)
        return b == 0;
    return a->c < b->c;
}

bool Directed_MST(int st, int nV, int nE, int &ret) {
    int i, j, k, idx[MAXN], cnt, id, vis[MAXN];
    int in[MAXN];//pre
    COST * least[MAXN];

    ret = 0;
    while(true) {
        for(i=0; i<nV; ++i) least[i] = NULL, in[i] = -1;
        for(i=0; i<nE; ++i) {
            int &u = edge[i].u, &v = edge[i].v;
            if(u != v && cmp(edge[i].c, least[v])) {
                least[v] = edge[i].c;
                in[v] = u;
            }
        }
        for(i=0; i<nV; ++i) if(i != st && in[i] == -1) return false;
        Clr(vis, -1);
        Clr(idx, -1);
        for(i=cnt=id=0, in[st]=-1; i<nV; ++i) if(i!=st) {
            ret += least[i]->c;
            least[i] -> use();
            if(~vis[i]) continue;
            for(j=i, k=-1; ; j=in[j]) {
                if(~vis[j]) {
                    if(vis[j] == cnt) k = j;
                    break;
                }
                vis[j] = cnt;
            }
            if(~k) {
                for(j=in[k]; j!=k; j=in[j]) idx[j] = id;
                idx[k] = id++;
            }
            ++cnt;
        }
        if(!id) break;
        //重标号
        for(i=0; i<nV; ++i) if(idx[i] == -1) idx[i] = id++;
        for(i=0; i<nE; ++i) {
            int v = edge[i].v;
            edge[i].u = idx[edge[i].u];
            edge[i].v = idx[edge[i].v];
            if(edge[i].u != edge[i].v)
                edge[i].c = new COST(edge[i].c, least[v]);
        }
        nV = id;
        st = idx[st];
    }
    return true;
}

int main()
{
    //freopen("input.txt","r",stdin);
    //freopen("output.txt","w",stdout);

    int n, m;
    scanf("%d%d", &n, &m);
    edgel = 0;
    for(int i=0; i<m; i++)
    {
        int x, y, z;
        scanf("%d%d%d", &x, &y, &z);
        --x, --y;
        edge[edgel++] = EDGE(x, y, z, i+1);
    }
    int ans;
    if(!Directed_MST(0, n, m, ans))
        puts("-1");
    else
    {
        printf("%d\n", ans);
        for(int i=cov.size()-1; i>=0; --i)
            cov[i] -> push();
        vector<int> tv;
        for(int i=0; i<cov.size(); ++i)
        {
            COST * c = cov[i];
            if(c -> id != -1 && c -> c > 0 && c -> used > 0)
                tv.push_back(c -> id);
        }
        sort(tv.begin(), tv.end());
        for(int i=0; i<tv.size(); ++i) printf("%d ", tv[i]);
        printf("\n");
    }
    return 0;
}
