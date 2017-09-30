const int MAX_LOG = 18;
const int N = 1 << MAX_LOG;

int p[MAX_LOG][N], dep[N], Time;
vector <int> children[N];

int L[N], R[N];
int lowbit_id[N], Log[N];

void log_init() {
    Log[0] = -1;
    for (int i = 1; i < N; i++) {
        Log[i] = Log[i >> 1] + 1;
    }
    lowbit_id[0] = -1;
    for (int i = 1; i < N; i++) {
        lowbit_id[i] = Log[lowbit(i)];
    }
}

void dfs(int u, int f) {
    dep[u] = dep[f] + 1;
    L[u] = ++Time;
    p[0][u] = f;
    rep(i, SIZE(children[u])) {
        int v = children[u][i];
        if (v != f) {
            dfs (v, u);
        }
    }
    R[u] = Time;
}

void init_lca(int n) {
    dfs(0, 0);
    for (int i = 1; i < MAX_LOG; i++) {
        rep (j, n) {
            p[i][j] = p[i - 1][p[i - 1][j]];
        }
    }
}

//get u's dth ancestor
inline int jump(int u, int d) {
    if (d < 0) {
        return u;
    }
    while (d) {
        u = p[lowbit_id[d]] [u];
        d -= 1 << lowbit_id[d];
    }
    return u;
}

int get_lca(int a,int b) {
    if(dep[a] > dep[b]) {
        swap(a,b);
    }
    int d = dep[b] - dep[a];
    b = jump(b, d);
    if(a == b) return a;
    for (int i = MAX_LOG - 1; i >= 0; i--) if (p[i][a] != p[i][b]) {
        a = p[i][a];
        b = p[i][b];
    }
    return p[0][a];
}
