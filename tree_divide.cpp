#include <cstdio>
#include <vector>
#include <cstring>
#include <algorithm>
const int N = 100010;
const int M = 200010;
int pool[N * 40];
int pool_size;
struct Bit
{
        void allocate(int size) {
                n = size;
                c = pool + pool_size;
                std::fill(c, c + n, 0);
                pool_size += size;
        }
        void insert(int x, int v) {
                for(x++; x < n; x += x&-x) {
                        c[x] += v;
                }
        }
        int query(int x) {
                x = std::min(x + 1, n - 1);
                int ret = 0;
                for(; x > 0; x -= x &-x) {
                        ret += c[x];
                }
                return ret;
        }
        int n;
        int *c;
}bit[2 * N], *node[2 * N];
int bit_size;
struct Center
{
        int u, v, d;
        Bit* bit;
        Bit* sub_bit;
        Center(int u, int v, int d, Bit* bit, Bit* sub_bit): u(u), v(v), d(d), bit(bit), sub_bit(sub_bit) {
        }
};
std::vector<Center>  edge[N];
int first[N], pnt[M], nxt[M], E, mx[N], w[N], n, solved[N], queue[N], parent[N], depth[N], size[N];
void init(int n)
{
        E = 0;
        std::fill(first, first + n + 1, -1);
        std::fill(solved, solved + n + 1, 0);
        bit_size = pool_size = 0;
        for(int i = 1; i <= n; i++) {
                edge[i].clear();
        }
}
void add_edge(int a, int b)
{
        pnt[E] = b;
        nxt[E] = first[a];
        first[a] = E++;
}
int bfs(int root, int fa)
{
        int head = 0, tail = 0;
        queue[tail++] = root;
        parent[root] = fa;
        while(head < tail) {
                int u = queue[head++];
                for(int i = first[u]; ~i; i = nxt[i]) {
                        int v = pnt[i];
                        if(v != parent[u] && !solved[v]) {
                                queue[tail++] = v;
                                parent[v] = u;
                                depth[v] = depth[u] + 1;
                        }
                }
        }
        return tail;
}
int get_root(int root, int &sz)
{
        int tail = bfs(root, -1);
        for(int i = 0; i < tail; i++) {
                int u = queue[i];
                size[u] = mx[u] = 1;
        }
        for(int i = tail - 1; i >= 1; i--) {
                int u = queue[i];
                size[parent[u]] += size[u];
                mx[parent[u]] = std::max(mx[parent[u]], size[u]);
        }
        for(int i = 0; i < tail; i++) {
                int u = queue[i];
                mx[u] = std::max(mx[u], tail - size[u]);
                if(mx[u] < mx[root]) {
                        root = u;
                }
        }
        sz = tail;
        return root;
}
void solve(int root)
{
        int sz = -1;
        depth[root] = 0;
        root = get_root(root, sz);
        solved[root] = 1;
        Bit &bitu = bit[bit_size];
        edge[root].push_back(Center(root, -1, 0, &bitu, &bitu));
        bitu.allocate(sz + 1);
        bitu.insert(0, w[root]);
        bit_size++;
        for(int i = first[root]; ~i; i = nxt[i]) {
                int v = pnt[i];
                if(solved[v]) {
                        continue;
                }
                depth[v] = 1; 
                int tail = bfs(v, root);
                bit[bit_size].allocate(tail + 2);
                for(int i = 0; i < tail; i++) {
                        int u = queue[i];
                        edge[u].push_back(Center(root, v, depth[u], &bitu, &bit[bit_size]));
                        bitu.insert(depth[u], w[u]);
                        bit[bit_size].insert(depth[u], w[u]);
                }
                bit_size++;
                solve(v);
        }
}
int main()
{
        int q, a, b;
        while(scanf("%d%d", &n, &q) == 2) {
                init(n);
                for(int i = 1; i <= n; i++) {
                        scanf("%d", &w[i]);
                }
                for(int i = 1; i < n; i++) {
                        scanf("%d%d", &a, &b);
                        add_edge(a, b);
                        add_edge(b, a);
                }
                solve(1);
                while(q--) {
                        char op[5];
                        scanf("%s%d%d", op, &a, &b);
                        if(op[0] == '?') {
                                int ret = 0;
                                for(int i = 0; i < (int)edge[a].size(); i++) {
                                        Center root = edge[a][i];
                                        if(b >= root.d) {
                                                ret += root.bit->query(b - root.d);
                                                if(~root.v) {
                                                        ret -= root.sub_bit->query(b - root.d);
                                                }
                                        }
                                }
                                printf("%d\n", ret);
                        } else {
                                b -= w[a];
                                for(int i = 0; i < (int)edge[a].size(); i++) {
                                        Center root = edge[a][i];
                                        root.bit->insert(root.d, b);
                                        if(~root.v) {
                                                root.sub_bit->insert(root.d, b);
                                        }
                                }
                                w[a] += b;
                        }
                }
        }
        return 0;
}
