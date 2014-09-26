#include <cstdio>
#include <vector>
#include <cstring>
#include <algorithm>
using std::pair;
using std::make_pair;
const int N = 100010;
const int M = 200010;

struct Node* null;
 
struct Node 
{
        Node* c[2];
        Node* fa;
        int val, depth, mx;
        int id;
        int mi_id;
        inline void setc(int d, Node *s) {
                c[d] = s;
                s->fa = this;
        }
        inline bool d() {
                return fa->c[1] == this;
        }
        void up() {
                mx = val; mi_id = id;
                if(c[0]->mx > mx || c[0]->mx == mx && c[0]->mi_id < mi_id) {
                        mx = c[0]->mx;
                        mi_id = c[0]->mi_id;
                }
                if(c[1]->mx > mx || c[1]->mx == mx && c[1]->mi_id < mi_id) {
                        mx = c[1]->mx;
                        mi_id = c[1]->mi_id;
                }
        }
        void clear() {
                c[0] = c[1] = fa = null;
        }
}node[N * 40], *pool;
struct SplayTree
{
        Node* root;
        void rot(Node* x) 
        {
                Node* p = x->fa;
                int f = x->d();
                p->fa->setc(p->d(), x);
                p->setc(f, x->c[!f]);
                x->setc(!f, p);
                p->up();
        }
        void splay(Node* x, Node* goal) 
        {
                while(x->fa != goal) {
                        if(x->fa->fa == goal) {
                                rot(x);
                        } else {
                                x->d() == (x->fa->d()) ? rot(x->fa): rot(x);
                                rot(x);
                        }
                }
                x->up();
                if(goal == null) {
                        root = x;
                }
        }
        void insert(int depth, int value, int ind)
        {
                if(root == null) {
                        Node* tmp = pool++;
                        tmp->clear();
                        tmp->val = value;
                        tmp->id = ind;
                        tmp->mi_id = ind;
                        tmp->depth = depth;
                        tmp->mx = value;
                        root = tmp;
                        return ;
                }
                Node* x = root;
                Node* p = null;
                int c = -1;
                while(x != null) {
                        p = x;
                        if(x->depth == depth) {
                                break;
                        }
                        if(x->depth < depth) {
                                x = x->c[1];
                                c = 1;
                        } else {
                                x = x->c[0];
                                c = 0;
                        }
                }
                if(x != null) {
                        splay(x, null);
                        if(value > x->val || value == x->val && ind < x->id) {
                                x->id = ind;
                                x->val = value;
                                x->up();
                        }
                        return ;
                }
                Node* tmp = pool++;
                tmp->clear();
                tmp->val = value;
                tmp->id = ind;
                tmp->depth = depth;
                p->setc(c, tmp);
                splay(tmp, null);
        }
        pair<int, int> query(int d)
        {
                int ret = -1, mx = -1;
                Node* x = root;
                Node* p = null;
                while(x != null) {
                        p = x;
                        if(x->depth <= d) {
                                if(x->c[0]->mx > mx || x->c[0]->mx == mx && x->c[0]->mi_id < ret) {
                                        mx = x->c[0]->mx;
                                        ret = x->c[0]->mi_id;
                                }
                                if(x->val > mx || x->val == mx && x->id < ret) {
                                        mx = x->val;
                                        ret = x->id;
                                }
                                x = x->c[1];
                        } else {
                                x = x->c[0];
                        }
                }
                if(p != null)
                        splay(p, null);
                return make_pair(ret, mx);
        }
        void init()
        {
                root = null;
        }
}spt[N * 40];
int spt_size;
struct Center
{
        int d;
        SplayTree* root;
        Center(int d, SplayTree* root): d(d), root(root) {
        }
};
std::vector<Center>  edge[N];
int first[N], pnt[M], nxt[M], weight[M], E, mx[N], w[N], n, solved[N], queue[N], parent[N], depth[N], size[N];
void init(int n)
{
        E = 0;
        std::fill(first, first + n + 1, -1);
        std::fill(solved, solved + n + 1, 0);
        pool = node;
        spt_size = 0;
        for(int i = 1; i <= n; i++) {
                edge[i].clear();
        }
}
void add_edge(int a, int b, int c)
{
        pnt[E] = b;
        nxt[E] = first[a];
        weight[E] = c;
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
                                depth[v] = depth[u] + weight[i];
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
        SplayTree *_tree = &spt[spt_size];
        edge[root].push_back(Center(0, _tree));
        _tree->init();
        _tree->insert(0, w[root], root);
        spt_size++;
        for(int i = first[root]; ~i; i = nxt[i]) {
                int v = pnt[i];
                if(solved[v]) {
                        continue;
                }
                depth[v] = weight[i]; 
                int tail = bfs(v, root);
                for(int i = 0; i < tail; i++) {
                        int u = queue[i];
                        edge[u].push_back(Center(depth[u], _tree));
                        _tree->insert(depth[u], w[u], u);
                }
                solve(v);
        }
}
int get_val(){
        int ret(0),sign(1);
        char c;
        while((c=getchar())==' '||c=='\n'||c=='\r');if(c=='-')    sign=-1;
        else    ret=c-'0';
        while((c=getchar())!=' '&&c!='\n'&&c!='\r')    if(c=='-')    sign=-1;
        else    ret=ret*10+c-'0';
        return ret*sign;
}
int main()
{
        null = new Node();
        null->clear();
        null->mx = -1;
        int q, a, b, c;
        while(scanf("%d", &n) == 1) {
                init(n);
                for(int i = 1; i <= n; i++) {
                        //scanf("%d", &w[i]);
                        w[i]=get_val();
                }
                for(int i = 1; i < n; i++) {
                        a=get_val(),b=get_val(),c=get_val();
                        //scanf("%d%d%d", &a, &b, &c);
                        add_edge(a, b, c);
                        add_edge(b, a, c);
                }
                solve(1);
                scanf("%d", &q);
                while(q--) {
                        //scanf("%d%d", &a, &b);
                        a=get_val(),b=get_val();
                        int ret = 0, mx = 0;
                        for(int i = 0; i < (int)edge[a].size(); i++) {
                                Center root = edge[a][i];
                                if(b >= root.d) {
                                        pair<int, int> tmp = root.root->query(b - root.d);
                                        if(tmp.second > mx || tmp.second == mx && tmp.first < ret) {
                                                ret = tmp.first;
                                                mx = tmp.second;
                                        }
                                }
                        }
                        printf("%d\n", ret);
                }
        }
        return 0;
}
