//http://www.spoj.com/problems/QTREE/
#include <cstdio>
#include <algorithm>
const int N = 10010;
const int inf = ~0u >> 2;
struct Node* null;
struct Node
{
        Node* ch[2];
        Node* fa;
        int mx;
        int val;
        inline bool isroot() {
                return fa == null || fa->ch[0] != this && fa->ch[1] != this;
        }
        inline bool d() {
                return fa->ch[1] == this;
        }
        inline void setc(int s, Node* who) {
                who->fa = this;
                ch[s] = who;
        }
        inline void rotate() {
                Node* y = fa;
                bool f = d();
                y->setc(f, ch[!f]);
                fa = y->fa;
                if(!y->isroot())
                        fa->ch[y->d()] = this;
                this->setc(!f, y), y->up();
        }
        void splay() {
                for(; !isroot(); ) {
                        if(!fa->isroot()) (fa->d() ^ d()) ? rotate() : fa->rotate();
                        rotate();
                }
                up();
        }
        inline void up() {
                mx = std::max(ch[0]->mx, ch[1]->mx);
                mx = std::max(mx, val);
        }
        void init(int v, Node* f) {
                fa = f;
                val = mx = v;
                ch[0] = ch[1] = null;
        }
        void print() {
                if(this != null) {
                        if(ch[0] != null) ch[0]->print();
                        printf("%d\n", val);
                        if(ch[1] != null) ch[1]->print();
                }
        }
};
struct Link_Cut_Tree 
{
        Node node[N];
        Node* access(Node* u) {
                Node* v = null;
                for(; u != null; u = u->fa) {
                        u->splay();
                        u->ch[1] = v;
                        u->up();
                        v = u;
                }
                return v;
        }
        int ask(Node* x, Node* y) {
                access(x);
                for(x = null; y != null; y = y->fa) {
                        y->splay();
                        if(y->fa == null) {
                                return std::max(y->ch[1]->mx, x->mx);
                        }
                        y->ch[1] = x, y->up(), x = y;
                }
        }
        int ask(int a, int b) {
                return ask(node + a, node + b);
        }
        void change(int u, int v) {
                change(node + u, v);
        }
        void change(Node* x, int v) {
                x->val = v;
                x->splay();
        }
}lct;
int head[N], nxt[2 * N], pnt[2 * N], E, id[2 * N], mp[N], wi[N * 2];
void add(int a, int b, int c, int idx)
{
        id[E] = idx;
        pnt[E] = b;
        nxt[E] = head[a];
        wi[E] = c;
        head[a] = E++;
}
void dfs(int u, int fa)
{
        if(fa == -1) lct.node[u].init(-inf, null);
        for(int i = head[u]; i != -1; i = nxt[i]) {
                int v = pnt[i], idx = id[i], w = wi[i];
                if(v == fa) {
                        continue;
                }
                lct.node[v].init(w, lct.node + u);
                mp[idx] = v;
                dfs(v, u);
        }
}
int main()
{
        null = new Node();
        null->ch[0] = null->ch[1] = null->fa = null;
        null->mx = null->val = -inf;
        int n, t, a, b, c;
        scanf("%d", &t);
        while(t--) {
                scanf("%d", &n);
                E = 0;
                std::fill(head, head + n + 1, -1);
                for(int i = 1; i < n; i++) {
                        scanf("%d%d%d", &a, &b, &c);
                        add(a, b, c, i);
                        add(b, a, c, i);
                }
                dfs(1, -1);
                char op[10];
                while(scanf("%s", op) == 1) {
                        if(op[0] == 'D') {
                                break;
                        }
                        if(op[0] == 'Q') {
                                scanf("%d%d", &a, &b);
                                printf("%d\n", lct.ask(a, b));
                        } else {
                                scanf("%d%d", &a, &b);
                                lct.change(mp[a], b);
                        }
                }
        }
        return 0;
}
