#pragma comment(linker, "/STACK:1000000000,1000000000")
#include <cstdio>
#include <queue>
#include <algorithm>
#include <cstring>
const int N = 100010;

struct Node;
Node* null;
int K;

struct Node {
    Node *c[2], *p;
    Node() {
        c[0] = c[1] = p = null;
        isRoot = sz = cnt = val = 0;
    }
    Node(int _isRoot, int _sz, int _cnt, int _val) {
        c[0] = c[1] = p = null;
        isRoot = _isRoot;
        sz = _sz;
        cnt = _cnt;
        val = _val;
    }
    int isRoot, sz, cnt, val;
    int fill;
    inline void setc(int d,Node *s) {
        c[d] = s;
        s->p = this;
    }
    inline bool d() {
        return p->c[1] == this;
    }
    void up() {
        sz = c[0]->sz + c[1]->sz + cnt;
        fill = c[0]->fill + c[1]->fill + (cnt == K);
    }
    void clear() {
        sz = 1;
        isRoot = fill = 0;
        c[0] = c[1] = p = null;
    }
}pool[N], *C ;

Node* root[N];

/*
 * 从高到低依次修改
 */
void rot(Node *x) {
    Node *p = x->p;
    int f = x->d();
    std::swap(x->isRoot, p->isRoot);
    p->p->setc(p->d(), x);
    p->setc(f, x->c[!f]);
    x->setc(!f, p);
    p->up();
}

void splay(Node *x) {
    while(!x->isRoot) {
        if(x->p->isRoot) {
            rot(x);
        } else {
            bool f = x->p->d();
            x->d() == f ? rot(x->p): rot(x);
            rot(x);
        }
    }
    x->up();
}

Node* insert(Node *x, Node *y) { //将y节点插入某棵splay树
    if(x->val == y->val) {
        x->cnt += y->cnt;
        return x;
    }
    int d = y->val > x->val;
    if(x->c[d] == null) {
        x->setc(d, y);
        return y;
    }
    return insert(x->c[d], y);
}
void debug(Node *x) {
    if(x->c[0] != null) debug(x->c[0]);
    printf("col=%d %d %d\n",x,x->val,x->cnt);
    if(x->c[1] != null) debug(x->c[1]);
}

void merge(Node *&A, Node *B) { // merge tree_B to  tree_A
    if(A->sz <= B->sz) {
        std::swap(A, B);
    }
    std::queue<Node*> Q;
    Q.push(B);
    while(!Q.empty()) {
        Node *fr = Q.front(); Q.pop();
        if(fr->c[0] != null) Q.push(fr->c[0]);
        if(fr->c[1] != null) Q.push(fr->c[1]);
        fr->clear();
        Node *tmp = insert(A, fr);
        splay(tmp);
        A = tmp;
    }
}


struct Edge {
    int v;
    Edge *next;
}ee[N*2], *E;
Edge *list[N];
void addEdge(int a, int b) {
    Edge *e = new(E++)Edge();
    e->v = b;
    e->next = list[a];
    list[a] = e;
}
int ans[N], col[N], n;

void dfs(int u,int f) {
    for(Edge *e = list[u]; e; e = e->next) {
        int v = e->v;
        if(v == f) continue;
        dfs(v, u);
        merge(root[u], root[v]);
    }
    ans[u] = root[u]->fill;
}

void init() {
    C = pool; E = ee;
    null = new (C++)Node();
    for(int i = 0; i < n; i++) {
        list[i] = NULL;
        scanf("%d",&col[i]);
    }
    for(int i = 0; i < n; i++) { 
        root[i] = new (C++) Node(1, 1, 1, col[i]);
        root[i]->up();
        //        debug(root[i]);
    }
}

int main() {
    int t, ca = 1;
    scanf("%d",&t);
    while(t--) {
        scanf("%d%d",&n,&K);
        init();
        for(int i = 0, a, b; i < n - 1; i++) {
            scanf("%d%d",&a,&b);
            a --; b --;
            addEdge(a, b);
            addEdge(b, a);
        }
        dfs(0, -1);
        int q, a;
        if(ca > 1) puts("");
        printf("Case #%d:\n",ca++);
        scanf("%d",&q);
        while(q--) {
            scanf("%d",&a);
            a--;
            printf("%d\n",ans[a]);
        }
    }
    return 0;
}
