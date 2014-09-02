struct Node;
Node* null;

struct Node {
    Node *c[2], *p;
    Node() {
        c[0] = c[1] = p = null;
        isRoot = sz;
    }
    Node(int _isRoot, int _sz, int _cnt, int _val) {
        c[0] = c[1] = p = null;
        isRoot = _isRoot;
        sz = _sz;
    }
    int isRoot, sz;
    int val;
    inline void setc(int d,Node *s) {
        c[d] = s;
        s->p = this;
    }
    inline bool d() {
        return p->c[1] == this;
    }
    void up() {
        sz = c[0]->sz + c[1]->sz + 1;
    }
    void clear() {
        sz = 1;
        isRoot = 0;
        c[0] = c[1] = p = null;
    }
}pool[N], *C ;

Node* idx[N];
Node* root;
void rot(Node *x) {
    Node *p = x->p;
    int f = x->d();
    std::swap(x->isRoot, p->isRoot);
    p->p->setc(p->d(), x);
    p->setc(f, x->c[!f]);
    x->setc(!f, p);
    p->up();
}

void splay(Node *x, Node *goal) {
    while(x->p != goal) {
        if(x->p->p == goal) {
            rot(x);
        } else {
            bool f = x->p->d();
            x->d() == f ? rot(x->p): rot(x);
            rot(x);
        }
    }
    x->up();
    if(goal == null) root = x;
}
void slect(int k, Node *goal) {
    Node *x = root;
    while(x->c[0]->sz + 1 != k) {
        if(k < x->c[0]->sz + 1) {
            x = x->c[0];
        } else {
            k -= x->c[0]->sz + 1;
            x = x->c[1];
        }
    }
    splay(x, goal);
}
int a[N];
void build(Node *&x, int l, int r, Node *fa) {
    if(l > r)  {
        return ;
    }
    int m = (l + r) >> 1;
    x = new(C++) Node();
    x->p = fa;
    x->sz = 1;
    idx[a[m]] = x;
    x->val = a[m];
    build(x->c[0], l, m - 1, x);
    build(x->c[1], m + 1, r, x);
    x->up();
}
void debug(Node *x) {

    if(x->c[0] != null) debug(x->c[0]);
    printf("%d ", x->val);
    if(x->c[1] != null) debug(x->c[1]);
}
void init() {
    res.clear();
    C = pool;
    null = new (C++) Node();
    root = new(C++) Node();
    root->isRoot = 1;
    root->setc(1, new(C++) Node());
    build(root->c[1]->c[0], 1, n, root->c[1]);
    root->c[1]->up(); root->up();
}

// cut the number between L and R, and insert them after position p
void cut_insert(int L, int R, int p) {
    slect(L - 1, null);
    slect(R + 1, root);
    Node *stree = root->c[1]->c[0];
    root->c[1]->c[0] = null;
    slect(p, null);
    slect(p + 1, root);
    root->c[1]->setc(0, stree);
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

void delRoot() {
    Node* t = root;
    if(t->c[1]!=null) {
        root = t->c[1];
        slect(1, null);
        root->setc(0,t->c[0]);
    } else root = root->c[0];
    root->fa = null;
    if(root!=null) root->up();
}





#include <cstdio>
#include <algorithm>
const int N = 100110;
 
struct Node;
Node* null;
 
struct Node 
{
        Node* c[2];
        Node* fa;
        Node() {
                c[0] = c[1] = fa = null;
        }
        Node* link;
        int sz, val;
        inline void setc(int d, Node *s) {
                c[d] = s;
                s->fa = this;
        }
        inline bool d() {
                return fa->c[1] == this;
        }
        void up() {
                sz = c[0]->sz + c[1]->sz + 1;
        }
        void clear() {
                sz = 1;
                c[0] = c[1] = fa = null;
        }
};
 
void rot(Node* x) 
{
        Node* p = x->fa;
        int f = x->d();
        p->fa->setc(p->d(), x);
        p->setc(f, x->c[!f]);
        x->setc(!f, p);
        p->up();
}
 
void splay(Node* &root, Node* x, Node* goal) 
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
 
void slect(Node* &root, int k, Node* goal) 
{
        Node* x = root;
        while(x->c[0]->sz + 1 != k) {
                if(k < x->c[0]->sz + 1) {
                        x = x->c[0];
                } else {
                        k -= x->c[0]->sz + 1;
                        x = x->c[1];
                }
        }
        splay(root, x, goal);
}
 
int get_rank(Node* &root, Node* x)
{
        splay(root, x, null);
        return x->c[0]->sz + 1;
}
 
Node pool[2 * N], *tail = pool;
Node* root;      // root ： 位置信息的平衡树
Node* value[N]; // 每个值的平衡树，保存位置信息平衡树中的相应节点
 
void show(Node* root)
{
        if(root == null) {
                return ;
        }
        show(root->c[0]);
        printf("rank=%d\n", get_rank(::root, root->link));
        show(root->c[1]);
}
 
void debug(Node *x) 
{
        if(x->c[0] != null) debug(x->c[0]);
        printf("x = %d ", x->val);
        if(x->c[1] != null) debug(x->c[1]);
}
 
void erase(Node* &root, Node* x) 
{
        splay(root, x, null);
        Node *t = root;
        if(t->c[1] != null) {
                root = t->c[1];
                slect(root, 1, null);
                root->setc(0, t->c[0]);
        } else {
                root = root->c[0];
        }
        root->fa = null;
        if(root != null) {
                root->up();
        }
}
 
Node* seq_find(int pos)
{
        Node* p = root;
        while(true) {
                if(p->c[0]->sz + 1 < pos) {
                        pos -= p->c[0]->sz + 1;
                        p = p->c[1];
                } else if(p->c[0]->sz + 1 == pos) {
                        return p;
                } else {
                        p = p->c[0];
                }
        }
}
 
void insert(Node* &root, Node* x)
{
        Node* now = root;
        Node* pre = root;
        int c = -1;
        while(now != null) {
                pre = now;
                if(get_rank(::root, x->link) < get_rank(::root, now->link)) {
                        now = now->c[0], c = 0;
                } else {
                        now = now->c[1], c = 1;
                }
        }
        x->clear();
        pre->setc(c, x);
        splay(root, x, null);
}
 
void shift(int l, int r)
{
        if(l == r) {
                return ;
        }
        Node* x = seq_find(l);
        Node* y = seq_find(r);
        erase(root, y);
        y->clear();
        splay(root, x, null);
        if(x->c[0] == null) {
                x->setc(0, y);
                x->up();
        } else {
                Node* prev = x->c[0];
                while(prev->c[1] != null) {
                        prev = prev->c[1];
                }
                splay(root, prev, x);
                prev->setc(1, y);
                prev->up(); x->up();
        }
        erase(value[y->val], y->link);
        insert(value[y->val], y->link);
}
 
int get_cnt(Node* &root, int pos)
{
        Node* p = root;
        Node* q = null;
        int cnt = 0;
        while(p != null) {
                q = p;
                if(pos < get_rank(::root, p->link)) {
                        p = p->c[0];
                } else {
                        cnt += p->c[0]->sz + 1;
                        p = p->c[1];
                }
        }
        if(q != null) {
                splay(root, q, null);
        }
        return cnt;
}
 
void init(int n)
{
        null = tail++;
        null->sz = 0; null->c[0] = null->c[1] = null->fa = null;
        null->val = 0;
        root = null;
        for(int i = 1; i <= n; i++) {
                value[i] = null;
        }
}
int main()
{
        int n;
        static int a[N];
        scanf("%d", &n);
        init(n);
        for(int i = 1; i <= n; i++) {
                scanf("%d", &a[i]);
        }
        for(int i = 1; i <= n; i++) {
                Node* p = tail++;
                Node* q = tail++;
                p->clear(); q->clear();
                p->setc(0, root);
                p->val = a[i];
                p->up();
                root = p;
                q->link = p, p->link = q;
                insert(value[a[i]], q);
        }
        int q, last = 0;
        scanf("%d", &q);
        while(q--) {
                int op, l, r, k; scanf("%d%d%d", &op, &l, &r);
                l = (l + last - 1) % n + 1;
                r = (r + last - 1) % n + 1;
                if(l > r) std::swap(l, r);
                if(op == 1) {
                        shift(l, r);
                } else {
                        scanf("%d", &k);
                        k = (k + last - 1) % n + 1;
                        last = get_cnt(value[k], r) - get_cnt(value[k], l - 1);
                        printf("%d\n", last);
                }
        }
}
