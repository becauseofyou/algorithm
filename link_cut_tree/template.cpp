#include <bits/stdc++.h>
using namespace std;
#define lc c[0]
#define rc c[1]
#define me this
const int N = 10010;
const int inf = ~0u >> 2;
struct Data {
    int maximum;
    int value;
    Data(int maximum = -inf, int value = -inf) : maximum(maximum), value(value) {
    }
};
Data operator + (Data a, Data b) {
    int mx = max(a.maximum, b.maximum);
    mx = max(mx, a.value);
    return Data(mx, a.value);
}
struct Node* null;
struct Node {
    Node* c[2], *p;
    Data data;
    int idx;
    Node(int _value = 0, Node* _p = null, int _idx = 0): p(_p), idx(_idx), data(_value, _value){
        lc = rc = null;
    }
    int S(){return p->lc==me?0:p->rc==me?1:-1;}
    void T(int s, Node* w){w->p=me;c[s]=w;}
    void R() {
        int a=S(),b=p->S(); 
        Node* y=p; 
        y->T(a,c[!a]),p=y->p;
        p->c[b]=(~b)?me:null;
        T(!a,y),y->U();
    }
    void F() {}
    void Splay() {
        int a,b;
        while(~(a=S())) {
            b=p->S();
            (~b) ? ((a^b)?R():p->R()) : F();
            R();
        }
        U();
    }
    void U();
    void Print() ;
};
void Node::U() {
    data.maximum = data.value;
    data = data + lc->data;
    data = data + rc->data;
}
void Node::Print() {
    if(lc != null) {
        lc->Print();
    }
    printf("idx=%d\n", idx);
    if(rc != null) {
        rc->Print();
    }
}

Node node[N];
Node* Access(Node *u) { // return the splay tree associate with root
    Node *v = null;
    while(u != null) {
        u->Splay();
        u->rc = v, u->U(), v = u, u = u->p;
    }
    return v;
}
int Ask(Node *u, Node *v) {
    Access(u);
    Node *pre = null;
    while(v != null) {
        v->Splay();
        if (v->p == null) { //  merge the information under lca 
            return max(v->rc->data.maximum, pre->data.maximum);
        }
        v->rc = pre, v->U(), pre = v, v = v->p;
    }
}
void Change(Node *x, int value) {
    x->data.value = value;
    x->Splay();
}
int Ask(int u, int v) {
    return Ask(node + u, node + v);
}
void Change(int u, int value) {
    Change(node + u, value);
}
vector < pair<int, int> > edge[N]; 
int weight[N], point_to[N];
bool used[N];
int q[N];

bool read (int &x) {
        int c = getchar (); int sign = 1;
        while (~c && c < '0' || c > '9') { if (c == '-') sign = -1; c = getchar (); }
        for (x = 0; ~c && '0' <= c && c <= '9'; c = getchar ()) x = x * 10 + c - '0';
        x *= sign;
        return ~c;
}

int main () {
    null = new Node(-inf, null, -1);
    int t, a, b, c, n;
    scanf("%d", &t);
    while (t--) {
        scanf("%d", &n);
        for (int i = 1; i <= n; i++) {
            edge[i].clear();
            used[i] = false;
        }
        for (int i = 0; i < n - 1; i++) {
            //scanf("%d%d%d", &a, &b, &c);
            read(a),read(b),read(c);
            edge[a].emplace_back(make_pair(b, i + 1));
            edge[b].emplace_back(make_pair(a, i + 1));
            weight[i + 1] = c;
        }
        used[1] = true; node[1] = Node(-inf, null, 1);
        int *s = q, *e = q;
        *e++ = 1;
        while(s != e) {
            int u = *s++; 
            for (pair<int, int> it : edge[u]) {
                int v = it.first;
                if(used[v]) {
                    continue;
                }
                used[v] = true; 
                *e++ = v;
                int id = it.second;
                int w = weight[id];
                point_to[id] = v;
                node[v] = Node(w, node + u, v);
            }
        }
        char op[10];
        while(scanf("%s", op) == 1) {
            if(strcmp(op, "DONE") == 0) {
                break;
            }
            read(a), read(b);
            //scanf("%d%d", &a, &b);
            if(op[0] == 'Q') {
                printf("%d\n", Ask(a, b));
            } else {
                Change(point_to[a], b);
            }
            /*
            for (int i = 1; i <= n; i++) { // print every splay tree to debug
                (node + i)->Splay();
                (node + i)->Print();
                puts("");
            }
            puts("ended");
            */
        }
    }
    return 0;
}