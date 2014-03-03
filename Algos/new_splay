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
  Node *c[2], *fa;
  int isRoot, sz, cnt, val;
  int fill;
  void setc(int d,Node *s) {
    c[d] = s;
    s->fa = this;
  }
  bool d() {
    return fa->c[1] == this;
  }
  void up() {
    sz = c[0]->sz + c[1]->sz + cnt;
    fill = c[0]->fill + c[1]->fill + (cnt == K);
  }
  void clear() {
    sz = 1;
    isRoot = fill = 0;
    c[0] = c[1] = fa = null;
  }
}pool[N];

int top;
Node* root[N];

Node* newnode() {
  Node *tmp;
  tmp = &pool[top++];
  tmp->clear();
  return tmp;
}

void rot(Node *x,int f) {
  Node *y = x->fa;
  std::swap(x->isRoot, y->isRoot);
  y->setc(!f, x->c[f]);
  x->fa = y->fa;  
  if(y->fa != null) y->fa->c[y->d()] = x;
  x->setc(f, y);
  y->up();
}

void splay(Node *x) {
  while(!x->isRoot) {
    if(x->fa->isRoot) {
      rot(x, !x->d());
    } else {
      bool f = x->fa->d();
      x->d() == f ? rot(x->fa, !f): rot(x, f);
      rot(x, !f);
    }
  }
  x->up();
}

Node* insert(Node *&x, Node *y, Node *fa) {
  if(x == null) {
    x = y;
    x->fa = fa;
    return x;
  }
  Node *ret;
  if(x->val == y->val) {
    x->cnt += y->cnt;
    ret = x;
  }else if(x->val > y->val) {
    ret = insert(x->c[0], y, x);
  } else{
    ret = insert(x->c[1], y, x);
  }
  x->up();
  return ret;
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
    Node *tmp = insert(A, fr, null);
    splay(tmp);
    A = tmp;
  }
}

int firstEdge[N], nextEdge[N*2], pointTo[N*2];
int totEdge, ans[N], col[N], n;

void addEdge(int a,int b) {
  pointTo[totEdge] = b;
  nextEdge[totEdge] = firstEdge[a];
  firstEdge[a] = totEdge ++;
}

void dfs(int u,int f) {
  for(int iter = firstEdge[u]; iter != -1; iter = nextEdge[iter]) {
    int v = pointTo[iter];
    if(v == f) continue;
    dfs(v, u);
    merge(root[u], root[v]);
  }
  ans[u] = root[u]->fill;
}

void init() {
  top = 0; 
  null = new Node();
  null->sz = 0; null->c[0] = null->c[1] = NULL; null->cnt = 0;
  null->isRoot = 0;
  totEdge = 0;
  std::fill(firstEdge, firstEdge + n, -1);
  for(int i = 0; i < n; i++) {
    scanf("%d",&col[i]);
  }
  for(int i = 0; i < n; i++) { 
    root[i] = newnode();
    root[i]->cnt = 1;
    root[i]->val = col[i];
    root[i]->isRoot = 1;
    root[i]->up();
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
