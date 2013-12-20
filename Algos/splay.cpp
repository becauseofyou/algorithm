struct node 
{
  node *c[2],*fa;
  int val;
  int sz;
  int belong;
  void setc(int d,node *s) {
    c[d]  = s;
    s->fa = this;
  }
  bool d() {
    return fa->c[1]==this;
  }
  void up(){
    sz = c[0]->sz + c[1]->sz + 1;
  }
  void clear(node *null) {
    sz = 1; 
    c[0] = c[1] = fa = null;
  }
}NODE[N],*null=NODE;
node* ID[N];
int top;
struct Tree
{
  node* root;
  void init(int v) {
    node* x = &NODE[++top];
    x->belong = top;
    x->val = v;
    x->clear(null);
    root = x;
  }
  void rot(node* x,int f) {
    node* y = x->fa;
    y->setc(!f,x->c[f]);
    x->fa = y->fa;
    if(y->fa!=null) y->fa->c[y->d()] = x;
    x->setc(f,y);
    y->up();
  }
  void splay(node *x,node *goal) {
    while(x->fa!=goal) {
      if(x->fa->fa == goal) rot(x,!x->d());
      else {
        bool f = x->fa->d();
        x->d() == f ? rot(x->fa,!f): rot(x,f);
        rot(x,!f);
      }
    }
    x->up();
    if(goal == null) root = x;
  }
  void RTO(int k,node *goal) {
    node *x = root;
    while(x->c[0]->sz+1 != k) {
      if(x->c[0]->sz+1 > k) x = x->c[0];
      else {
        k -= x->c[0]->sz + 1;
        x = x->c[1];
      }
    }
    splay(x,goal);
  }
  void delRoot() {
    node* t = root;
    if(t->c[1]!=null) {
      root = t->c[1];
      RTO(1,null);
      root->setc(0,t->c[0]);
    } else root = root->c[0];
    root->fa = null;
    if(root!=null) root->up();
  }
  void change(node* x,int v){
    splay(x,null);
    delRoot();
    x->val = v; 
    x->clear(null);
    insert(root,x);
  }
  void insert(node* &x,node *y,node* fa=null) {// 如何减少对父亲指针的赋值
    if(x == null) {
      x = y;
      y->fa = fa;
      return ;
    }
    if(x->val >= y->val) insert(x->c[0],y,x);
    else insert(x->c[1],y,x);
    x->up();
  }
  int find(node *x,int k) {
    if(x->c[0]->sz + 1 == k) return x->val;
    else if(x->c[0]->sz + 1 > k) return find(x->c[0],k);
    else return find(x->c[1],k-x->c[0]->sz-1);
  }
  void debug(node *x) {
    printf("%d lc:%d rc:%d\n",x->val,x->c[0]->val,x->c[1]->val);
    if(x->c[0]!=null) debug(x->c[0]);
    if(x->c[1]!=null) debug(x->c[1]);
  }
};
