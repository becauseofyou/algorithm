//太长了，压缩一下
#include <cstdio>
#include <algorithm>
const int N = 400010;
const int inf = ~0u >> 2;
struct Node* null;
struct Node
{
        Node* ch[2];
        Node* fa;
        int mi, val, rev;
        inline int sgn(){return fa->ch[0]==this?0:fa->ch[1]==this?1:-1;}
        inline void setc(int s, Node* who){who->fa=this;ch[s]=who;}
        inline void rotate() {
                int a=sgn(),b=fa->sgn(); 
                Node* y=fa; y->setc(a,ch[!a]),fa=y->fa;
                if(~b) fa->ch[b]=this;
                this->setc(!a,y),y->up();
        }
        inline void splay() {
                go();
                for(int a,b;~(a=sgn());rotate()) 
                        if(~(b=fa->sgn()))
                                (a^b)?rotate():fa->rotate();
                up();
        }

        inline void up() {
                mi = std::min(ch[0]->mi, ch[1]->mi);
                mi = std::min(mi, val);
        }
        void init(int v, Node* f) {
                fa = f;
                val = mi = v;
                ch[0] = ch[1] = null;
                rev = 0;
        }
        void flip() {
                std::swap(ch[0], ch[1]);
                rev ^= 1;
        }
        void push() {
                if(rev) {
                        ch[0]->flip(), ch[1]->flip(), rev = 0;
                }
        }
        void print() {
                if(this != null) {
                        if(ch[0] != null) ch[0]->print();
                        printf("%d\n", val);
                        if(ch[1] != null) ch[1]->print();
                }
        }
        void go() {
                if(~sgn()) fa->go();
                push();
        }
};
struct Link_Cut_Tree 
{
        Node node[N]; 
        Node* access(Node* u) {
                Node* v = null;
                for(; u != null; v = u, u = u->fa) {
                        u->splay();
                        u->ch[1] = v;
                        u->up();
                }
                return v;
        }
        int ask(Node* x, Node* y) {
                access(x);
                for(x = null; y != null; y = y->fa) {
                        y->splay();
                        if(y->fa == null)  {
                                return std::min(y->ch[1]->mi, x->mi);
                        }
                        y->ch[1] = x, y->up(), x = y;
                }
        }
        int ask(int a, int b) {
                return ask(node + a, node + b);
        }
        void cut(int v) {
                Node* u = node + v;
                access(u);
                u->splay();
                u->ch[0]->fa = null;
                u->ch[0] = null;
                u->up();
        }

        //切掉u跟v之间的边
        void cut(Node* u, Node* v) {
                access(u), v->splay();
                if(v->fa == u) {
                        v->fa = null;
                } else {
                        access(v), u->splay(), u->fa = null;
                }
        }
        //连接u跟v之间的边
        void link(Node* u, Node* v) {
                make_root(u);
                u->fa = v;
        }
        //将u变成树根
        void make_root(Node* u) {
                Node* tmp = access(u);
                tmp->flip();
                u->splay();
        }
        bool judge(Node* u, Node* v) {
                while(u->fa != null) u = u->fa;
                while(v->fa != null) v = v->fa;
                return u == v;
        }
}lct;
