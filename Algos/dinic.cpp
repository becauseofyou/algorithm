#include <cstdio>
#include <algorithm>

const int INF = ~0u >> 2;
const int MAX_N = 20010;
const int MAX_E = 500010;

int head[MAX_N], level[MAX_N], cur[MAX_N], pre[MAX_N];
int nxt[MAX_E], pnt[MAX_E], cap[MAX_E], E;
int n;

void add_edge(int from,int to,int c) {
	pnt[E] = to;
	cap[E] = c;
	nxt[E] = head[from];
	head[from] = E++;
}
 
int Q[MAX_N], sign;
bool bfs(int s,int t){
	sign = t;
	std::fill(level, level + n, -1);
	int fr = 0, ed = 0;
	Q[ed++] = ed; level[s] = 0;
	while(fr < ed && level[t] == -1) {
		int u = Q[fr++];
		for(int i = head[u]; i != -1; i = nxt[i]) {
			if(cap[i] > 0 && level[pnt[i]] < 0) {
				level[pnt[i]] = level[u] + 1;
				Q[ed++] = pnt[i];
			}
		}
	}
	return level[t] != -1;
}

void push(int t,int &flow) {
	int mi = INF, p = pre[t];
	for(int p = pre[t]; p != -1; p = pre[pnt[p^1]]) {
		mi = std::min(mi,cap[p]);
	}
	for(int p = pre[t]; p != -1; p = pre[pnt[p^1]]) {
		cap[p] -= mi;
		if(!cap[p]) sign = pnt[p^1];
		cap[p^1] += mi;
	}
	flow += mi;
}

void dfs(int u,int t, int &flow) {
	if(u == t) {
		push(t,flow);
		return ;
	}
	for(int &i = cur[u]; i != -1; i = nxt[i]) {
		if(cap[i] > 0 && level[u] < level[pnt[i]]) {
			pre[pnt[i]] = i;
			dfs(pnt[i],t,flow);
			if(level[sign] < level[u]) return ;
			sign = t;
		}
	}
}

int dinic(int s,int t) {
	pre[s] = -1;
	int flow = 0;
	while(bfs(s,t)) {
	  std::copy(head,head+n,cur);
		dfs(s,t,flow);
	}
	return flow;
}

int main() {
	int m,a,b,w;
	while(scanf("%d%d",&n,&m)!=EOF) {
		int s = 0, t = n+1;
		std::fill(head,head+t+1,-1);
		E = 0;
		for(int i = 1; i <= n; i++) {
			scanf("%d%d",&a,&b);
			add_edge(s,i,a); add_edge(i,s,a);
			add_edge(i,t,b); add_edge(t,i,b);
		}
		for(int i = 0; i < m; i++) {
			scanf("%d%d%d",&a,&b,&w);
			add_edge(a,b,w); 
			add_edge(b,a,w);
		}
		n = t + 1;
		printf("%d\n",dinic(s,t));
	}
	return 0;
}
