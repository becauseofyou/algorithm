struct Graph {
        int head[N];
        int pnt[M];
        int nxt[M];
        int w[M];
        int E;
        int son[N];
        int mx[N];
        int node[N];
        int n;
        void init(int n){
                this->n = n;
                E = 0;
                std::fill(head, head + n + 1, -1);
        }
        void add_edge(int a, int b, int c) {
                pnt[E] = b;
                w[E] = c;
                nxt[E] = head[a];
                head[a] = E++;
        }
        int tot, cnt;
        int get_root(int u) {
                tot = 0;
                dfs(u, -1);
                int ans = -1, tmp, opt = N;
                for(int i = 0; i < tot; i++) {
                     //   printf("hehe %d\n", mx[node[i]]);
                        tmp = tot - son[node[i]] - 1;
                        if(tmp < mx[node[i]]) {
                                tmp = mx[node[i]];
                        }
                        if(tmp < opt) {
                                opt = tmp;
                                ans = node[i];
                        }
                }
                return ans;
        }
        void dfs(int u, int f) {
                node[tot++] = u;
                son[u] = 1;
                mx[u] = 0;
                for(int i = head[u]; i != -1; i = nxt[i]) if(pnt[i] - f) {
                        if(solved[pnt[i]]) {
                                continue;
                        }
                        dfs(pnt[i], u);
                        son[u] += son[pnt[i]];
                        if(son[pnt[i]] > mx[u]) {
                                mx[u] = son[pnt[i]];
                        }
                }
        }
        pii p[N];
        void get_dist(int u, int f, int dis, int d) {
                p[cnt++] = mp(d, dis);
                for(int i = head[u]; i != -1; i = nxt[i]) if(pnt[i] - f) {
                        int to = pnt[i];
                        if(solved[to]) continue;
                        get_dist(to, u, dis + w[i], d + 1);
                }
        }
        int map[N];
        int count[N];
        bool solved[N];
        int used[N];
        int tim;
        int cc;
        int dd[N];
        void solve(int u, int f) {
                u = get_root(u);
                ++tim;
            //    printf("u = %d\n", u);
                cc = 0;
                for(int i = head[u]; i != -1; i = nxt[i]) if(pnt[i] - f){
                        if(solved[pnt[i]]) continue;
                 //       printf("to = %d\n", pnt[i]);
                        cnt = 0;
                        get_dist(pnt[i], u, w[i], 1);
                        update();
                }
                for(int i = 0; i < cc ; i++) {
                    map[dd[i]] = 0;
                    count[dd[i]] = 0;
                }
                //printf("ret=%d count=%I64d\n", ret, ret_count);
                solved[u] = true;
                for(int i = head[u]; i != -1; i = nxt[i]) if(pnt[i] - f){
                        if(!solved[pnt[i]]) {
                                solve(pnt[i], u);
                        }
                }
        }
        void update() {
                for(int j = 0; j < cnt; j++) {
                       // printf("%d %d\n", p[j].first, p[j].second);
                        if(K > (p[j].first + 1)) {
                                int len = K - p[j].first - 1;
                                int tmp = p[j].second + map[len];
                                if(tmp > ret) {
                                        ret = tmp;
                                        ret_count = count[len];
                                } else if(tmp == ret) {
                                        ret_count += count[len];
                                }
                        }
                        if(p[j].first == K - 1) {
                                if(p[j].second > ret) {
                                        ret = p[j].second;
                                        ret_count = 1;
                                } else if(p[j].second == ret) {
                                        ret_count += 1;
                                }
                        }
                }
                for(int j = 0; j < cnt; j++) {
                        if(used[p[j].first] != tim) {
                            used[p[j].first] = tim;
                            dd[cc++] = p[j].first;
                        }
                        if(p[j].second > map[p[j].first]) {
                                map[p[j].first] = p[j].second;
                                count[p[j].first] = 1;
                        } else if(p[j].second == map[p[j].first]) {
                                count[p[j].first] += 1;
                        }
                }
        }
        int ret;
        int ret_count;
        void ac() {
                tim = 0;
                std::fill(used, used + n + 1,0);
                std::fill(map, map + n + 1, 0);
                std::fill(count, count + n + 1, 0);
                std::fill(solved, solved + n + 1, false);
                ret = 0;
                ret_count = 0;
                solve(1, -1);
                printf("%d %d\n", ret, ret_count);
        }
}
