
const int N = 150 * 70 + 10;
const int CD = 26;
//template based
int sz, ch[N][CD], val[N], fail[N], Q[N], ID[512];


//problem needed
char s[1000010];
int flag[N], id[N];

void Init() {
    fail[0] = 0;
    for (int i = 0; i < 26; i++) {
        ID[i + 'a'] = i;
    }
}

void Reset() {
    memset (ch[0], 0, sizeof(ch));
    sz = 1;
}

void Insert(char *a, int idx) {
    int p = 0;
    for(; *a; a++)  {
        int c = ID[*a];
        if (!ch[p][c])  {
            memset (ch[sz], 0, sizeof(ch[sz]));
            val[sz] = 0;
            flag[sz] = 0;
            ch[p][c] = sz++;
        }
        p = ch[p][c];
    }
    val[p] = 1;
    id[p] = idx;
}

void Construct() {
    int *s = Q, *e = Q;
    for (int i = 0; i < CD; i++) {
        if (ch[0][i]) {
            fail[ ch[0][i] ] = 0;
            *e ++ = ch[0][i];
        }
    }
    while (s != e) {
        int u = *s++;
        for (int i = 0; i < CD ;i++){
            int &v = ch[u][i];
            if(v) {
                *e++ = v;
                fail[v] = ch[ fail[u] ][i];
            } else  {
                v = ch[ fail[u] ][i];
            }
        }
    }
}

