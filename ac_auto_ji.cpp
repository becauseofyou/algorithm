const int M = 240000;
const int CD = 26;
int ch[M][CD];
int val[M];
int fail[M];
int Q[M];
int used[M];
int sz;
int ID[128];
void init(){
    fail[0]=0;
    for(int i=0;i<CD;i++) ID[i+'a']=i;
}
void Reset(){
    memset(ch[0],0,sizeof(ch));
    sz=1;
}
void Insert(char *a){
    int p=0;
    for(;*a;a++)  {
        int c=ID[*a];
        if(!ch[p][c])  {
            memset(ch[sz],0,sizeof(ch[sz]));
            val[sz]=0;
            used[sz]=false;
            ch[p][c]=sz++;
        }
        p=ch[p][c];
    }
    val[p]++;
}
void Construct(){
    int *s=Q,*e=Q;
    for(int i=0;i<CD;i++) {
        if(ch[0][i]){
            fail[ ch[0][i] ] = 0;
            *e ++ = ch[0][i];
        }
    }while(s!=e){
        int u = *s++;
        for(int i = 0; i < CD ;i++){
            int &v = ch[u][i];
            if(v){
                *e++ = v;
                fail[v] = ch[ fail[u] ][i];
            } else  {
                v=ch[ fail[u] ][i];
            }
        }
    }
}
int AC(char *s){
    int p=0,ans=0;
    for(;*s;s++){
       int cur=ch[p][ID[*s]];
       while(cur && !used[cur]){
           used[cur]=true;
           if(val[cur]) ans+=val[cur];
           cur=fail[cur];
       }
       p=ch[p][ID[*s]];
    }
    return ans;
}
