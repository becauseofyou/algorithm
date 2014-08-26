struct Mancher {  
    char str[maxn];//start from index 1  
    int p[maxn];  
    char s[maxn];  
    int n;  
    void checkmax(int &ans,int b){  
        if(b>ans) ans=b;  
    }  
    inline int min(int a,int b){  
        return a<b?a:b;  
    }  
    void kp(){  
        int i;  
        int mx = 0;  
        int id;  
        for(i=1; i<n; i++){  
            if( mx > i )  
                p[i] = min( p[2*id-i], p[id]+id-i );  
            else  
                p[i] = 1;  
            for(; str[i+p[i]] == str[i-p[i]]; p[i]++) ;  
            if( p[i] + i > mx ) {  
                mx = p[i] + i;  
                id = i;  
            }  
        }  
    }  
    void pre()  
    {  
        int i,j,k;  
        n = strlen(s);  
        str[0] = '$';  
        str[1] = '#';  
        for(i=0;i<n;i++)  
        {  
            str[i*2 + 2] = s[i];  
            str[i*2 + 3] = '#';  
        }  
        n = n*2 + 2;  
        str[n] = 0;  
    }  
    void solve() // 求出所有的最长回文子串所在的区间  
    {  
        int & tot = M::n;  
        tot = 0;  
        for(int i = 2; i < n; i++)  
        {  
            if(i%2&&p[i]==1) continue;  
            if(i%2)  
            {  
                M::in[tot++] = M::node(i/2-p[i]/2+1,i/2+p[i]/2);  
                //printf("%d %d\n",i/2-p[i]/2+1,i/2+p[i]/2);  
            }  
            else   
            {  
                M::in[tot++] = M::node(i/2-(p[i]/2-1),i/2+(p[i]/2-1));  
                //printf("%d %d\n",i/2-(p[i]/2-1),i/2+(p[i]/2-1));  
            }  
        }  
    }  
}task1;  
