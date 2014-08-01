#include <cstdio>
#include <cstring>
#include <algorithm>
const int N = 300010;
int sa[N],X[N],Y[N],b[N],a[N],h[N],r[N];
bool comp(int *r,int a,int b,int le)
{
        return r[a] == r[b] && r[a+le] == r[b+le];
}
void sort(int *Rank,int *Id,int n,int m)
{
        std::fill(b,b+m,0);
        for(int i = n-1; i >= 0; i--) b[Rank[i]]++;
        for(int i = 1; i < m; i++) b[i] += b[i-1];
        for(int i = n-1; i >= 0; i--) sa[--b[Rank[Id[i]]]] = Id[i];
        //之所以倒过来：如果第一关键字相同，应该把第二关键字大的排在后面
        //这里的Id[]已经是第二关键字从小到大的下标了
}
void calh(int n)
{
        for(int i = 1; i <= n; i++) r[sa[i]] = i;
        int height = 0;
        for(int i = 0; i < n; i++){
                if(height) height--;
                int j = sa[r[i]-1];
                while(a[j+height]==a[i+height]) height++;
                h[r[i]] = height;
        }
}
void suffix(int n,int m=500)
{
        int *Rank = X, *Id = Y, p = 1;
        for(int i = 0; i < n; i++) Rank[i] = a[i], Id[i] = i;
        sort(Rank,Id,n,m);
        for(int j = 1; p < n; j <<= 1){
                p = 0;
                for(int i = n-j; i < n; i++) Id[p++] = i;
                //suffix(n-j) -> suffix(n-1)都没有第二关键字，应该排在前面
                for(int i = 0; i < n; i++) if(sa[i] >= j) Id[p++] = sa[i] - j;
                // Id[] 为第二关键字从小到大的下标 , Id[i]处的第二关键字小于 Id[i+1]处的第二关键字
                // 即s[Id[i]+j] <= s[Id[i+1]+j]
                sort(Rank,Id,n,p);
                std::swap(Rank,Id);
                Rank[sa[0]] = 0, p = 1;
                for(int i = 1; i < n; i++)
                        Rank[sa[i]] = comp(Id,sa[i-1],sa[i],j) ? p-1 : p++;
                m = p;
        }
        calh(n-1);
}
