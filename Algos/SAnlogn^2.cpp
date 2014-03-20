#include <cstdio>
#include <cstring>
#include <algorithm>
const int N = 50010;
int sa[N],r[N],h[N],a[N],Rank[N];
char s[N];
int len;
bool cmpsa(int i,int j){
  return r[i] == r[j] ? r[i+len] < r[j+len]: r[i] < r[j];
}
void calh(int n){	
  int height = 0;
  for(int i = 0; i < n; i++){
    if(height) height--;
    int j = sa[r[i]-1];
    while(a[j+height] == a[i+height]) height++;
    h[r[i]] = height;
  }
}
void suffix(int n,int m = 128){
  int i;
  for(i = 0; i < n; i++) r[i] = a[i], sa[i] = i;
  for(len = 1; len < n; len <<=1){
    std::sort(sa,sa+n,cmpsa);
    Rank[sa[0]] = 0; int tot = 1;
    for(i = 1; i < n; i++)
      Rank[sa[i]] = cmpsa(sa[i-1],sa[i]) ? tot++ : tot-1;
    for(i = 0; i < n; i++) r[i] = Rank[i];
    if(tot == n) break;
  }
  calh(n-1);
}
