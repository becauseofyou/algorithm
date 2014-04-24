#include <cstdio>
#include <ctime>
#include <cstring>
#include <algorithm>
const int N = 10000010;
int a[N];
void quick_sort(int a[],int left,int right)  {
	int mid = (left + right) / 2;
	int midValue = a[mid];
	int i = left, j = right;
	do {
		while(a[i] < midValue) i++;
		while(a[j] > midValue) j--;
		if(i <= j) {
			std::swap(a[i], a[j]);
			i++; j--;
		}
	}while(i <= j);
	if(left < j)  quick_sort(a, left, j);
	if(i < right) quick_sort(a, i, right);
}
int tmp[N];
void merge_sort(int a[],int left, int right) {
	if(left == right) return ;
	int mid = (left + right) / 2;
	merge_sort(a, left, mid);
	merge_sort(a, mid + 1, right);
	int i = left, j = mid + 1, tot = 0;
	while(i <= mid || j <= right)  {
		if(i <= mid) {
			if(j > right || a[i] <= a[j] ){
				tmp[tot++] = a[i];
				i++;
			}
		}
		if(j <= right) {
			if(i > mid || a[j] <= a[i]) {
				tmp[tot++] = a[j];
				j++;
			}
		}
	}
	for(int i = left; i <= right; i++) a[i] = tmp[i-left];
}
const double C = CLOCKS_PER_SEC;
int main(){
	srand(time(NULL));
	int n;
	while(scanf("%d",&n)!=EOF) {
		double s1 = clock();
		for(int i = 0; i < n; i++) {
			a[i] = rand() % 10;
		}
		quick_sort(a, 0, n - 1);
	/*	for(int i = 0; i < n; i++)
		{
			printf("%d ",a[i]);
		}
		puts("");*/
		double t1 = clock();

		double s2 = clock();
		for(int i = 0; i < n; i++) {
			a[i] = rand();
		}
		merge_sort(a, 0, n - 1);
/*		for(int i = 0; i < n; i++)
		{
			printf("%d ",a[i]);
		}
		puts("");*/
		double t2 = clock();
		printf("quicksort: %.8f   mergesort: %.8f\n",(t1 - s1) / C, (t2 - s2) / C);
	}
	return 0;
}
