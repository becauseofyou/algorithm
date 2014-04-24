#define zero(x) ((x>0? x:-x)<1e-15)
const int N = 100;
double a[N][N];
double b[N][N];
int g[110][110];
double det(double a[][N],int n) {
    int i, j, k, sign = 0;
    double ret = 1, t;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            b[i][j] = a[i][j];
    for (i = 0; i < n; i++) {
        if (zero(b[i][i])) {
            for (j = i + 1; j < n; j++)
                if (!zero(b[j][i]))
                    break;
            if (j == n)
                return 0;
            for (k = i; k < n; k++)
                swap(b[i][k], b[j][k]);
            sign++;
        }
        ret *= b[i][i];
        for (k = i + 1; k < n; k++)
            b[i][k] /= b[i][i];
        for (j = i + 1; j < n; j++)
            for (k = i + 1; k < n; k++)
                b[j][k] -= b[j][i] * b[i][k];
    }
    if (sign & 1)
        ret = -ret;
    return ret;
}
int main() {
    /* 0 based
     * a[i][i] = degree[i];
     * if(g[i][j]) {
     *     a[i][j]--;
     * }
     * ans = det(a, n - 1);
     */
}

