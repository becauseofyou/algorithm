typedef vector<int> vec;
typedef vector<vec> mat;

mat mul(const mat& a, const mat& b) {
    mat res(a.size(), vec(b[0].size(), 0));
    for (int r = 0; r < a.size(); r++)
        for (int c = 0; c < b[0].size(); c++)
            for (int k = 0; k < b.size(); k++)
                res[r][c] = (res[r][c] + a[r][k] * b[k][c] % M) % M;

    return res;
}

mat pow(const mat& a, int b) {
    mat ans(a.size(), vec(a.size(), 0));
    for (int i = 0; i < a.size(); i++) ans[i][i] = 1;

    mat base = a;
    while (b) {
        if (b & 1) ans = mul(ans, base);
        base = mul(base, base);
        b >>= 1;
    }

    return ans;
}
