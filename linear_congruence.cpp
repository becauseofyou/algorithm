int extgcd(int a, int b, int &x, int &y) {
        int d = a;
        if(b != 0) {
                d = extgcd(b, a % b, y, x);
                y -= (a / b) * x;
        } else {
                x = 1;
                y = 0;
        }
        return d;
}
int mod_inverse(int a, int m) {
        int x, y;
        extgcd(a, m, x, y);
        return (m + x % m) % m;
}
int gcd(int a, int b) {
        return !b ? a : gcd(b, a % b);
}
//A[i] * x % M[i] = B[i];
std::pair<int, int> linear_congruence(const std::vector<int> &A, const std::vector<int> &B, const std::vector<int> &M) {
        int x = 0, m = 1;
        for(int i = 0; i < A.size(); i++) {
                int a = A[i] * m, b = B[i] - A[i] * x, d = gcd(M[i], a);
                if(b % d != 0) return std::make_pair(0, -1); // no solutioin
                int t = b / d * mod_inverse(a / d, M[i] / d) % (M[i] / d);
                x = x + m * t;
                m *= M[i] / d;
        }
        x = (x + m)  % m;
        return std::make_pair(x % m, m);
}
