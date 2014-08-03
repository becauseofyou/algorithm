     {
      int mul = 1;
        for(int i = 0; i < 3; i++) mul = mul * p[i];
        int ret = 0;
        for(int i = 0; i < 3; i++) {
            int x, y;
            int g = exgcd(mul / p[i], p[i], x, y);
            ret += mul / p[i] *  x   * ans[i];
            ret %= mul;
        }
        printf("%d\n", (ret % mul + mul) % mul);
    }
    int exgcd(int a, int b, int &x, int &y) {
        int gcd = a;
        if(b == 0) {
            x = 1; y = 0;
        } else {
            gcd = exgcd(b, a % b, x, y);
            x -= a / b * y;
            std::swap(x, y);
        }
        return gcd;
    }
