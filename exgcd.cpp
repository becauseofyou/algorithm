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
