bool read (int &x) {
        int c = getchar (); int sign = 1;
        while (~c && c < '0' || c > '9') { if (c == '-') sign = -1; c = getchar (); }
        for (x = 0; ~c && '0' <= c && c <= '9'; c = getchar ()) x = x * 10 + c - '0';
        x *= sign;
        return ~c;
}
