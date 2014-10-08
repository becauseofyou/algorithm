bool read (int &x) {
        int c = getchar (); int sign = 1;
        while (~c && c < '0' || c > '9') { if (c == '-') sign = -1; c = getchar (); }
        for (x = 0; ~c && '0' <= c && c <= '9'; c = getchar ()) x = x * 10 + c - '0';
        x *= sign;
        return ~c;
}
double get_val()
{
    double  ret(0),sign(1);
    int flag=0;
    double tmp=0.1;
    char c;
    while((c=getchar())==' '||c=='\n'||c=='\r');
    if(c=='-')
        sign=-1;
    else
        ret=c-'0';
    while((c=getchar())!=' '&&c!='\n'&&c!='\r')
    {
        if(c=='.')
        {
            flag=1;
            continue;
        }
        if(!flag)
        ret=ret*10+c-'0';
        else 
        {
            ret+=(c-'0')*tmp;
            tmp*=0.1;
        }
    }
    return ret*sign;
}
