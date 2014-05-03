struct Exp
{
	bool error;
	int tot,top;
	int num[N];
	char op[N];
	void ini()
	{
		error=false;
		tot=0;
		top=1;
		op[1]='(';
	}
	bool prior(char a,char b)
	{
		if(b=='+'||b=='-')
			return a!='(';
		return a=='*'||a=='/';
	}
	int cal(char c,int a,int b)
	{
		if(c=='+') return a+b;
		if(c=='-') return a-b;
		if(c=='*') return a*b;
		if(b!=0) return a/b;
		error=true;
		return 0;
	}
	bool digit(char ch)
	{
		return ch>='0'&&ch<='9';
	}
	int solve(char *s,int len)
	{
		s[len++]=')';
		for(int i=0;i<len;i++)
		{
			if(s[i]=='(') op[++top]=s[i];
			else if(s[i]==')')
			{
				while(top>0&&op[top]!='(')
				{
					num[tot-1]=cal(op[top],num[tot-1],num[tot]);
					tot--;
					top--;
				}
				top--;
			}
			else if(s[i]=='-'&&(i==0||s[i-1]=='('))
			{
				num[++tot]=0;
				op[++top]='-';
			}
			else if(digit(s[i]))
			{
				int t=s[i]-'0';
				for(i++;digit(s[i]);i++)
				{
					t=t*10+s[i]-'0';
				}
				num[++tot]=t;
				i--;
			}
			else
			{
				while(top>0&&prior(op[top],s[i]))
				{
					num[tot-1]=cal(op[top],num[tot-1],num[tot]);
					tot--;
					top--;
				}
				op[++top]=s[i];
			}
		}
		return num[1];
	}
}E;
{
   E.ini();
   int res=E.solve(s,len);//s为表达式
   return (res==tmp && !E.error);
}
