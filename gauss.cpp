#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <cstring>
const double eps=  1e-8;
double A[20][20];
double x[20];
int equation;
int variable;
void gauss()
{
    int col=0, row=0;  
    for(; row<equation && col<variable; row++, col++)  
    {  
        int max_row=row;  
        for(int i=row+1; i<equation; i++)  
            if(fabs(A[i][col])>fabs(A[max_row][col]))  
                max_row=i;  
        //如果 col 那列元素最大是 0，表明这一列全部是 0，处理下一列  
        if(fabs(A[max_row][col])<eps)  
        {  
            row--;  
            continue;  
        }  
        //如果不是同一行，交换元素  
        if(max_row!=row)  
        {  
            for(int i=col; i<=variable; i++)  
                std::swap(A[max_row][i], A[row][i]);  
        }  
        //枚举要删去的行  
        for(int i=0; i<equation; i++)   if(i!=row)
        {  
            if(fabs(A[i][col])==0) continue;  
            double ta=A[row][col];  
            double tb=A[i][col];  
            for(int j=col; j<=variable; j++)  
                A[i][j]=A[i][j]-A[row][j]*tb/ta;  
        }  
    }
    for(int i = row;i < equation; i++) if(A[i][variable]!=0) return -1;
    for(int i = 0; i < row; i++) 
    {
        x[i] = A[i][variable] / A[i][i];
    }
}
