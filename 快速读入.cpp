inline int read()    
{    
    char ch;  
    bool flag = false;  
    int a = 0;    
    while(true) {
            ch = getchar();
            if(ch >= '0' && ch <= '9' || ch == '-') {
                    break;
            }
    }
    (ch == '-') ? (flag = true) : (a = a * 10 + ch - '0');
    while(true) {  
            ch = getchar(); 
            if(ch >= '0' && ch <= '9') {
                    a = a * 10 + ch - '0';
            } else {
                    break;
            }
    }     
    if(flag) {  
            a = -a;  
    }  
    return a;    
}    
