#include <string>
#include <vector>
#include <map>
#include <list>
#include <iterator>
#include <set>
#include <queue>
#include <iostream>
#include <sstream>
#include <stack>
#include <deque>
#include <cmath>
#include <memory.h>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <algorithm>
#include <utility>
#include <ctime>
using namespace std;
typedef long long LL;
const int MAX = 67;
const int BASE = 1000000000;
struct INT
{
        int len;
        int D[MAX];
        INT() {
                len = 1;
                D[0] = 0;
        }
        INT(int x) { 
                len = 1;
                D[0] = x;
        }
        INT operator+(INT x) {
                INT res(0);
                int u = 0;
                for (int i = 0; i < max(x.len, len) || u; ++i) {
                        if (i == res.len) {
                                res.D[res.len ++] = 0;
                        }
                        int t = (i >= len ? 0 : D[i]) + (i >= x.len ? 0 : x.D[i]) + u;
                        res.D[i] = t % BASE;
                        u = t / BASE;
                }
                while (res.len > 1 && res.D[res.len-1] == 0) {
                        -- res.len;
                }
                return res;
        }
        INT operator*(INT x) {
                INT res(0);
                for (int i = 0; i < len; ++i) {
                        LL u = 0;
                        for (int j = 0; j < x.len || u; ++j) {
                                if (i+j == res.len) {
                                        res.D[res.len ++] = 0;
                                }
                                LL t = LL(D[i]) * (j >= x.len ? 0 : x.D[j]) + res.D[i+j] + u;
                                res.D[i+j] = t % BASE;
                                u = t / BASE;
                        }
                }
                while (res.len > 1 && res.D[res.len-1] == 0) {
                        -- res.len;
                }
                return res;
        }
        void print() {
                printf("%d", D[len-1]);
                for (int i = len-2; i >= 0; --i) {
                        printf("%09d", D[i]);
                }
                printf("\n");
        }
        void scan() {
                string s;
                cin >> s;
                len = 0;
                int pos = (int)s.size() - 1;
                while (pos >= 0) {
                        string t;
                        stringstream ss;
                        if (pos-8 >= 0) {
                                t = s.substr(pos-8, 9);
                        } else {
                                t = s.substr(0, pos+1);
                        }
                        ss << t;
                        ss >> D[len ++];
                        pos -= 9;
                }
        }
        bool operator<(INT& x) {
                if (len < x.len) {
                        return 1;
                }
                if (len > x.len) {
                        return 0;
                }
                for (int i = len-1; i >= 0; --i) {
                        if (D[i] < x.D[i]) {
                                return 1;
                        }
                        if (D[i] > x.D[i]) {
                                return 0;
                        }
                }
                return 0;
        }
        bool operator==(INT x) {
                if (len != x.len) {
                        return 0;
                }
                for(int i = 0; i < len; i++) {
                        if (D[i] != x.D[i]) {
                                return 0;
                        }
                }
                return 1;
        }
        INT operator/(int x) {
                INT res = *this;
                LL u = 0;
                for (int i = res.len-1; i >= 0; --i) {
                        LL t = (u * BASE + res.D[i]);
                        res.D[i] = t / x;
                        u = t % x;
                }
                while (res.len > 1 && res.D[res.len-1] == 0) {
                        -- res.len;
                }
                return res;
        }
        INT operator-(int x) {
                INT res = *this;
                LL u = 0;
                for (int i = 0; i < res.len; ++i) {
                        res.D[i] -= u * BASE + (i == 0 ? x : 0);
                        if (res.D[i] < 0) {
                                res.D[i] += BASE;
                                u = 1;
                        }
                        else {
                                u = 0;
                        }
                }
                while (res.len > 1 && res.D[res.len-1] == 0) {
                        -- res.len;
                }
                return res;
        }
};
