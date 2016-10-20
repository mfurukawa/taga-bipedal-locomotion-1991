#pragma once
/* (c) Masahiro Furukawa
*
* m.furukawa@ist.osaka-u.ac.jp
*
* Created Oct 11, 2016
*
* Description
* - gauss-jordan method to obtain inverce matrix and solution.
*
* referece:
* http://www.yamamo10.jp/yamamoto/lecture/2006/5E/Linear_eauations/gaussj.pdf
*/

#define GAUSS_JORDAN_MAXN 8

int gauss_jordan(int n, long double a[][GAUSS_JORDAN_MAXN + 10], long double b[]);
long double lfabs(long double x);
