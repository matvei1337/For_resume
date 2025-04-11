#include "header.h"

double residual1(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*fp)(double, double), int p, int k){
    int N = (nx+1)*(ny+1), i, j, i1, i2;
    thread_rows(N, p, k, i1, i2);
    double res = -1;
    for(int l = i1; l < i2; ++l){
        l_to_ij(nx, ny, i, j, l);
        if(i == nx || j == ny) continue;
        double val1 = fabs(fp(a + (i + 2./3)*hx, c + (j + 1./3)*hy) - (x[l] + x[l + 1] + x[l + 1 + nx + 1]) / 3);
        double val2 = fabs(fp(a + (i + 1./3)*hx, c + (j + 2./3)*hy) - (x[l] + x[l + nx + 1] + x[l + 1 + nx + 1]) / 3);
        res = std::max(val1, val2);
    }
    reduce_sum(p, &res, 1, &max);    
    return res;
} 

double residual2(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*fp)(double, double), int p, int k){
    int N = (nx+1)*(ny+1), i, j, i1, i2;
    thread_rows(N, p, k, i1, i2);
    double sum = 0;
    for(int l = i1; l < i2; ++l){
        l_to_ij(nx, ny, i, j, l);
        if(i == nx || j == ny) continue;
        double val1 = fabs(fp(a + (i + 2./3)*hx, c + (j + 1./3)*hy) - (x[l] + x[l + 1] + x[l + 1 + nx + 1]) / 3);
        double val2 = fabs(fp(a + (i + 1./3)*hx, c + (j + 2./3)*hy) - (x[l] + x[l + nx + 1] + x[l + 1 + nx + 1]) / 3);
        sum += val1 + val2;
    }
    sum = red_sum_det(p, k, sum);    
    return hx*hy*sum/2;
} 

double residual3(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*fp)(double, double), int p, int k){
    int N = (nx+1)*(ny+1), i, j, i1, i2;
    thread_rows(N, p, k, i1, i2);
    double res = -1;
    for(int l = i1; l < i2; ++l){
        l_to_ij(nx, ny, i, j, l);
        res = std::max(res, fabs(fp(a + i*hx, c + j*hy) - x[l]));
    }
    reduce_sum(p, &res, 1, &max);    
    return res;
}

double residual4(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*fp)(double, double), int p, int k){
    int N = (nx+1)*(ny+1), i, j, i1, i2;
    thread_rows(N, p, k, i1, i2);
    double sum = 0;
    for(int l = i1; l < i2; ++l){
        l_to_ij(nx, ny, i, j, l);
        sum += fabs(fp(a + i*hx, c + j*hy) - x[l]);
    }
    sum = red_sum_det(p, k, sum);    
    return hx*hy*sum;
}
