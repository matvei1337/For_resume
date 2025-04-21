#include"functions.h"
#define MAX_ITER 20000

void three_diag(double *a, int n, double *x,double eps){
    double scalar = 0;
    double h = 0;
    int i,j,k;

    for(i = 0; i < n-2; i++){
        scalar = 0;
        for(j = i+2; j < n; j++){
            scalar += a[j*n+i]*a[j*n+i];
            x[j] = a[j*n+i];
        }

        h = sqrt(a[(i+1)*n+i]*a[(i+1)*n+i] + scalar);
        if(h < 1e-50*eps){
            a[(i+1)*n+i] = 0;
            a[(i+2)*n+i] = 0;
            continue;
        }
        if(scalar < 1e-50*eps){
            a[(i+2)*n+i] = 0;
            continue;
        }

        x[i+1] = a[(i+1)*n+i] - h;
        scalar = sqrt(x[i+1]*x[i+1] + scalar);

        if(scalar > 1e-50*eps){
            x[i+1] /= scalar;
            for(j = i+2; j < n; j++) x[j] /= scalar;

            for(j = i+1; j < n; j++){
                scalar = 0;
                for(k = i + 1; k < n; k++) scalar += x[k] * a[k*n+j];
                for(k = i+1; k < n; k++) a[k*n+j] -= 2*scalar*x[k];
            }

            for(j = 0; j < n; j++) {
                scalar = 0;
                for(k = i+1; k < n; k++)  scalar += x[k]*a[j*n+k];
                for(k = i+1; k < n; k++)  a[j*n+k] -= 2*scalar*x[k];
            }
            a[(i+1)*n+i] = h;
        }
    }
}




int solve(double* a, int n, double* x, double eps, double *x1, double *x2, int& its){
    double sdvig = 0;
    double h;
    double scalar;
    double b,c,D;
    int iter = 0;
    int i = 0,j = 0,k = 0;
    for(k = n; k > 2; k--) {

        while(fabs(a[(k-1)*n+k-2]) > eps) {
            sdvig = a[(k-1)*n+k-1] + a[(k-1)*n+k-2]/2;
            for(i = 0; i < k; i++) {
                a[i*n+i] -= sdvig;
            }


            for(i = 0; i < k-1; i++){
                scalar = a[(i+1)*n+i] * a[(i+1)*n+i];
                if(scalar < 1e-35*eps){
                    x1[i] = 0;
                    x2[i] = 0;
                    continue;
                }
                h = sqrt(a[i*n+i]*a[i*n+i] + scalar) * ((a[i*n+i] > 0) ? -1 : 1);
                a[i*n+i] -= h;
                scalar = sqrt(a[i*n+i]*a[i*n+i] + scalar);
                x1[i] = a[i*n+i] / scalar;
                x2[i] = a[(i+1)*n+i] / scalar;

                for(j = i+1; j < k; j++){
                    scalar = a[i*n+j]*x1[i] + a[(i+1)*n+j]*x2[i];
                    a[i*n+j] -= 2*scalar*x1[i];
                    a[(i+1)*n+j] -= 2*scalar*x2[i];
                }
                a[i*n+i] = h;
                a[(i+1)*n+i] = 0;
            }


            for(i = 0; i < k-1; i++){
                if (fabs(x1[i]) < 1e-35*eps && fabs(x2[i]) < 1e-35*eps) continue;

                for(j = 0; j < i+2; j++){
                    scalar = a[j*n+i]*x1[i] + a[j*n+i+1]*x2[i];
                    a[j*n+i] -= 2*scalar*x1[i];
                    a[j*n+i+1] -= 2*scalar*x2[i];
                }
            }

            for(i = 0; i < k; i++) {
                a[i*n+i] += sdvig;
            }

            iter++;
            if(iter > MAX_ITER){
                its = iter;
                return -1;
            }
        }
    }
    if(n > 1){
        b = -(a[0*n+0]+a[1*n+1]);
        c = a[0*n+0]*a[1*n+1] - a[0*n+1]*a[1*n+0];
        D = b*b-4*c;
        if(D < 0){
            its = iter;
            return -1;
        }
        x[1] = (b>0) ? (-b-sqrt(D))/2 : (-b+sqrt(D))/2;
        x[0] = c/x[1];
    }
    for(i = 2; i < n; i++) x[i] = a[i*n+i];
    its = iter;
    return 0;
}
