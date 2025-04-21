#include"functions.h"

double f(int k, int n, int i, int j){
    if(k == 1) return n - std::max(i,j);
    else if(k == 2){
        if(i == j) return 2.0;
        else if(abs(i - j) == 1) return -1.0;
        else return 0.0;
    }
    else if(k == 3){
        if(i == j && j < n-1 && i < n-1) return 1.0;
        else if(j == n-1) return i+1;
        else if(i == n-1) return j+1;
        else return 0.0;
    }
    else return 1.0/(1.0*(i+j+1));
}

int read(double* a, int n, const char *name){
    int i = 0, j = 0;

    FILE* file;
    file = fopen(name, "r");
    if (!file){
        std::cerr << "Cann't read file";
        return -1;
    }

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if(fscanf(file, "%lf", a+i*n+j) != 1){
                fclose(file);
                return -1;
            }
        }
    } 

    fclose(file);
    return 1;
}

void print(double* a, int n, int m, int r){
    int N = (n < r) ? n : r;
    int M = (m < r) ? m : r;
    int i = 0, j = 0;
    for(i = 0; i < N; i++){
        for(j = 0; j < M; j++){
            printf(" %10.3e",a[i*n+j]);
        }
        printf("\n");
    }
}

void init_a(double* a, int n, int k){
    int i = 0, j = 0;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            a[i*n+j] = f(k,n,i,j);
        }
    }
}



double norma1(double* a, double* x, int n){
    int i = 0;
    double norm1 = 0.0, norm2 = 0.0;
    for(i = 0; i < n; i++){
        norm1 += a[i*n+i];
        norm2 += x[i];
    }
    return fabs(norm1 - norm2);
}

double norma2(double* a, double* x, int n){
    int i = 0, j = 0;
    double norm1 = 0.0, norm2 = 0.0;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            norm1 += a[i*n+j]*a[j*n+i];
        }
        norm2 += x[i]*x[i];
    }
    norm1 = sqrt(norm1);
    norm2 = sqrt(norm2);
    return fabs(norm1 - norm2);
}

double matrix_norma(double* a, int n){
    int i = 0, j = 0;
    double norma = 0.0, s;
    for(j = 0; j < n; j++){
        s = 0.0;
        for(i = 0; i < n; i++){
            s += fabs(a[i*n+j]);
        }
        norma = std::max(s, norma);
    }
    return norma;
}

