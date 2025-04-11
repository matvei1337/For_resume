#include <stdio.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <stddef.h>
#include <cstring>
#include <math.h>


class Args{
public:
    double a;
    double b;
    double c;
    double d;
    double eps;
    int nx;
    int ny;
    int k;
    int maxit;
    int p;
    int* I;
    double* A;
    double* B;
    double* x;
    double* r;
    double* u;
    double* v;
    double (*f)(double, double);
    int its = 0;
    double r1 = 0;
    double r2 = 0;
    double r3 = 0;
    double r4 = 0;
    double t1 = 0;
    double t2 = 0;
};

int init_red_sum(int p);
double red_sum_det(int p, int k, double s);
void free_res();

void matrix_multiply_vector(int n, double* A, int* I, double* x, double* y, int p, int k);
void preconditioner_for_msr(int n, double* A, int* I, double* v1, double* v2, int flag, int p, int k);
int step(int n, double* A, int* I, double* x, double* r, double* u, double* v, double prec, int p, int k);
int min_res_msr(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, double eps, int maxit, int p, int k);
int min_res_msr_full(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, double eps, int maxit, int maxsteps, int p, int k);

void thread_rows(int n, int p, int k, int& i1, int& i2);
double scalar_multiply(int n, double* x, double* y, int p, int k);
void multply_sub_vector(int n, double* x, double* y, double tau, int p, int k);
void ij_to_l(int nx, int, int i, int j, int& l);
void l_to_ij(int nx, int, int& i, int& j, int l);
int get_diag(int nx, int ny, int i, int j, int* I_ij = nullptr);
int get_len_diag(int nx, int ny);
int alloc_msr(int nx, int ny, double** p_A, int** p_I);
void fill_I(int nx, int ny, int* I);
int check_symm(int nx, int ny, int* I, double* A, double eps, int p, int k);
void fill_A_ij(int nx, int ny, double hx, double hy, int i, int j, double* A_diag, double* A_in_diag);
double F_IJ(int nx, int ny, double hx, double hy, double a, double с, int i, int j, double (*f)(double, double));

void fill_A(int nx, int ny, double hx, double hy, int* I, double* A, int p, int k);
void fill_B(int nx, int ny, double hx, double hy, double a, double c, double* B, double (*f)(double, double), int p, int k);
void solve_r(int n, int* I, double* U, double* b, double* x, double w, int p, int k);
void solve_l(int n, int* I, double* U, double* b, double* x, double w, int p, int k);

double residual1(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*fp)(double, double), int p, int k);
double residual2(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*fp)(double, double), int p, int k);
double residual3(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*fp)(double, double), int p, int k);
double residual4(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*fp)(double, double), int p, int k);

void* solution(void* ptr);


template<class T>
void sum(T* r, T* a, int n){
    for(int i = 0; i < n; ++i){
        r[i] += a[i];
    }    
}

template<class T>
void max(T* r, T* a, int n){
    for(int i = 0; i < n; ++i){
        r[i] = std::max(r[i], a[i]);
    }  
}

template<class T>
void reduce_sum(int p, T* a = nullptr, int n = 0, void (*func)(T*, T*, int) = &sum){
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static T* r = nullptr;
    int i;

    if(p <= 1) return;

    pthread_mutex_lock(&m);

    if(r == nullptr) r = a;
    else func(r, a, n);

    t_in++;
    if(t_in >= p){
        t_out = 0;
        pthread_cond_broadcast(&c_in);
        
    }
    else{
        while(t_in < p){
            pthread_cond_wait(&c_in, &m);
        }
    }

    if(r != a){
        for(i = 0; i < n; ++i){
            a[i] = r[i];
        }
    } 

    t_out++;
    if(t_out >= p){
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);        
    }
    else{
        while(t_out < p){
            pthread_cond_wait(&c_out, &m);
        }
    }

    pthread_mutex_unlock(&m);
}
