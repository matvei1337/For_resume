#include "header.h"
#define MAXSTEP 300

double get_cpu_time(){
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec*1e-6;
}

void* solution(void* ptr) {
    Args* args = (Args*)ptr;
    double a = args->a, b = args->b, c = args->c, d = args->d,  eps = args->eps;
    int nx = args->nx, ny = args->ny, k = args->k, maxit = args->maxit, p = args->p;
    int* I = args->I;
    double* A = args->A;
    double* B = args->B;
    double* x = args->x;
    double* r = args->r;
    double* u = args->u;
    double* v = args->v;
    double (*f)(double, double) = args->f;

    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus - 1 - (k % n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();
    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

    double hx = (b - a)/nx;
    double hy = (d - c)/ny;
    int N = (nx + 1)*(ny + 1);
    
    fill_A(nx, ny, hx, hy, I, A, p, k);
    fill_B(nx, ny, hx, hy, a, c, B, f, p, k);
    

    args->t1 = get_cpu_time();
    args->its = min_res_msr_full(N, A, I, B, x, r, u, v, eps, maxit, MAXSTEP, p, k);
    args->t1 = get_cpu_time() - args->t1;


    args->t2 = get_cpu_time();

    args->r1 = residual1(nx, ny, a, c, hx, hy, x, f, p, k);
    args->r2 = residual2(nx, ny, a, c, hx, hy, x, f, p, k);
    args->r3 = residual3(nx, ny, a, c, hx, hy, x, f, p, k);
    args->r4 = residual4(nx, ny, a, c, hx, hy, x, f, p, k);

    args->t2 = get_cpu_time() - args->t2;

    reduce_sum<int>(p);
    return nullptr;
}
