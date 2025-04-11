#include "header.h"

double f0(double, double){
    return 1;
}

double f1(double x, double){
    return x;
}

double f2(double, double y){
    return y;
}

double f3(double x, double y){
    return x + y;
}

double f4(double x, double y){
    return sqrt(x*x + y*y);
}

double f5(double x, double y){
    return x*x + y*y;
}

double f6(double x, double y){
    return exp(x*x - y*y);
}

double f7(double x, double y){
    return 1.0/(25*(x*x + y*y) + 1);
}


int main(int argc, char* argv[]) {
    double a, b, c, d, eps;
    int nx, ny, k, maxit, p;
    double (*fp)(double, double);
    double r1, r2, r3, r4, t1, t2;
    if (argc != 11) {
        printf("Usage: %s a b c d nx ny k eps maxit p", argv[0]);
        return -1;
    }
    if(sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1
        || sscanf(argv[3], "%lf", &c) != 1 || sscanf(argv[4], "%lf", &d) != 1
        || sscanf(argv[5], "%d", &nx) != 1 || sscanf(argv[6], "%d", &ny) != 1
        || sscanf(argv[7], "%d", &k) != 1 || sscanf(argv[8], "%lf", &eps) != 1
        || sscanf(argv[9], "%d", &maxit) != 1 || sscanf(argv[10], "%d", &p) != 1){
        printf("Wrong input!\n");
        return 0;
    }
    
    int* I = nullptr;
    double* A = nullptr;
    if(alloc_msr(nx, ny, &A, &I)) return -2;
    init_red_sum(p);
    int N = (nx + 1) * (ny + 1);
    double* B = new double[N];
    double* x = new double[N];
    double* r = new double[N];
    double* u = new double[N];
    double* v = new double[N];

    fill_I(nx, ny, I);

    memset(x, 0, N*sizeof(double));

    if(k == 0) fp = f0;
    else if(k == 1) fp = f1;
    else if(k == 2) fp = f2;
    else if(k == 3) fp = f3;
    else if(k == 4) fp = f4;
    else if(k == 5) fp = f5;
    else if(k == 6) fp = f6;
    else fp = f7;

    Args* args = new Args[p];
    pthread_t* t_id = new pthread_t[p];
        
    for(int i = 1; i < p; ++i){
        args[i].a = a;
        args[i].b = b;
        args[i].c = c;
        args[i].d = d;
        args[i].eps = eps;
        args[i].nx = nx;
        args[i].ny = ny;
        args[i].maxit = maxit;
        args[i].p = p;
        args[i].k = i;
        args[i].I = I;
        args[i].A = A;
        args[i].B = B; 
        args[i].x = x;
        args[i].r = r;
        args[i].u = u;
        args[i].v = v;
        args[i].f = fp;
        pthread_create(t_id + i, 0, &solution, args + i); 
    }

    args[0].a = a;
    args[0].b = b;
    args[0].c = c;
    args[0].d = d;
    args[0].eps = eps;
    args[0].nx = nx;
    args[0].ny = ny;
    args[0].maxit = maxit; 
    args[0].p = p;
    args[0].k = 0;
    args[0].I = I;
    args[0].A = A;
    args[0].B = B;
    args[0].x = x;
    args[0].r = r;
    args[0].u = u;
    args[0].v = v;
    args[0].f = fp;
    solution(args + 0);

    for(int i = 1; i < p; ++i){
        pthread_join(t_id[i], 0);
    }    

    int its = args[0].its;
    r1 = args[0].r1;
    r2 = args[0].r2;
    r3 = args[0].r3;
    r4 = args[0].r4;
    t1 = args[0].t1;
    t2 = args[0].t2;

    printf(
    "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\n\
          It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
          argv[0], 3, r1, r2, r3, r4, t1, t2, its, eps, k, nx, ny, p);

    free_res();
    delete[] I;
    delete[] A;
    delete[] B;
    delete[] x;
    delete[] r;
    delete[] u;
    delete[] v;
    delete[] args;
    delete[] t_id;
    return 0;
}
