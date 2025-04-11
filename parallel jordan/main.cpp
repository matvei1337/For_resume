#include"functions.h"
using namespace std;


int main(int argc, char *argv[]){
    double* a;
    double* b;
    double* x;
    int task = 17;
    int n,m,p,r,s;
    int k;
    double r1, r2, t = 0;
    const char* filename = nullptr;
    //разобрали аргументы командной строки
    if( !(argc == 6 || argc == 7) || sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1 || sscanf(argv[3], "%d", &p) != 1
           || sscanf(argv[4], "%d", &r) != 1 || sscanf(argv[5], "%d", &s) != 1){
        printf("Usage: %s n m p r s (file)\n", argv[0]);
        return -1;
    }
    if(argc == 7){
        filename = argv[6];
        if(s!=0) return -1;
    }
    //выделяем память но не обращаемся к ней
    a = new double[n*n];
    b = new double[n];
    x = new double[n];
    Args* arg = new Args[p];
    pthread_t* tid = new pthread_t[p];

    for(k = 0; k < p; k++){
        arg[k].a = a;
        arg[k].b = b;
        arg[k].x = x;
        arg[k].k = k;
        arg[k].m = m;
        arg[k].n = n;
        arg[k].p = p;
        arg[k].s = s;
        arg[k].r = r;
        arg[k].filename = filename;
    }

    double elapsed = get_full_time();
    for(k = 1;k < p;k++){
        if(pthread_create(tid + k,0,thread_func,arg+k)){
            printf("Can not create thread number: %d\n",k);
            abort();
        }
    }
    thread_func(arg+0);
    for(k = 1; k < p; k++){
        pthread_join(tid[k], 0);
    }
    elapsed = get_full_time() - elapsed;

    for (k = 0; k < p; k++) {
        if(arg[k].status < 0) {
            printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
            argv[0], task, -1., -1., 0., 0., s, n, m, p);
            delete[] a;
            delete[] b;
            delete[] x;
            delete[] arg;
            delete[] tid;
            return -1;
        }
    }
    for(k = 0; k < p; k++){
        printf("Thread %d = %.2f\n",k,arg[k].t_cpu);
    }

    if(argc == 6) init_a_prosto(a, n, s);
    else read(a, n, filename);
    init_b_prosto(b, a, n);

    t = clock();
    r1 = norma1(a, x, b, n);
    r2 = norma2(x, n);
    t = (clock()-t)/CLOCKS_PER_SEC;

    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
     argv[0], task, r1, r2, elapsed, t, s, n, m, p);
    delete[] a;
    delete[] b;
    delete[] x;
    delete[] arg;
    delete[] tid;
    return 0;

}
