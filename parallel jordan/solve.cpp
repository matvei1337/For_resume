#include"functions.h"




void swap_rows(double* a,int n,int i, int j){
    int k;
    for(k = 0; k < n; k++){
        std::swap(a[i*n+k],a[j*n+k]);
    }
}

int reverse(double* C,double* D, int n,double norm){
    double leader;
    int ind_lead;
    double tmp;
    int i,j,k;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            D[i*n+j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (k = 0; k < n; k++) {
        leader = C[k*n+k];
        ind_lead = k;
        for(j = k; j < n; j++){
            if(fabs(C[j*n+k]) > fabs(leader)){
                leader = C[j*n+k];
                ind_lead = j;
            }
        }
        if (fabs(leader) < EPS*norm)  return -1;
        swap_rows(C,n,k,ind_lead);
        swap_rows(D,n,k,ind_lead);

        tmp = C[k*n+k];

        for (j = k; j < n; j++) C[k*n+j] /= tmp;
        for (j = 0; j < n; j++) D[k*n+j] /= tmp;

        for (i = 0;i < n; i++){
            if(i==k) continue;
            tmp = C[i*n+k];
            for(j = k; j < n; j++) C[i*n+j] -= tmp * C[k*n+j];
            for(j = 0; j < n; j++) D[i*n+j] -= tmp * D[k*n+j];
        }
    }

    return 0;
}

void multiply(double* A, double* B,double* C, int na, int ma,int nb, int mb){
    int v3 = na%3;
    int h3 = mb%3;
    int i = 0, j = 0, k = 0;
    for(i = 0; i < na; i++){
        for(j = 0; j < mb; j++){
            C[i*mb+j] = 0.0;
        }
    }
    double s, tmp;
    double s00, s01, s02,
           s10, s11, s12,
           s20, s21, s22;
    for(i = 0; i < v3; i++){
        for(j = 0; j < h3; j++){
            s = 0;
            for(k = 0; k < nb; k++){
                s += A[i*ma+k]*B[k*mb+j];
            }
            C[i*mb+j] = s;
        }
        for(;j < mb;j+=3){
            s00 = 0;
            s01 = 0;
            s02 = 0;
            for(k = 0; k < ma; k++){
                tmp = A[i*ma+k];
                s00 += tmp*B[k*mb+j];
                s01 += tmp*B[k*mb+j+1];
                s02 += tmp*B[k*mb+j+2];
            }
            C[i*mb+j] += s00;
            C[i*mb+j+1] += s01;
            C[i*mb+j+2] += s02;
        }
    }
    for(;i < na; i+=3){
        for(j = 0; j < h3; j++){
            s00 = 0;
            s10 = 0;
            s20 = 0;
            for(k = 0; k < ma; k++){
                tmp = B[k*mb+j];
                s00 += A[i*ma+k]*tmp;
                s10 += A[(i+1)*ma+k]*tmp;
                s20 += A[(i+2)*ma+k]*tmp;
            }
            C[i*mb+j] += s00;
            C[(i+1)*mb+j] += s10;
            C[(i+2)*mb+j] += s20;
        }
        for(;j < mb; j+=3){
            s00 = 0;
            s01 = 0;
            s02 = 0;
            s10 = 0;
            s11 = 0;
            s12 = 0;
            s20 = 0;
            s21 = 0;
            s22 = 0;
            for(k = 0; k < ma; k++){
                s00 += A[i*ma+k]*B[k*mb+j];
                s01 += A[i*ma+k]*B[k*mb+j+1];
                s02 += A[i*ma+k]*B[k*mb+j+2];
                s10 += A[(i+1)*ma+k]*B[k*mb+j];
                s11 += A[(i+1)*ma+k]*B[k*mb+j+1];
                s12 += A[(i+1)*ma+k]*B[k*mb+j+2];
                s20 += A[(i+2)*ma+k]*B[k*mb+j];
                s21 += A[(i+2)*ma+k]*B[k*mb+j+1];
                s22 += A[(i+2)*ma+k]*B[k*mb+j+2];
            }
            C[i*mb+j] += s00;
            C[i*mb+j+1] += s01;
            C[i*mb+j+2] += s02;
            C[(i+1)*mb+j] += s10;
            C[(i+1)*mb+j+1] += s11;
            C[(i+1)*mb+j+2] += s12;
            C[(i+2)*mb+j] += s20;
            C[(i+2)*mb+j+1] += s21;
            C[(i+2)*mb+j+2] += s22;
        }
    }
}


void swap_block_rows(double* a,int n,int m,int i, int j, double* C, double* D){
    int s;
    int k = n / m;

    for(s = 0; s < k; s++){
        get_block(a, C, n, m, i, s);
        get_block(a, D, n, m, j, s);
        set_block(a, D, n, m, i, s);
        set_block(a, C, n, m, j, s);
    }
}

void swap_block_columns(double* a,int n,int m,int i, int j, double* C, double* D){
    int s;
    int k = n / m;

    for(s = 0; s < k; s++){
        get_block(a, C, n, m, s, i);
        get_block(a, D, n, m, s, j);
        set_block(a, D, n, m, s, i);
        set_block(a, C, n, m, s, j);
    }
}



int solve(double* a, double* b, double* x, int n, int m, int p, int k,double* b1,double* b2,double* b3,double norm,double* copy_row,double* copy_b){
    int d = n / m;
    int l = n - d * m;
    int h = l ? d + 1 : d;
    int u,v;
    int I,J;

    for (int t = 0; t < d; ++t) {
        if(t%p != k){
            copy_matrix(copy_row, a+t*m*n, m, n);
            copy_matrix(copy_b, b+t*m, m, 1);
        }
        reduce_sum<int>(p);

        if(t%p == k) get_block(a,b1,n,m,t,t);
        else get_block(copy_row,b1,n,m,0,t);
        int res = reverse(b1,b2,m,norm);
        if (res == -1) {
            return -1;
        }

        reduce_sum<int>(p);

        for(u = t+1; u < d; u++){
            if(t%p == k){
                get_block(a,b1,n,m,t,u);
                multiply(b2,b1,b3,m,m,m,m);
                set_block(a,b3,n,m,t,u);
            }
            else{
                get_block(copy_row,b1,n,m,0,u);
                multiply(b2,b1,b3,m,m,m,m);
                set_block(copy_row,b3,n,m,0,u);
            }
        }

        if(l > 0){
            if(t%p == k){
                get_block(a,b1,n,m,t,d);
                multiply(b2,b1,b3,m,m,m,l);
                set_block(a,b3,n,m,t,d);
            }
            else{
                get_block(copy_row,b1,n,m,0,d);
                multiply(b2,b1,b3,m,m,m,l);
                set_block(copy_row,b3,n,m,0,d);
            }
        }

        if(t%p == k){
            multiply(b2,b+m*t,b3,m,m,m,1);
            for(u = 0; u < m; u++){
                b[t*m+u] = b3[u];
            }
        }
        else{
            multiply(b2,copy_b,b3,m,m,m,1);
            for(u = 0; u < m; u++){
                copy_b[u] = b3[u];
            }
        }
        reduce_sum<int>(p);

        for(int q = k; q < h; q+=p){
            if(q == t) continue;
            get_block(a,b1,n,m,q,t);

            int multiplier_rows = q < d ? m : l;
            for(v = t+1; v < h; v++){
                if(t%p == k) get_block(a,b2,n,m,t,v);
                else  get_block(copy_row,b2,n,m,0,v);
                int block_cols = v < d ? m : l;
                multiply(b1,b2,b3,multiplier_rows,m,m,block_cols);
                get_block(a,b2,n,m,q,v);
                for (I = 0; I < multiplier_rows; ++I) {
                    for (J = 0; J < block_cols; ++J) {
                        b2[block_cols*I + J] -= b3[block_cols*I + J];           
                    }
                }
                set_block(a,b2,n,m,q,v);
            }

            if(t%p == k) multiply(b1,b+m*t,b3,multiplier_rows,m,m,1);
            else multiply(b1,copy_b,b3,multiplier_rows,m,m,1);
            for(u = 0; u < multiplier_rows; u++){
                b[q*m+u] -= b3[u];
            }
        }
        reduce_sum<int>(p);
    }

    if(l != 0){
        int res = 0;
        if(d%p == k){
            get_block(a,b1,n,m,d,d);
            res = reverse(b1, b2, l,norm); 
            if (res == 0) {
                multiply(b2,b+m*d,b3,l,l,l,1);
                for(u = 0; u < l; u++){
                    b[d*m+u] = b3[u];
                }
            }
        }

        reduce_sum(p, &res, 1);
        if (res < 0) { 
            return -1;
        }
        if(d%p != k){
            copy_matrix(copy_b, b+d*m, m, 1);
        }
        reduce_sum<int>(p);
        for (int j = k; j < d; j+=p) { 
            get_block(a,b1,n,m,j,d);
            if(d%p == k) multiply(b1,b+m*d,b3,m,l,l,1);
            else multiply(b1,copy_b,b3,m,l,l,1);
            for(u = 0; u < m; u++){
                b[j*m+u] -= b3[u];
            }
        }
    }
    reduce_sum<int>(p);

    for (int i = k; i < h; i += p) {
        int to = (i < d) ? m : l;
        for (int j = 0; j < to; ++j) {
            x[i*m + j] = b[i*m + j];       
        }   
    }
    return 0;

}



void* thread_func(void* ptr){
    Args* arg = (Args*)ptr;
    double* a = arg->a;
    double* b = arg->b;
    double* x = arg->x;

    int n = arg->n;
    int p = arg->p;
    int k = arg->k;
    int m = arg->m;
    int s = arg->s;
    int r = arg->r;
    const char* name = arg->filename;

    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus-1-(k%n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();
    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

    for (int i = k*m; i < n; i += p*m) {
        int h = (i + m < n) ? m : n - i;
        memset(a+i*h,0,h*n*sizeof(double));
        memset(b + i, 0, h*sizeof(double));
        memset(x + i, 0, h*sizeof(double));
    }

    reduce_sum<int>(p);
    double* block1 = new double[m*m];
    double* block2 = new double[m*m];
    double* block3 = new double[m*m];
    double* copy_row = new double[m*n];
    double* copy_b = new double[m];

    if(name != nullptr){
        int res = 0;
        if(k == 0) res = read(a,n,name);
        reduce_sum(p,&res,1);
        if(res < 0){
            arg->status = res;
            delete[] block1;
            delete[] block2;
            delete[] block3;
            delete[] copy_row;
            delete[] copy_b;
            return nullptr;
        }
    }
    else init_a(a,n,m,s,k,p);
    init_b(a,b,n,m,k,p);

    double norm = 0;
    if(k == 0){
        print(a,n,n,r);

        norm = matrix_norma(a,n);
    }
    reduce_sum(p,&norm,1);

    double t = get_cpu_time();
    int res = solve(a,b,x,n,m,p,k,block1,block2,block3,norm,copy_row,copy_b);
    t = get_cpu_time() - t;
    reduce_sum(p,&res,1);
    //unsigned int g = t*1.2/2;
    //sleep(g);
    arg->status = res;
    arg->t_cpu = t;
    if (res < 0) {
        delete[] block1;
        delete[] block2;
        delete[] block3;
        delete[] copy_row;
        delete[] copy_b;
        return nullptr;
    }
    if(k == 0){
        print(x, 1, n, r);
    }
    delete[] block1;
    delete[] block2;
    delete[] block3;
    delete[] copy_row;
    delete[] copy_b;
    reduce_sum<int>(p);
    return nullptr;

}

