#include "header.h"

#define FUNC(I, J) do {ij_to_l(nx, ny, I, J, k); if (I_ij) { I_ij[m] = k; } m++;} while (0)
#define F(I, J) (f(a + (I)*hx, с + (J)*hy))

static double* results = nullptr;

int init_red_sum(int p){
    results = new double[p];
    if(results == nullptr) return -1;
    return 0;
}

double red_sum_det(int p, int k, double s){
    double sum = 0; int l;
    results[k] = s;
    reduce_sum<int>(p);
    for(l = 0; l < p; ++l){
        sum += results[l];
    }
    reduce_sum<int>(p);
    return sum;
}

void free_res(){
    delete[] results;
}


void matrix_multiply_vector(int n, double* A, int* I, double* x, double* y, int p, int k){
    int i, i1, i2, l, J; double s;
    thread_rows(n, p, k, i1, i2);
    for(i = i1; i < i2; ++i){
        s = A[i] * x[i];
        l = I[i+1] - I[i];
        J = I[i];
        for(int j = 0; j < l; ++j){
            s += A[J + j] * x[I[J + j]];
        }
        y[i] = s;
    }
}

void preconditioner_for_msr(int n, double* A, int* I, double* v1, double* v2, int flag, int p, int k){
    double w = 1;
    if(flag) solve_l(n, I, A, v2, v1, w, p, k);
    else solve_r(n, I, A, v2, v1, w, p, k);
    reduce_sum<int>(p);
}

int step(int n, double* A, int* I, double* x, double* r, double* u, double* v, double prec, int p, int k){
    matrix_multiply_vector(n, A, I, v, u, p, k);
    double c1 = scalar_multiply(n, u, r, p, k);
    double c2 = scalar_multiply(n, u, u, p, k);
    if(c1 < prec || c2 < prec) return 1;
    double tau = c1 / c2;
    multply_sub_vector(n, x, v, tau, p, k);
    multply_sub_vector(n, r, u, tau, p, k);
    return -1;
}

int min_res_msr(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, double eps, int maxit, int p, int k){
    double prec, b_norm2;
    int it;
    b_norm2 = scalar_multiply(n, b, b, p, k);
    prec = b_norm2 * eps * eps;
    matrix_multiply_vector(n, A, I, x, r, p, k);
    multply_sub_vector(n, r, b, 1., p, k);
    
    for(it = 0; it < maxit; ++it){
        preconditioner_for_msr(n, A, I, v, r, 0, p, k);
        if(step(n, A, I, x, r, u, v, prec, p, k) == 1) break;

        matrix_multiply_vector(n, A, I, x, u, p, k);
        multply_sub_vector(n, u, b, 1., p, k);
        preconditioner_for_msr(n, A, I, v, u, 1, p, k);
        if(step(n, A, I, x, r, u, v, prec, p, k) == 1) break;
    }

    if(it >= maxit) return -1;

    return it;
}

int min_res_msr_full(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, double eps, int maxit, int maxstep, int p, int k){
    int step, ret, its = 0;
    for(step = 0; step < maxstep; ++step){
        ret = min_res_msr(n, A, I, b, x, r, u, v, eps, maxit, p, k);
        if(ret >= 0){
            its += ret;
            break;
        }
        its += maxit;
    }

    if(step >= maxstep) return -1;

    return its;
}

void thread_rows(int n, int p, int k, int& i1, int& i2){
    i1 = n*k;
    i1 /= p;
    i2 = n*(k + 1);
    i2 /= p; 
}

double scalar_multiply(int n, double* x, double* y, int p, int k){
    int i1, i2, i;
    double s = 0;
    thread_rows(n, p, k, i1, i2);
    for(i = i1; i < i2; ++i){
        s += x[i] * y[i];
    }
    s = red_sum_det(p, k, s);
    return s;
}

void multply_sub_vector(int n, double* x, double* y, double tau, int p, int k){
    int i, i1, i2;
    thread_rows(n, p, k, i1, i2);
    for(i = i1; i < i2; ++i){
        x[i] -= tau * y[i];
    }
    reduce_sum<int>(p);
}

void ij_to_l(int nx, int, int i, int j, int& l){
    l = i + j * (int)(nx + 1);
}

void l_to_ij(int nx, int, int& i, int& j, int l){
    j = l / (nx + 1);
    i = l - j * (nx + 1);
}


int get_diag(int nx, int ny, int i, int j, int* I_ij){
    int m = 0; int k;
    if(i < nx) FUNC(i+1, j);
    if(j > 0) FUNC(i, j-1);
    if(i > 0 && j > 0) FUNC(i-1, j-1);
    if(i > 0) FUNC(i - 1, j);
    if(j < ny) FUNC(i, j+1);
    if(i < nx && j < ny) FUNC(i+1, j+1);
    return m;
}

int get_len_diag(int nx, int ny){
    int m = 0;
    for(int i = 0; i <= nx; ++i){
        for(int j = 0; j <= ny; ++j){
            m += get_diag(nx, ny, i, j);
        }
    }    
    return m;
}

int alloc_msr(int nx, int ny, double** p_A, int** p_I){
    int diag_len = (nx+1)*(ny+1);
    int off_diag = get_len_diag(nx, ny);
    int len = diag_len + off_diag + 1;

    double* A = nullptr;
    int* I = nullptr;
    A = new double[len];
    if(A == nullptr) return 1;
    I = new int[len];
    if(I == nullptr) return 2;
    *p_A = A;
    *p_I = I;
    return 0;  
}

void fill_I(int nx, int ny, int* I){ 
    int N = (nx+1) * (ny+1);
    int i, j, l, m, r = N + 1;
    for(l = 0; l < N; ++l){
        l_to_ij(nx, ny, i, j, l);
        I[l] = r;
        m = get_diag(nx, ny, i, j, I + r);
        r += m;
    }
    I[l] = r;
}

int check_symm(int nx, int ny, int* I, double* A, double eps, int p, int k){
    int l1, l2, l;
    int N = (nx+1)*(ny+1);
    l1 = (k * N) / p;
    l2 = ((k + 1) * N) / p;
    int err = 0;

    for(l = l1; l < l2; ++l){
        int m = I[l + 1] - I[l];
        double* A_in_diag = A + I[l];
        for(int q = 0; q < m; ++q){
            double a_ij = A_in_diag[q];
            int j = I[I[l] + q];
            int m2 = I[j + 1] - I[j];
            int q2;
            for(q2 = 0; q2 < m2; ++q2){
                if(I[I[j] + q2] == l) break;
            }
            if(q2 >= m2) err++;
            else if(fabs(A[I[j] + q2] - a_ij) > eps) err++;
        }
    }
    reduce_sum<int>(p, &err, 1);
    return err;
}

void fill_A_ij(int nx, int ny, double hx, double hy, int i, int j, double* A_diag, double* A_in_diag){
    double s = hx*hy;
    if(i > 0 && i < nx && j > 0 && j < ny){
        *A_diag = s / 2;
        for(int l = 0; l < 6; ++l) A_in_diag[l] = s / 12;
    }
    else if(j == 0 && i > 0 && i < nx){
        *A_diag = 3 * s/12;
        A_in_diag[0] = 1*s/24;
        A_in_diag[1] = 1*s/24;
        A_in_diag[2] = 2*s/24;
        A_in_diag[3] = 2*s/24;
    }
    else if(j == ny && i > 0 && i < nx){
        *A_diag = 3 * s/12;
        A_in_diag[0] = 1*s/24;
        A_in_diag[1] = 2*s/24;
        A_in_diag[2] = 2*s/24;
        A_in_diag[3] = 1*s/24;
    }
    else if(i == 0 && j > 0 && j < ny){
        *A_diag = 3 * s/12;
        A_in_diag[0] = 2*s/24;
        A_in_diag[1] = 1*s/24;
        A_in_diag[2] = 1*s/24;
        A_in_diag[3] = 2*s/24;
    }
    else if(i == nx && j > 0 && j < ny){
        *A_diag = 3 * s/12;
        A_in_diag[0] = 1*s/24;
        A_in_diag[1] = 2*s/24;
        A_in_diag[2] = 2*s/24;
        A_in_diag[3] = 1*s/24;
    }
    else if(i == 0 && j == 0){
        *A_diag = 2 * s/12;
        A_in_diag[0] = 1*s/24;
        A_in_diag[1] = 1*s/24;
        A_in_diag[2] = 2*s/24;
    }
    else if(i == nx && j == ny){
        *A_diag = 2 * s/12;
        A_in_diag[0] = 1*s/24;
        A_in_diag[1] = 2*s/24;
        A_in_diag[2] = 1*s/24;
    }
    else if(i == 0 && j == ny){
        *A_diag = s/12;
        A_in_diag[0] = 1*s/24;
        A_in_diag[1] = 1*s/24;
    }
    else if(i == nx && j == 0){
        *A_diag = s/12;
        A_in_diag[0] = 1*s/24;
        A_in_diag[1] = 1*s/24;
    }
}

double F_IJ(int nx, int ny, double hx, double hy, double a, double с, int i, int j, double (*f)(double, double)){
    double q = hx * hy / 192;
    if(i > 0 && i < nx && j > 0 && j < ny){
        return q * (36*F(i, j) + 20*(F(i+0.5, j) + F(i, j-0.5) + F(i-0.5,j-0.5) + F(i-0.5, j) + F(i, j+0.5) + F(i+0.5,j+0.5))
            + 4*(F(i+0.5, j-0.5) + F(i-0.5, j-1) + F(i-1,j-0.5) + F(i-0.5,j+0.5) + F(i+0.5,j+1) + F(i+1,j+0.5))
            + 2*(F(i+1,j) + F(i, j-1) + F(i-1,j-1) + F(i-1,j) + F(i,j+1) + F(i+1,j+1))
        );
    }
    
    if(i > 0 && i < nx && j == 0){
        return q * (18*F(i,j) + 10*(F(i+0.5,j) + F(i-0.5,j)) + 20*(F(i,j+0.5) + F(i+0.5,j+0.5))
            + 4*(F(i-0.5,j+0.5) + F(i+0.5,j+1) + F(i+1,j+0.5)) + 1*(F(i-1,j) + F(i+1,j)) + 2*(F(i,j+1) + F(i+1,j+1))
        );
    }

    if(i > 0 && i < nx && j == ny){
        return q * (18*F(i, j) + 10*(F(i+0.5,j) + F(i-0.5,j)) + 20*(F(i,j-0.5) + F(i-0.5,j-0.5))
            + 4*(F(i+0.5,j-0.5) + F(i-0.5,j-1) + F(i-1,j-0.5)) + 1*(F(i-1,j) + F(i+1,j)) + 2*(F(i,j-1) + F(i-1,j-1))
        );
    }

    if(i == 0 && j > 0 && j < ny){
        return q * (18*F(i,j) + 10*(F(i,j-0.5) + F(i,j+0.5)) + 20*(F(i+0.5,j) + F(i+0.5,j+0.5))
            + 4*(F(i+0.5,j-0.5) + F(i+0.5,j+1) + F(i+1,j+0.5)) + 1*(F(i,j-1) + F(i,j+1)) + 2*(F(i+1,j) + F(i+1,j+1))
        );     
    }

    if(i == nx && j > 0 && j < ny){
        return q * (18*F(i, j) + 10*(F(i,j-0.5) + F(i,j+0.5)) + 20*(F(i-0.5,j) + F(i-0.5,j-0.5))
            + 4*(F(i-0.5,j-1) + F(i-1,j-0.5) + F(i-0.5,j+0.5)) + 1*(F(i,j-1) + F(i,j+1)) + 2*(F(i-1,j) + F(i-1,j-1))
        );
    }

    if(i == 0 && j == 0){
        return q * (12*F(i,j) + 10*(F(i+0.5,j) + F(i,j+0.5)) + 20*F(i+0.5,j+0.5)
            + 4*(F(i+1,j+0.5) + F(i+0.5,j+1)) + 1*(F(i+1,j) + F(i,j+1)) + 2*F(i+1,j+1)
        );     
    }    

    if(i == nx && j == ny){
        return q * (12*F(i,j) + 10*(F(i-0.5,j) + F(i,j-0.5)) + 20*F(i-0.5,j-0.5)
            + 4*(F(i-0.5,j-1) + F(i-1,j-0.5)) + 1*(F(i,j-1) + F(i-1,j)) + 2*F(i-1,j-1)
        );   
    }

    if(i == 0 && j == ny){
        return q * (6*F(i, j) + 10*(F(i+0.5,j) + F(i,j-0.5)) + 4*F(i+0.5,j-0.5) + F(i+1,j) + F(i,j-1));
    }

    if (i == nx && j == 0) {
        return q * (6*F(i, j) + 10*(F(i-0.5,j) + F(i,j+0.5)) + 4*F(i-0.5,j+0.5) + F(i-1,j) + F(i,j+1));      
    }
    return 1e308;
}



void fill_A(int nx, int ny, double hx, double hy, int* I, double* A, int p, int k){
    int l1, l2, i, j;
    int N = (nx + 1)*(ny + 1);    
    l1 = (N * k) / p;
    l2 = (N * (k + 1)) / p;

    for(int l = l1; l < l2; ++l){
        l_to_ij(nx, ny, i, j, l);
        double* A_diag = A + l;
        double* A_in_diag = A + I[l];
        fill_A_ij(nx, ny, hx, hy, i, j, A_diag, A_in_diag);
    }

    reduce_sum<int>(p);
}

void fill_B(int nx, int ny, double hx, double hy, double a, double c, double* B, double (*f)(double, double), int p, int k){
    int l1, l2, i, j;
    int N = (nx + 1)*(ny + 1);    
    l1 = (N * k) / p;
    l2 = (N * (k + 1)) / p;

    for(int l = l1; l < l2; ++l){
        l_to_ij(nx, ny, i, j, l);
        B[l] = F_IJ(nx, ny, hx, hy, a, c, i, j, f);
    }

    reduce_sum<int>(p);
}

void solve_r(int n, int* I, double* U, double* b, double* x, double w, int p, int k){
    int i1, i2;
    thread_rows(n, p, k, i1, i2);

    for(int i = i2 - 1; i >= i1; --i){
        int m = I[i + 1] - I[i];
        double s = 0;
        for(int q = 0; q < m; ++q){
            int j = I[I[i] + q];
            if(j >= i1 && j < i2 && j > i) s += w * x[j] * U[I[i] + q];
        }
        x[i] = (b[i] - s) / U[i];
    }    
}

void solve_l(int n, int* I, double* U, double* b, double* x, double w, int p, int k){
    int i1, i2;
    thread_rows(n, p, k, i1, i2);

    for(int i = i2 - 1; i >= i1; --i){
        int m = I[i + 1] - I[i];
        double s = 0;
        for(int q = 0; q < m; ++q){
            int j = I[I[i] + q];
            if(j >= i1 && j < i2 && j < i) s += w * x[j] * U[I[i] + q];
        }
        x[i] = (b[i] - s) / U[i];
    }    
}
