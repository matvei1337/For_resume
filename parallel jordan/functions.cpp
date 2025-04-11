#include"functions.h"

double get_full_time()
{
    struct timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec / 1.e6;
}

double get_cpu_time()
{
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec /1.e6;
}


rows
void get_block(double* a, double* block, int n, int m, int row, int col){
    int i = 0, j = 0;
    int k = n / m;
    int l = n - k*m;
    int height = (row < k ? m : l);
    int width = (col < k ? m : l);

    double* bl = a + m * (row * n + col);
    for(i = 0; i < height; i ++){
        for(j = 0; j < width; j ++){
            block[i*width+j] = bl[i*n+j];
        }
    }
}

void get_block2(double *a, double *block, int n, int m, int col){
    int i = 0, j = 0;
    int k = n / m;
    int l = n - k * m;
    int width = col < k ? m : l;

    double *bl = a + col * m;
    for (i = 0; i < m; ++i){
        for (j = 0; j < width; ++j){
            block[i*width+j] = bl[i*n+j];
        }
    }
}




void set_block(double* a, double* block, int n, int m, int row, int col){
    int i = 0, j = 0;
    int k = n / m;
    int l = n - k*m;
    int height = (row < k ? m : l);
    int width = (col < k ? m : l);

    double* bl = a + m * (row * n + col);
    for(i = 0; i < height; i ++){
        for(j = 0; j < width; j ++){
            bl[i*n+j] = block[i*width+j];
        }
    }
}

void set_block2(double *a, double *block, int n, int m, int row){
    int i = 0, j = 0;
    int k = n / m;
    int l = n - k * m;
    int cols = row < k ? m : l;

    double *bl = a + row * m;
    for (i = 0; i < m; ++i){
        for (j = 0; j < cols; ++j){
            bl[i*n+j] = block[i*cols+j];
        }
    }
}


void copy_matrix(double *a, double *b, int n, int m){
    int i, j;
    for (i = 0; i < n; ++i){
        for (j = 0; j < m; ++j) a[i*m+j] = b[i*m+j];
    }
}




double f(int s, int n, int i, int j){
    if(s == 1) return n - std::max(i,j);
    else if(s == 2) return std::max(i,j)+1;
    else if(s == 3) return fabs(i-j);
    else return 1.0/(1.0*(i+j+1));
}



int read(double *a, int n, const char *name) {
    FILE* file;
    double buf;
    if(!(file = fopen(name, "r"))) {
        return -1;
    }
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(fscanf(file, "%lf", &buf) != 1) {
            fclose(file);
            return -1;
            }
        a[i * n + j] = buf;
        }
    }
    fclose(file);
    return 0;
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


void init_a(double* a,int n,int m, int s, int k, int p) {
    for (int j = k*m; j < n; j += p*m) {
        int h = (j + m < n) ? m : n - j;
        for (int i = 0; i < n; ++i) {
            for (int v = 0; v < h; ++v) {
                a[i*n + j + v] = f(s, n, i, j + v);
            }   
        }
    }
    reduce_sum<int>(p);
}

void init_b(double* a,double* b,int n,int m, int k, int p) {
    for (int i = k*m; i < n; i += p*m) { 
        int h = (i + m < n) ? m : n - i;
        for (int t = i; t < i + h; ++t) {
            double sum = 0;
            for (int k = 0; k <= (n-1) / 2; ++k) {
                sum += a[n * t + 2*k];
            }
            b[t] = sum; 
        }
    }   
    reduce_sum<int>(p);
}



double norma1(double* a, double* x, double* b, int n){
    int i = 0, j = 0;
    double norm1 = 0.0, norm2 = 0.0, s;
    for(i = 0; i < n; i++){
        s = 0.0;
        for(j = 0; j < n; j++){
            s += a[i*n+j]*x[j];
        }
        norm1 += fabs(s-b[i]);
        norm2 += fabs(b[i]);
    }
    return norm1/norm2;
}

double norma2(double* x, int n){
    int i = 0;
    double norm = 0.0;
    for(i = 0; i < n; i++){
        norm += fabs(x[i] - (i+1)%2);
    }
    return norm;
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

void init_a_prosto(double* a, int n, int s){
    int i = 0, j = 0;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            a[i*n+j] = f(s,n,i,j);
        }
    }
}


void init_b_prosto(double* b, double* a, int n){
    int i = 0, j = 0;
    double s;
    for(i = 0; i < n; i++){
        s = 0.0;
        for(j = 0; j < n; j+=2){
            s += a[i*n+j]; 
        }
        b[i] = s;
    }
}
