#include"functions.h"
using namespace std;

int main(int argc, char *argv[]){
	double* a;
    double* x;
    double* x1;
    double* x2;

	int n, m, k,res;
    char* filename;
    int its = 0;
    double norma_mat,eps;
    double r1 = -1.0, r2 = -1.0, t1 = 0.0, t2 = 0.0;

	if((argc != 5 && argc != 6) || sscanf(argv[1],"%d",&n) != 1 || sscanf(argv[2],"%d",&m) != 1 || sscanf(argv[3],"%lf",&eps) != 1 || sscanf(argv[4],"%d",&k) != 1){
        printf("Usage: %s n m eps k [file]\n", argv[0]);
        return -1;
    }
    filename = nullptr;
    if(argc == 6){
    	filename = argv[5];
    	if(k!=0) return -1;
    }
    
    a = new double[n*n];
    if(!a){
        delete[] a;
    	cout << "Bad alloc\n";
    	return -1;
    }
    int flag = 0;
    if(argc == 5) init_a(a, n, k);
    else{
        flag = read(a, n, filename);
        if(flag != 1){
            cerr << "Error file\n";
            delete[] a;
            return -1;
        }
    }
    print(a, n, n, m);
    printf("\n");
    norma_mat = matrix_norma(a,n);
    eps = eps*norma_mat;

    x = new double[n];
    x1 = new double[n];
    x2 = new double[n];

    t1 = clock();
    three_diag(a,n,x,eps);
    t1 = (clock()-t1)/CLOCKS_PER_SEC;

    t2 = clock();
    res = solve(a,n,x,eps,x1,x2,its);
    t2 = (clock()-t2)/CLOCKS_PER_SEC;
    print(a, n, n, m);


    if(flag == 0) init_a(a, n, k);
    else read(a, n, filename);

    if(res < 0) printf("Cannot solve");
    else{
        r1 = norma1(a, x, n);
        r1 = r1/norma_mat;
        r2 = norma2(a, x, n);
        r2 = r2/norma_mat;
        print(x, 1, n, m);
    }

    printf("\n");
    printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
    argv[0], r1, r2, its, its/n, t1, t2);
    
    delete[] a;
    delete[] x;
    delete[] x1;
    delete[] x2;
	return 0;
}
