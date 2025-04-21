#include <string.h>
#include "gui.h"
#include "solver.h"

static double length(double x, double y){
    return sqrt(x*x + y*y);
}

void Gui_data::read_data(char* argv[]){
    sscanf(argv[1], "%lf", &a);
    sscanf(argv[2], "%lf", &b);
    sscanf(argv[3], "%lf", &c);
    sscanf(argv[4], "%lf", &d);
    sscanf(argv[5], "%d", &nx);
    sscanf(argv[6], "%d", &ny);
    sscanf(argv[7], "%d", &mx);
    sscanf(argv[8], "%d", &my);
    sscanf(argv[9], "%d", &m);
    sscanf(argv[10], "%lf", &eps);
    sscanf(argv[11], "%d", &maxit);
    sscanf(argv[12], "%d", &p);
}

void Gui_data::realloc_data(){
    if(A){
        delete[] A;
    }
    if(I){
        delete[] I;
    }
    if(B){
        delete[] B;
    }
    if(x){
        delete[] x;
    }
    if(u){
        delete[] u;
    }
    if(v){
        delete[] v;
    }
    if(r){
        delete[] r;
    }

    I = nullptr;
    A = nullptr;
    alloc_msr(nx, ny, &A, &I);
    init_red_sum(p);
    int N = (nx + 1) * (ny + 1);
    B = new double[N];
    x = new double[N];
    r = new double[N];
    u = new double[N];
    v = new double[N];

    fill_I(nx, ny, I);
    memset(x, 0, N*sizeof(double));
}

double Gui_data::find_min_max(double& abs_min, double& abs_max){
    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double value1, value2;
    for(int i = 0; i < mx; ++i){
        for(int j = 0; j < my; ++j){
            value1 = f(a + (i + 1.0/3.0)*hx, c + (j + 2.0/3.0) * hy);
            value2 = f(a + (i + 2.0/3.0)*hx, c + (j + 1.0/3.0) * hy);
            abs_min = std::min(abs_min, value1);
            abs_max = std::max(abs_max, value1);
            abs_min = std::min(abs_min, value2);
            abs_max = std::max(abs_max, value2);
        }
    }
    return std::max(abs_min, abs_max);
}

double Gui_data::pf(Point p){
    int i, j, l;
    double hnx = (b - a) / nx;
    double hny = (d - c) / ny;
    double w1, w2;

    i = (p.x - a) / hnx;
    j = (p.y - c) / hny;
    ij_to_l(nx, ny, i, j, l);

    if(std::fabs((p.x - a) - i * hnx) < EPS && fabs((p.y - c) - j * hny) < EPS){
        return x[l];
    }

    double d2 = length((p.x - a) - (i + 1) * hnx, (p.y - c) - j*hny);
    double d3 = length((p.x - a) - i * hnx, (p.y - c) - (j + 1)*hny);

    if(d2 < d3){
        w1 = ((i + 1) * hnx - (p.x - a)) / hnx;
        w2 = (p.y - c - j * hny) / hny;
        return  w1 * x[l] + w2 * x[l + 1 + (nx + 1)] + (1 - w1 - w2) * x[l + 1];
    }

    w1 = (p.x - a - i * hnx) / hnx;
    w2 = ((j + 1) * hny - (p.y - c)) / hny;
    return w1 * x[l + 1 + (nx + 1)] + w2 * x[l] + (1 - w1 - w2) * x[l + (nx + 1)];
}

void Gui_data::pfind_min_max(double& abs_min, double& abs_max){
    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double value1, value2;
    for(int i = 0; i < mx; ++i){
        for(int j = 0; j < my; ++j){
            value1 = pf({a + (i + 1.0/3.0)*hx, c + (j + 2.0/3.0) * hy});
            value2 = pf({a + (i + 2.0/3.0)*hx, c + (j + 1.0/3.0) * hy});
            abs_min = std::min(abs_min, value1);
            abs_max = std::max(abs_max, value1);
            abs_min = std::min(abs_min, value2);
            abs_max = std::max(abs_max, value2);
        }
    }
}

void Gui_data::residual_min_max(double& abs_min, double& abs_max){
    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double f_val, pf_val, value1, value2;
    for(int i = 0; i < mx; ++i){
        for(int j = 0; j < my; ++j){
            f_val = f(a + (i + 1.0/3.0)*hx, c + (j + 2.0/3.0) * hy);
            pf_val = pf({a + (i + 1.0/3.0)*hx, c + (j + 2.0/3.0) * hy});
            value1 = std::fabs(f_val - pf_val);

            f_val = f(a + (i + 2.0/3.0)*hx, c + (j + 1.0/3.0) * hy);
            pf_val = pf({a + (i + 2.0/3.0)*hx, c + (j + 1.0/3.0) * hy});
            value2 = std::fabs(f_val - pf_val);

            abs_min = std::min(abs_min, value1);
            abs_max = std::max(abs_max, value1);
            abs_min = std::min(abs_min, value2);
            abs_max = std::max(abs_max, value2);
        }
    }
}

Gui_data::~Gui_data(){
    if(A){
        delete[] A;
    }
    if(I){
        delete[] I;
    }
    if(B){
        delete[] B;
    }
    if(x){
        delete[] x;
    }
    if(u){
        delete[] u;
    }
    if(v){
        delete[] v;
    }
    if(r){
        delete[] r;
    }
}