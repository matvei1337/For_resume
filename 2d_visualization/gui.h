#pragma once
#define EPS 1e-15

class Point{
public:
    double x;
    double y;
    Point(): x(0), y(0){}
    Point(double x, double y): x(x), y(y){}
};

class Gui_data{
public:
    double *A = nullptr;
    int *I = nullptr;
    double *B = nullptr;
    double *x = nullptr;
    double *u = nullptr;
    double *v = nullptr;
    double *r = nullptr;

    double a;
    double b;
    double c;
    double d;
    double eps;
    int nx;
    int ny;
    int mx;
    int my;
    int m;
    int maxit;
    int p;
    double (*f)(double, double);

    void read_data(char* argv[]);
    void realloc_data();
    double find_min_max(double& abs_min, double& abs_max);
    double pf(Point p);
    void pfind_min_max(double& abs_min, double& abs_max);
    void residual_min_max(double& abs_min, double& abs_max);
    ~Gui_data();
};