#include "calc_func.h"

// Chebyshev method

void CalculateAlpha(int n, double a, double b, const QVector<double>& x, const QVector<double>& y, QVector<double>* alpha) {

    for (int j = 0; j < n; ++j) {
        FillAlpha(n, a, b, x[j], y[j], alpha);
    }

    (*alpha)[0] /= n;
    for (int i = 1; i < n; ++i) {
        (*alpha)[i] *= 2;
        (*alpha)[i] /= n;
    }    
}

void FillAlpha(int n, double a, double b, double x, double y, QVector<double>* alpha) {
    double T0 = 1;
    double T1 = (2*x - (b + a)) / (b - a);
    (*alpha)[0] += y;
    if (n == 1) {
        return;
    }

    (*alpha)[1] += T1 * y;
    double var = T1;
    double res = 0;
    for (int k = 2; k < n; ++k) {
        res = 2*var*T1 - T0;
        T0 = T1;
        T1 = res;
        (*alpha)[k] += res * y;
    }
}

double Pf(int n, double x, double a, double b, const QVector<double>& alpha) {
    double T0 = 1;
    double T1 = (2*x - (b + a)) / (b - a);
    double Pf_x = alpha[0] * T0;
    if (n == 1) {
        return Pf_x;
    } else if (n == 2) {
        return Pf_x + alpha[1] * T1;
    }

    Pf_x += alpha[1] * T1;
    double var = T1;
    double res = 0;
    for (int k = 2; k < n; ++k) {
        res = 2*var*T1 - T0;
        T0 = T1;
        T1 = res;
        Pf_x += alpha[k] * res;
    }

    return Pf_x;
}

// Spline method

double difference(double yi, double yj, double xi, double xj) {
    return (yj - yi) / (xj - xi);
}

void construct_matrix(int n, double (*df) (double), const QVector<double>& x, const QVector<double>& y,
                      double* left_d, double* main_d, double* right_d, double* b) {

    main_d[0] = 1; main_d[n - 1] = 1;
    right_d[0] = 0; left_d[n - 2] = 0;
    for (int i = 1; i < n - 1; ++i) {
        left_d[i - 1] = x[i+1] - x[i];
        main_d[i] = 2*(x[i+1]-x[i-1]);
        right_d[i] = x[i] - x[i-1];
    }


    b[0] = df(x[0]);
    b[n-1] = df(x[n-1]);
    for (int i = 1; i < n - 1; ++i) {
        b[i] = 3*difference(y[i-1], y[i], x[i-1], x[i])*(x[i+1]-x[i]) + 3*difference(y[i], y[i+1], x[i], x[i+1])*(x[i]-x[i-1]);
    }
}

void solution(int n, double* main_d, double* left_d, double* right_d, double* b, double* d) {
    for (int i = 0; i < n - 1; ++i) {
        right_d[i] /= main_d[i];
        b[i] /= main_d[i];

        main_d[i + 1] -= right_d[i] * left_d[i];
        b[i + 1] -= b[i] * left_d[i];
    }

    d[n - 1] = b[n - 1] / main_d[n - 1];

    for (int i = n - 2; i >= 0; --i) {
        d[i] = b[i] - right_d[i] * d[i + 1];
    }
}

int calc_i(double pt, int n, const QVector<double>& x) {
    auto lower = std::lower_bound(x.begin(), x.end(), pt);
    if (lower == x.begin()) {
        return 0;
    } else if (lower == x.end()) {
        return n - 2;
    }

    return std::distance(x.begin(), lower) - 1;
}

void calc_coeff(int i, const QVector<double>& x, const QVector<double>& y, double* d, Spline* spline) {
    spline->c1 = y[i];
    spline->c2 = d[i];
    spline->c3 = (3*difference(y[i], y[i+1], x[i], x[i+1]) - 2*d[i] - d[i+1]) / (x[i+1] - x[i]);
    spline->c4 = (d[i] + d[i+1] - 2*difference(y[i], y[i+1], x[i], x[i+1])) / ((x[i] - x[i+1]) * (x[i] - x[i+1]));
}

double Sf(double pt, int n, const QVector<double>& x, const QVector<double>& y, double* d) {
    int i = calc_i(pt, n, x);
    Spline spline;
    calc_coeff(i, x, y, d, &spline);
    return spline.c1 + spline.c2*(pt-x[i]) + spline.c3*(pt-x[i])*(pt-x[i]) + spline.c4*(pt-x[i])*(pt-x[i])*(pt-x[i]);
}
