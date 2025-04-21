#include <QVector>

// Chebyshev method
void CalculateAlpha(int n, double a, double b,
    const QVector<double>& x, const QVector<double>& y, QVector<double>* alpha);
void FillAlpha(int n, double a, double b, double x, double y, QVector<double>* alpha);
double Pf(int n, double x, double a, double b, const QVector<double>& alpha);

// Spline method

struct Spline {
    double c1;
    double c2;
    double c3;
    double c4;
};

double difference(double yi, double yj, double xi, double xj);
void construct_matrix(int n, double (*df) (double), const QVector<double>& x, const QVector<double>& y,
                      double* left_d, double* main_d, double* right_d, double* b);
void solution(int n, double* main_d, double* left_d, double* right_d, double* b, double* d);
int calc_i(double pt, int n, const QVector<double>& x);
void calc_coeff(int i, const QVector<double>& x, const QVector<double>& y, double* d, Spline* spline);
double Sf(double pt, int n, const QVector<double>& x, const QVector<double>& y, double* d);
