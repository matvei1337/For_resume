#pragma once
#include<stdio.h>
#include<iostream>
#include<time.h>
#include<cmath>

void three_diag(double *a, int n, double *x,double eps);
int solve(double* a, int n, double* x, double eps, double *x1, double *x2, int& its);


double f(int k, int n, int i, int j);
int read(double* a, int n, const char *name);
void print(double* a, int l, int n, int r);
void init_a(double* a, int n, int k);

double norma1(double* a, double* x, int n);
double norma2(double* a, double* x, int n);
double matrix_norma(double* a, int n);
