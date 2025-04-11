#pragma once
#include<stdio.h>
#include<sched.h>
#include<sys/sysinfo.h>
#include<sys/resource.h>
#include<sys/time.h>
#include<unistd.h>
#include<pthread.h>
#include<time.h>
#include<cmath>
#include <cstring>
#define EPS 1e-15

struct Args{	
	double* a = nullptr;
	double* b = nullptr;
	double* x = nullptr;
	int n = 0;
	int m = 0;
	int p = 0;
	int k = 0;
	int r = 0;
    int s = 0; 
    const char* filename = nullptr;

	double t_cpu = 0;
	int status = 0;
	
};



void multiply(double* A, double* B,double* C, int na, int ma,int nb, int mb);
void swap_rows(double* a,int i, int j);
int reverse(double* C,double* D, int n,double norm);
int solve(double* a, double* b, double* x, int n, int m, int p, int k,double* b1,double* b2,double* b3,double norm,double* copy_row,double* copy_b);
void* thread_func(void* ptr);


void swap_block_rows(double* a,int n,int m,int i, int j, double* C, double* D);
void swap_block_columns(double* a,int n,int m,int i, int j, double* C, double* D);

double get_full_time();
double get_cpu_time();

double f(int s, int n, int i, int j);
int read(double* a, int n, const char *name);
void print(double* a, int l, int n, int r);

void init_a(double* a,int n,int m, int s, int k, int p);
void init_b(double* a,double* b,int n,int m, int k, int p);

void init_a_prosto(double* a, int n, int s);
void init_b_prosto(double* b, double* a, int n);


double norma1(double* a, double* x, double* b, int n);
double norma2(double* x, int n);
double matrix_norma(double* a, int n);

void get_block(double* a,double* block,int n, int m,int i, int j);
void set_block(double* a,double* block,int n, int m,int i, int j);
void get_block2(double *a, double *block, int n, int m, int row);
void set_block2(double *a, double *block, int n, int m, int row);

void copy_matrix(double *a, double *b, int n, int m);


template<class T>
void reduce_sum(int total_threads, T* a=nullptr, int n=0){
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t cond_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t cond_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0;
	static int t_out = 0;
	static T* res = nullptr;
	int i;

	pthread_mutex_lock(&mutex);
	if(res == nullptr) res = a;
	else{
		for(i = 0; i < n; i++){
			res[i] += a[i];
		}
	}
	t_in++;
	if(t_in >= total_threads){
		t_out = 0;
		pthread_cond_broadcast(&cond_in);
	}
	else{
		while(t_in < total_threads) pthread_cond_wait(&cond_in, &mutex);
	}
	if(res != a){
		for(i = 0; i < n; i++){
			a[i] = res[i];
		}
	}
	t_out++;
	if(t_out >= total_threads){
		t_in = 0;
		res = nullptr;
		pthread_cond_broadcast(&cond_out);
	}
	else{
		while(t_out < total_threads) pthread_cond_wait(&cond_out, &mutex);
	}
	pthread_mutex_unlock(&mutex);
}
