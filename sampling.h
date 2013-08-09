/* This file was automatically generated.  Do not edit! */
typedef struct __point__ __point__;
struct __point__ {
    double* x;
    int n;
};
typedef struct __point__ point;
typedef double(*density)(point);
void rejection_sample_dist(point p,density f,point lower,point upper,double f_max);
double funky_dist(point p);
int square(point p);
typedef int(*inclusion_oracle)(point);
void grid_walk(point start,int steps,double delta,inclusion_oracle in);
void ball_walk(point start,int steps,double delta,inclusion_oracle in);
void hit_and_run(point start,int steps,double epsilon,double diameter,inclusion_oracle in);
void hit_and_run_dist(point start,int steps,double epsilon,double diameter,inclusion_oracle in,density f);
double rejection_sample(density f,point dir,point start,double lbound,double a,double b,double m);
double find_argval_left(density f,point dir,point start,double epsilon,double lbound,double ubound);
double find_argval_right(density f,point dir,point start,double epsilon,double lbound,double ubound);
double find_argmax(density f,point dir,point start,double epsilon,double lbound,double ubound);
double max(double a,double b);
double boundary_distance(point start,point dir,double epsilon,double diameter,inclusion_oracle in);
void move_point(point dest,point src);
double sq_dist(point p1,point p2);
void sub_point(point dest,point src1,point src2);
void addmul_point(point dest,point src1,double lambda,point src2);
void add_point(point dest,point src1,point src2);
void random_dir(point v,double r);
void gaussians(double *dest1,double *dest2);
point copy_point(point p);
double uniform(double a,double b);
void print_point(point p);
#define INTERFACE 0
