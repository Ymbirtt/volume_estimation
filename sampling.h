/* This file was automatically generated.  Do not edit! */
typedef struct __point__{
    float* x;
    int n;
} point;
typedef int (*inclusion_oracle)(point);
typedef float (*density)(point);
void rejection_sample_dist(point p,density f,point lower,point upper,float f_max);
float funky_dist(point p);
int square(point p);
void grid_walk(point start,int steps,float delta,inclusion_oracle in);
void ball_walk(point start,int steps,float delta,inclusion_oracle in);
void hit_and_run(point start,int steps,float epsilon,float diameter,inclusion_oracle in);
void hit_and_run_dist(point start,int steps,float epsilon,float diameter,inclusion_oracle in,density f);
float rejection_sample(density f,point dir,point start,float lbound,float a,float b,float m);
float find_argval_left(density f,point dir,point start,float epsilon,float lbound,float ubound);
float find_argval_right(density f,point dir,point start,float epsilon,float lbound,float ubound);
float find_argmax(density f,point dir,point start,float epsilon,float lbound,float ubound);
float min(float a,float b);
float max(float a,float b);
float boundary_distance(point start,point dir,float epsilon,float diameter,inclusion_oracle in);
void move_point(point dest,point src);
float sq_dist(point p1,point p2);
void sub_point(point dest,point src1,point src2);
void addmul_point(point dest,point src1,float lambda,point src2);
void add_point(point dest,point src1,point src2);
void random_dir(point v,float r);
void gaussians(float *dest1,float *dest2);
point copy_point(point p);
float uniform(float a,float b);
void print_point(point p);
extern unsigned long int ZIGGURAT_SEED;
extern float ZIGGURAT_WN[128];
extern int ZIGGURAT_KN[128];
extern float ZIGGURAT_FN[128];
