#include"sampling.h"
#include"mt.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<time.h>

#define PI (3.1415926535897932384626338)

typedef int (*inclusion_oracle)(point);

double SCALE_FAC = HUGE_VAL;
//double GAMMA = HUGE_VAL;
double PENCIL_HEIGHT = HUGE_VAL;
inclusion_oracle PENCIL_ORACLE = NULL;
inclusion_oracle BALL_ORACLE;
double BALL_RADIUS;
double LAST_BALL_RADIUS;


int double_fac(int n){
    int i;
    int res=1;

    for (i=3;i<=n; i+=2){
        res *= i;
    }
    return res;
}

int fac(int n){
    int i;
    int res=1;

    for (i=2;i<=n;i++){
        res*=i;
    }
    return res;
}

double int_pow(double x, int y){
    double result = 1;

    while (y>0){
        if (y&1){
            result *= x;
        }

        x *= x;
        y = y>>1;
    }

    return result;
}

//http://en.wikipedia.org/wiki/Volume_of_an_n-ball
double volume_nball(double radius, int n){

    if (n&1){
        //formula for odd n
        return (1<<(((n-1)/2)+1)) * int_pow(PI,(n-1)/2) * int_pow(radius, n) / ((double) double_fac(n));
    } else {
        //formula for even n
        return int_pow(PI,n/2)*int_pow(radius,n)/((double)fac(n/2));
    }
}

//Initialise an array of values for the vector a
double* init_as(const int m, const int n){
    int i;
    double* a = malloc((m+1)*sizeof(double));
    double factor = (1-1/sqrt(n));

    a[1] = 6*n;

    for(i=2; i<=m; i++){
        a[i] = factor*a[i-1];
    }

    return a;
}

double* init_gammas(int m, int n, double* a){
    double* gamma = malloc(m*sizeof(double));
    int i;
    double sqrt_n = sqrt(n);

    for(i=0; i<m; i++){
        gamma[i] = max(1.0, a[i]/sqrt_n);
    }

    return gamma;
}

double exp_density(point p){
    return exp(SCALE_FAC * p.x[p.n-1]);
}

//Returns 0 iff p is not inside the pencil constructed over the pencil oracle
int in_pencil(point p){
    double sq_sum = 0;
    int i;
    
    //printf("Pencil!\n");
    
    if (p.x[p.n-1] < 0 || p.x[p.n-1] > PENCIL_HEIGHT){
        //printf("Failed from height\n");
        return 0;
    }
    
    for(i=0; i<p.n-1; i++){
        sq_sum += p.x[i] * p.x[i];
    }

    if(sq_sum > (p.x[p.n-1])*(p.x[p.n-1])){
        //printf("Failed from radius\n");
        return 0;
    }
    
    if (PENCIL_ORACLE(p)){
        return 1;
    } else {
        //printf("NOPE\n");
        return 0;
    }
}


double sq_norm(point p){
    int i;
    double res = 0;

    for(i=0; i<p.n; i++){
        res += p.x[i]*p.x[i];
    }
    return res;
}

int k_and_ball_i(point p){
    if (sq_norm(p) > BALL_RADIUS*BALL_RADIUS) return 0;
    return BALL_ORACLE(p);
}

int k_and_last_ball(point p){
    if (sq_norm(p)>LAST_BALL_RADIUS*LAST_BALL_RADIUS) return 0;
    return BALL_ORACLE(p);
}

double estimate_pencil(inclusion_oracle o, int n, double d, double epsilon, double theta){
    double vol = volume_nball(1,n);
    double radius_ratio = exp((double) 1/n);
    int stages = (int) n*ceil(log(d)) +1;
    int threads = 40000; //(int) (400 * n * log(n) / (epsilon * epsilon));
    int i,j,k;
    int points_in;
    point p;


    p.n = n;
    p.x = malloc(n*sizeof(double));

    for(i=0; i<n; i++) p.x[i] = 0;
    
    printf("Each stage uses %d threads\n", threads);

    PENCIL_ORACLE = o;
    PENCIL_HEIGHT = d;
    
    printf("Pencil height = %lf\n", PENCIL_HEIGHT);
    
    BALL_ORACLE = &in_pencil;
    LAST_BALL_RADIUS = 1;
    BALL_RADIUS = LAST_BALL_RADIUS * radius_ratio;

    for(i=1; i<stages; i++){
        printf("Stage %d of %d      \r", i, stages);
        fflush(stdout);
        points_in = 0;

        //BALL_RADIUS = exp((double) i/(n));
        //LAST_BALL_RADIUS = exp((double) (i-1)/(n));

        //printf("Querying ball of radius %lf\n", BALL_RADIUS);
        for(j=0; j<threads; j++){
            for(k=0; k<n-1; k++) p.x[k] = 0;
            p.x[n-1] = min(BALL_RADIUS/2.0, PENCIL_HEIGHT/2.0);
            hit_and_run(p,1<<5, theta, BALL_RADIUS*2, &k_and_ball_i);
            //print_point(p);
            //printf("\n");
            if (k_and_last_ball(p)){
                points_in++;
            }
        }

        LAST_BALL_RADIUS = BALL_RADIUS;
        BALL_RADIUS *= radius_ratio;
        
        printf("\n%d of %d points inside last ball\n", points_in, threads);
        vol *= (double) threads/points_in;
    }

    printf("\n");
    
    printf("vol = %lf\n", vol);
    
    return vol;

}

/*
double estimate_pencil(inclusion_oracle o, int n, double d, double epsilon){
    int m           = (int) (2*ceil(sqrt(n-1)*log((n-1)/epsilon)));
    int k           = 250; //(int) (8/(epsilon*epsilon)) * sqrt(n) * log(n/epsilon);
    double  delta   = pow(n-1,-10)*epsilon*epsilon;
    double* a       = init_as(m+1,n-1);
    double* gamma   = init_gammas(m+1, n-1, a);

    double  z       = fac(n-1) * volume_nball(1,n-1) * pow(6*(n-1),-(n));
    double  sum     = 0;
    point*  big_x   = malloc(k*sizeof(point));

    int i,j;

    PENCIL_HEIGHT = 2*d;
    PENCIL_ORACLE = o;

    printf("z starts at %lf\n", z);
    printf("ball volume is %lf\n", volume_nball(1,n-1));


    for(i=0; i<k; i++){
        big_x[i].x = malloc((n)*sizeof(double));
        big_x[i].n = n;
        for(j=0; j<n-1; j++){
            big_x[i].x[j] = 0;
        }
        big_x[i].x[n-1] = 1;
    }

    for(i=1; i<=m; i++){
        printf("Phase %d of %d        \n", i, m);
        //if (i!=1){
            SCALE_FAC = -a[i];//gamma[i];
            //GAMMA = gamma[i];
        //} else {
        //    SCALE_FAC = 1;
        //    GAMMA = 1;
        //}
        
        //printf("GAMMA = %lf\n", GAMMA);

        for(j=0; j<k; j++){
            //printf("Thread %d of %d\r", j, k);
            printf("Pre-transformed point");
            print_point(big_x[j]);
            printf("\n");
            //if (i!=1) big_x[j].x[n-1] = big_x[j].x[n-1]*gamma[i];
            hit_and_run_dist(big_x[j], 1<<5, delta, 4*d, &in_pencil, &exp_density);
            //if (i!=1) big_x[j].x[n-1] = big_x[j].x[n-1]/gamma[i];
            sum += exp((a[i]-a[i+1]) * big_x[j].x[n-1]);
        }
        z *= sum/k;
        sum = 0;
    }
    //GAMMA = 1;

    printf("pencil has volume %lf\n", z);
    return z;
}
*/

double estimate_volume(inclusion_oracle o, int n, double d, double epsilon){
    double pencil_volume;
    int i,j;
    point p;
    int points_in = 0;
    int total_points = (int) (1/(epsilon*epsilon));

    p.n = n+1;
    p.x = malloc(n*sizeof(double));

    for (i=0; i<p.n; i++){
        p.x[i] = 0;
    }

    pencil_volume = estimate_pencil(o, n+1, 2*d, epsilon, 0.00001);
    //pencil_volume = 23;
    //getchar();
    
    PENCIL_HEIGHT = 2*d;
    PENCIL_ORACLE = o;
    
    for (i=0; i < total_points; i++){
        //generate random point in o x [0,2d]

        for (j=0; j<p.n; j++){
            p.x[i] = uniform(-1,1);
        }

        //hit_and_run(p, (int) n*n*d*d, epsilon, 4*d, o);

        p.x[p.n-1] = uniform(0,2*d);

        if (in_pencil(p)){
            points_in++;
        }
    }

    printf("Cylinder is %lf the size of pencil\n", (double)total_points/points_in);

    printf("By cheating, I know the volume of the cylinder is %lf\n", 2*d*8);

    printf("The volume of the pencil should be %lf\n", 2*d*8*((double) points_in/total_points));

    return ((double) total_points/points_in) * pencil_volume / (2*d);
}

//The square [0,1]^5
/*
int square(point p){
    assert(p.n == 5);
    if(p.x[0] >= 0.0 && p.x[0] <= 1.0 &&
       p.x[1] >= 0.0 && p.x[1] <= 1.0 &&
       p.x[2] >= 0.0 && p.x[2] <= 1.0 &&
       p.x[3] >= 0.0 && p.x[3] <= 1.0 &&
       p.x[4] >= 0.0 && p.x[4] <= 1.0){
        return 1;
    } else {
        return 0;
    }
}
*/
int main (void){
    init_genrand(0);
    printf("%lf\n", estimate_volume(&square, 3, sqrt(3), 0.01));
}