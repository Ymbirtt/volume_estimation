#include"sampling.h"
#include"mt.h"
#include"ziggurat.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<time.h>


#define PI (3.1415926535897932384626338)

double BALL_RADIUS;
double LAST_BALL_RADIUS;
inclusion_oracle BALL_ORACLE;

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

double estimate_volume(inclusion_oracle o, int n, double d, double epsilon, double theta){
    double vol = volume_nball(1,n);
    double radius_ratio = exp((double) 1/n);
    int stages = (int) n*ceil(log(d)) +1;
    int threads = (int) (400 * n * log(n) / (epsilon * epsilon));
    int i,j,k;
    int points_in;
    point p;


    p.n = n;
    p.x = malloc(n*sizeof(double));

    for(i=0; i<n; i++) p.x[i] = 0;
    
    printf("Each stage uses %d threads\n", threads);

    BALL_ORACLE = o;
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
            for(k=0; k<n; k++) p.x[k] = 0;
            hit_and_run(p,1<<5, theta, BALL_RADIUS*2, &k_and_ball_i);
            if (k_and_last_ball(p)){
                points_in++;
            }
        }

        LAST_BALL_RADIUS = BALL_RADIUS;
        BALL_RADIUS *= radius_ratio;
        
        //printf("%d of %d points inside last ball\n", points_in, threads);
        vol *= (double) threads/points_in;
    }

    printf("\n");
    
    return vol;

}

int main(void){
    int i;
    double vol;
    double total = 0;
    
    ZIGGURAT_SEED = 123456789;
    
    init_genrand(clock());
    
    r4_nor_setup(ZIGGURAT_KN, ZIGGURAT_FN, ZIGGURAT_WN);
    
    for(i=0; i<1; i++){
        vol = estimate_volume(&square, 3, sqrt(3), 0.1, 0.00001);
        printf("%lf\n", vol);
        total += vol;
    }
    printf("%lf\n", total/1);
}