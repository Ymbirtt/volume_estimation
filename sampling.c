//Delicious Mersenne Twister PRNG
#include"mt.h"

//I don't see why I should have to justify any of these includes.
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<time.h>

//The probability of remaining stationary in a random walk
#define WAIT_PROB (0)
#define PI (3.1415926535897932384626338)
#define E  (2.718281828459045235360)

typedef struct __point__{
    double* x;
    int n;
} point;

//An inclusion oracle should take a point and return 0 iff the point is not
//inside the shape
typedef int (*inclusion_oracle)(point);

void print_point(point p){
    int i;
    for (i=0; i<p.n-1; i++){
        printf("%lf,", p.x[i]);
    }
    printf("%lf", p.x[i]);
}

//Returns a continuous, uniormly distributed random-ish number between a and b
//seeding function should be called beforehand
double uniform(double a, double b){
    return genrand_real1()*(b-a) + a;
}

//Creates a freshly allocated copy of the point p. Use with caution and collect
//your garbage
point copy_point(point p){
    point new_point;
    int i;

    new_point.n = p.n;
    new_point.x = malloc(p.n*sizeof(double));
    for(i=0; i<p.n; i++){
        new_point.x[i] = p.x[i];
    }
    return new_point;
}

//Stores two normal gaussian variables in dest1 and dest2
void gaussians(double* dest1, double* dest2){
    double rand1 = uniform(0,1);
    double rand2 = uniform(0,1);
    
    *dest1 = sqrt(-2*log(rand1)) * cos(2*PI*(rand2));
    *dest2 = sqrt(-2*log(rand1)) * sin(2*PI*(rand2));
}

//Store a uniformly selected random vector of magnitude r in v
//http://mathworld.wolfram.com/HyperspherePointPicking.html
//Probably worth bringing a Ziggurat library in, time permitting.
//http://people.sc.fsu.edu/~jburkardt/c_src/ziggurat/ziggurat.html
void random_dir(point v, double r){
    int i;
    double rand1, rand2;
    double square_sum = 0;

    for (i=0; i<v.n/2; i++){
        gaussians(&rand1, &rand2);
        
        v.x[2*i] = rand1;
        v.x[2*i+1] = rand2;

        square_sum += rand1*rand1 + rand2*rand2;
        
    }

    if (2*i != v.n){
        gaussians(&rand1, &rand2);  
        
        v.x[2*i] = rand1;

        square_sum += rand1*rand1;
    }

    square_sum = sqrt(square_sum);
    
    for (i=0; i<v.n; i++){
        v.x[i] = v.x[i] * r / square_sum;
    }

}

//Adds src1 and src2 and stores the result in dest
void add_point(point dest, point src1, point src2){
    int i;

    assert(src1.n == src2.n);
    assert(dest.n == src1.n);

    for (i=0; i<src1.n; i++){
        dest.x[i] = src1.x[i] + src2.x[i];
    }
}

//Subtracts src1-src2 and stores the result in dest
void sub_point(point dest, point src1, point src2){
    int i;

    assert(src1.n == src2.n);
    assert(dest.n == src1.n);

    for (i=0; i<src1.n; i++){
        dest.x[i] = src1.x[i] - src2.x[i];
    }
}

void ball_walk(point start, int steps, double delta, inclusion_oracle in){
    point query_point = copy_point(start);
    point dir;
    int i;

    dir.n = start.n;
    dir.x = malloc(dir.n*sizeof(double));

    //printf("%lf\n", delta);

    for (i=0; i<steps; i++){
        //print_point(start);
        //printf("\n");

        random_dir(dir, delta);
        
        add_point(query_point, start, dir);

        if (in(query_point)){
            add_point(start, start, dir);
        }
    }
    free(dir.x);

}

//Perturbs start according to a random walk in space after performing a particular
//number of steps in a uniform grid random walk, with steps of size delta
void grid_walk(point start, int steps, double delta, inclusion_oracle in){
    //The current location in the walk
    point query_point = copy_point(start);
    int i,j, k;
    //The directions in which the walk might proceed. This is implemented in a
    //slightly fudgy way, using forced sign-and-magnitude representation for ints,
    //meaning I can use a negative 0. The MSB is used as the sign bit. Since the
    //number of dimensions must be positive, using a signed type to represent it
    //ensures that this bit is available
    int* valid_dirs = malloc(start.n*sizeof(int));
    int msb_mask = 0;
    //The number of directions in which the walk might proceed
    int num_valid_dirs = 0;

    //printf("msb_mask = %d\n", msb_mask);
    msb_mask = ~(((unsigned int) ~msb_mask)>>1);
    //printf("msb_mask = %d\n", msb_mask);

    for (i=0; i<steps; i++){
        //print_point(start);
        //printf("\n");
        if(!(genrand_real1()<WAIT_PROB)){
            num_valid_dirs = 0;
            for(j=0; j<start.n; j++){
                //Can we move in the positive dimension by distance delta?
                query_point.x[j] += delta;
                if (in(query_point)){
                    valid_dirs[num_valid_dirs] = j;
                    num_valid_dirs++;
                }

                //How about the negative?
                query_point.x[j] -= 2*delta;
                if (in(query_point)){
                    valid_dirs[num_valid_dirs] = j^msb_mask;
                    num_valid_dirs++;
                }

                //And now reset the query point
                query_point.x[j] += delta;
            }

            //Pick a random valid direction in which to move
            k = genrand_int32() % num_valid_dirs;

            if(!(valid_dirs[k] & msb_mask)){
                //If the direction of travel's MSB is unset, go that way by a
                //positive amount
                //printf("%d\n", valid_dirs[k]);
                query_point.x[valid_dirs[k]] += delta;
                start.x[valid_dirs[k]] += delta;
            } else {
                //Otherwise, negative
                //printf("-%d\n", valid_dirs[k]^msb_mask);
                query_point.x[valid_dirs[k]^msb_mask] -= delta;
                start.x[valid_dirs[k]^msb_mask] -= delta;
            }
        } else {
            //printf("X\n");
        }
    }

    //For some reason this throws SIGABRTs on Cygwin
    //free(query_point.x);
    //free(valid_dirs);

}

//The square [0,1]^5
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

int main (int argc, char** argv) {
    int i;

    if (argc != 2){
        printf("NOPE NOPE NOPE");
        return 0;
    }

    point p;
    clock_t seed = clock();

    init_genrand(clock());

    p.n = 5;
    p.x = malloc(p.n*sizeof(double));
    for(i = 0; i<1000; i++){
        p.x[0] = 0.5;
        p.x[1] = 0.5;
        p.x[2] = 0.5;
        p.x[3] = 0.5;
        p.x[4] = 0.5;

        printf("%d,", 1<<atoi(argv[1]));
        printf("%d,", seed);
        ball_walk(p, 1<<atoi(argv[1]), 0.002, &square);

        print_point(p);
        printf("\n");
    }

    return 0;
}
