//Delicious Mersenne Twister PRNG
#include"mt.h"
//And ziggurat normal PRNG
#include"ziggurat.h"

#include"sampling.h"

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

/*
typedef struct __point__{
    float* x;
    int n;
} point;

//An inclusion oracle should take a point and return 0 iff the point is not
//inside the shape
typedef int (*inclusion_oracle)(point);
typedef float (*density)(point);
*/


float ZIGGURAT_FN[128];
int   ZIGGURAT_KN[128];
float ZIGGURAT_WN[128];
unsigned long int ZIGGURAT_SEED;




void print_point(point p){
    int i;
    for (i=0; i<p.n-1; i++){
        printf("%f,", p.x[i]);
    }
    printf("%f", p.x[i]);
}

//Returns a continuous, uniormly distributed random-ish number between a and b
//seeding function should be called beforehand
float uniform(float a, float b){
    return genrand_real1()*(b-a) + a;
}

//Creates a freshly allocated copy of the point p. Use with caution and collect
//your garbage
point copy_point(point p){
    point new_point;
    int i;

    new_point.n = p.n;
    new_point.x = malloc(p.n*sizeof(float));
    for(i=0; i<p.n; i++){
        new_point.x[i] = p.x[i];
    }
    return new_point;
}

//Stores two normal gaussian variables in dest1 and dest2
void gaussians(float* dest1, float* dest2){
    float rand1 = uniform(0,1);
    float rand2 = uniform(0,1);

    *dest1 = sqrt(-2*log(rand1)) * cos(2*PI*(rand2));
    *dest2 = sqrt(-2*log(rand1)) * sin(2*PI*(rand2));
}

//Store a uniformly selected random vector of magnitude r in v
//http://mathworld.wolfram.com/HyperspherePointPicking.html
//Probably worth bringing a Ziggurat library in, time permitting.
//http://people.sc.fsu.edu/~jburkardt/c_src/ziggurat/ziggurat.html
void random_dir(point v, float r){
    int i;
    float rand1, rand2;
    float square_sum = 0;
    /*
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
    }*/
    
    for (i=0; i<v.n; i++){
        rand1 = r4_nor(&ZIGGURAT_SEED, ZIGGURAT_KN, ZIGGURAT_FN, ZIGGURAT_WN);
        v.x[i] = rand1;
        //printf("%f\n", rand1);
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

    //assert(src1.n == src2.n);
    //assert(dest.n == src1.n);

    for (i=0; i<src1.n; i++){
        dest.x[i] = src1.x[i] + src2.x[i];
    }
}

//Adds src1+(lambda*src2) and stores the result in dest
void addmul_point(point dest, point src1, float lambda, point src2){
    int i;

    //assert(src1.n == src2.n);
    //assert(dest.n == src1.n);

    for (i=0; i<src1.n; i++){
        dest.x[i] = src1.x[i] + (lambda * src2.x[i]);
    }
}

//Subtracts src1-src2 and stores the result in dest
void sub_point(point dest, point src1, point src2){
    int i;

    //assert(src1.n == src2.n);
    //assert(dest.n == src1.n);

    for (i=0; i<src1.n; i++){
        dest.x[i] = src1.x[i] - src2.x[i];
    }
}

//Returns the squared distance between the given points
float sq_dist(point p1, point p2){
    int i;
    float dist = 0.0;

    //assert(p1.n == p2.n);

    for (i=0; i<p1.n; i++){
        dist += (p1.x[i] - p2.x[i]) * (p1.x[i] - p2.x[i]);
    }

    return dist;
}

//Copies the point from src to dest, assuming dest is already allocated
void move_point(point dest, point src){
    int i;

    //assert(src.n == dest.n);

    for(i=0; i<src.n; i++){
        dest.x[i] = src.x[i];
    }
}

//Returns the distance from the starting point to the boundary of the convex
//shape in the given direction, whose diameter is bounded above by the given,
//correct to error bounds +- sqrt(epsilon), assuming dir is of unit length
//If diameter is negative, the returned distance will be nagative and in the
//opposite direction from dir
float boundary_distance(point start, point dir, float epsilon, float diameter, inclusion_oracle in){
    point lower = copy_point(start);
    point upper;
    point query_point;
    float distance = 0;

    //printf("Finding boundary\n");

    //printf("starting at ");
    //print_point(start);
    //printf("\n");

    upper.n = start.n;
    upper.x = malloc(start.n*sizeof(float));

    query_point.n = start.n;
    query_point.x = malloc(start.n*sizeof(float));

    addmul_point(upper, start, diameter, dir);

    //assert(in(lower));
    //printf("Ending at ");
    //print_point(upper);
    //printf("\n");
    //printf("diameter = %f\n", diameter);
    
    if (in(upper)){
        printf("WARNING: Upper bound was inside shape - select larger diameter\n");
        print_point(upper);
        printf("\n");
        printf("Jumping in direction ");
        print_point(dir);
        printf("\n");
        printf("diameter = %f\n", diameter);
    }
    /*
    while(in(upper)){
        print_point(upper);
        printf("\n");
        printf("Jumping in direction ");
        print_point(dir);
        printf("\n");
        printf("diameter = %f\n", diameter);
        diameter *= 2;
        addmul_point(upper, start, diameter, dir);
    }
    */
    while (sq_dist(lower, upper) > epsilon){
        //printf("dist = %f\n", sq_dist(lower, upper));
        diameter = diameter/2;
        //printf("diameter = %f\n", diameter);
        addmul_point(query_point, lower, diameter, dir);
        if (in(query_point)){
            move_point(lower, query_point);
            distance += diameter;
        } else {
            move_point(upper, query_point);
        }
    }

    free(upper.x);
    free(lower.x);
    free(query_point.x);

    //printf("distance = %f\n", distance);

    //printf("======\n");

    return distance;
}

float max(float a, float b){
    if (a>b) return a;
    return b;
}

float min(float a, float b){
    if (a<b) return a;
    return b;
}

float find_argmax(density f, point dir, point start, float epsilon,
                    float lbound, float ubound){
    point q1,q2;

    q1.n = dir.n;
    q1.x = malloc(dir.n*sizeof(float));

    q2.n = dir.n;
    q2.x = malloc(dir.n*sizeof(float));

    //INVARIANT: The function increases to the right of lbound and decreases to
    //the right of ubound, implying that the maximum is in (lbound,ubound)
    while(ubound-lbound>epsilon){
        addmul_point(q1, start, lbound+(ubound-lbound)/2, dir);
        addmul_point(q2, q1, epsilon, dir);

        if (f(q1)<f(q2)){
            //f increases to the right of lbound+(ubound-lbound)/2
            //This should be the new lower bound
            lbound = lbound+(ubound-lbound)/2;
        } else {
            //f decreases to the right of lbound+(ubound-lbound)/2
            //This should be the new upper bound
            ubound = lbound+(ubound-lbound)/2;
        }
    }

    free(q1.x);
    free(q2.x);

    return (ubound+lbound)/2;
}



float find_argval_right(density f, point dir, point start, float epsilon,
                        float lbound, float ubound){
    float target;

    point p;

    p.n = start.n;
    p.x = malloc(p.n*sizeof(float));

    addmul_point(p, start, lbound, dir);

    target = f(p) * epsilon;

    while(ubound-lbound>epsilon){
        addmul_point(p,start,lbound+(ubound-lbound)/2, dir);
        if (f(p)>target){
            lbound = (lbound+ubound)/2;
        } else {
            ubound = (lbound+ubound)/2;
        }
    }

    free(p.x);

    return (ubound+lbound)/2;
}

float find_argval_left(density f, point dir, point start, float epsilon,
                        float lbound, float ubound){
    float target;

    point p;

    p = copy_point(start);

    p.n = start.n;
    p.x = malloc(p.n*sizeof(float));

    addmul_point(p, start, lbound, dir);

    target = f(p) * epsilon;

    while(ubound-lbound>epsilon){
        addmul_point(p,start,lbound+(ubound-lbound)/2, dir);
        if (f(p)<target){
            lbound = (lbound+ubound)/2;
        } else {
            ubound = (lbound+ubound)/2;
        }
    }

    free(p.x);

    return (ubound+lbound)/2;
}

float rejection_sample(density f, point dir, point start, float lbound,
                            float a, float b, float m){
    point p;
    float f_max;
    float y,r;

    if (a==b) return 0;

    p.n = dir.n;
    p.x = malloc(p.n*sizeof(float));

    addmul_point(p, start, m, dir);

    f_max = f(p);

    do {
        y = uniform(a, b);
        addmul_point(p, start, y, dir);
        r = uniform(0,1);
    } while(f(p) < r*f_max);

    //printf("Sampled %f\n", y);

    return y;
}


void hit_and_run_dist(point start, int steps, float epsilon,
                        float diameter, inclusion_oracle in, density f){
    point dir;
    //Upper and lower bounds on the distance in which we can walk in direction dir
    //Tight within an error of epsilon
    float ubound, lbound;
    //m maximises f(dir*m), a<m and b>m both satisfy f(dir*a) = epsilon*f(dir*m)
    float m,a,b;
    //The size of the step we eventually take
    float step;
    int i;

    dir.n = start.n;
    dir.x = malloc(dir.n*sizeof(float));

    //printf("Start: ");
    //print_point(start);
    //printf("\n");
    //assert(in(start));
    
    
    for (i=0; i<steps; i++){
        //print_point(start);
        //printf("\n");
        random_dir(dir, 1);
        
        //printf("Jumping in direction ");
        //print_point(dir);
        //printf("\n");

        ubound = boundary_distance(start, dir, epsilon,  diameter, in);
        lbound = boundary_distance(start, dir, epsilon, -diameter, in);

        m = find_argmax(f, dir, start, epsilon, lbound, ubound);

        a =  find_argval_left(f, dir, start, epsilon, lbound, m);
        b = find_argval_right(f, dir, start, epsilon, m, ubound);

        //printf("Boundaries are [%f,%f], max at %f\n", a, b, m);

        step = rejection_sample(f, dir, start, lbound, a ,b, m);

        addmul_point(start, start, step, dir);
        //assert(in(start));
    }

    
    //printf("Finish: ");
    //print_point(start);
    //printf("\n");

    //assert(in(start));

}

//Performs a hit and run walk for the given number of steps
void hit_and_run(point start, int steps, float epsilon,
                    float diameter, inclusion_oracle in){
    point dir;
    //Upper and lower bounds on the distance in which we can walk in direction dir
    //Tight within an error of epsilon
    float ubound;
    float lbound;
    float dist;
    int i;

    dir.n = start.n;
    dir.x = malloc(dir.n*sizeof(float));

    assert(in(start));
    
    for (i=0; i<steps; i++){


        random_dir(dir, 1);

        ubound = boundary_distance(start, dir, epsilon,  diameter, in);
        lbound = boundary_distance(start, dir, epsilon, -diameter, in);
        
        //print_point(dir);
        //printf("\n");
        dist = uniform(lbound,ubound);
        //printf("%f from [%f,%f]\n", dist, lbound, ubound);
        addmul_point(start, start, dist, dir);
        if (!in(start)){
            printf("WARNING: Escaped shape\n");
            print_point(dir);
            printf("\n");
            printf("%f from [%f,%f]\n", dist, lbound, ubound);
            print_point(start);
            printf("\n");
            print_radius(start);
        }
        //assert(in(start));
    }
}

void ball_walk(point start, int steps, float delta, inclusion_oracle in){
    point query_point = copy_point(start);
    point dir;
    int i;

    dir.n = start.n;
    dir.x = malloc(dir.n*sizeof(float));

    //printf("%f\n", delta);

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
void grid_walk(point start, int steps, float delta, inclusion_oracle in){
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
        print_point(start);
        printf("\n");
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

//The square [-1,1]^n
int square(point p){
    //assert(p.n == 5);
    int i;
    
    for(i=0; i<2; i++){
        if (p.x[i]<-1 || p.x[i] >1) return 0;
    }
    
    if (p.x[3]<-1 || p.x[3] >19) return 0;
    return 1;
}

float funky_dist(point p){
    float sq_sum = 0;
    int i;

    if (!square(p)) return 0;

    for (i=0; i<p.n; i++){
        sq_sum -= (p.x[i])*(p.x[i]);
    }

    return exp(sq_sum);
    //return 1;
}

void rejection_sample_dist(point p, density f, point lower, point upper, float f_max){
    float r;
    int i;

    do{
        for (i=0; i<p.n; i++){
            p.x[i] = uniform(lower.x[i],upper.x[i]);
        }
        r = uniform(0,1);
    } while (f(p) < r*f_max);
}

//int main (int argc, char** argv) {
//    int i;

    /*

    float lbound;
    float ubound;
    float m;
    point p,dir;

    p.n = 5;
    p.x = malloc(p.n*sizeof(float));

    dir.n = 5;
    dir.x = malloc(dir.n*sizeof(float));

    init_genrand(clock());

    p.x[0] = 0.75;
    p.x[1] = 0.75;
    p.x[2] = 0.75;
    p.x[3] = 0.75;
    p.x[4] = 0.75;

    random_dir(dir,1);

    printf("Heading in direction");
    print_point(dir);
    printf("\n");

    ubound = boundary_distance(p, dir, 0.001,  sqrt(5), &square);
    lbound = boundary_distance(p, dir, 0.001, -sqrt(5), &square);

    printf("Boundaries are [%f,%f]\n", lbound, ubound);

    m = find_argmax(&funky_dist, dir, p, 0.001, lbound, ubound);

    printf("Function is maximised at %f\n", m);
    */
/*
    if (argc != 2){
        printf("NOPE NOPE NOPE");
        return 0;
    }

    point p;

    point l;
    point u;

    clock_t seed = clock();

    init_genrand(clock());

    p.n = 5;
    p.x = malloc(p.n*sizeof(float));

    u.n = 5;
    u.x = malloc(p.n*sizeof(float));

    l.n = 5;
    l.x = malloc(p.n*sizeof(float));

    l.x[0] = 0.0;
    l.x[1] = 0.0;
    l.x[2] = 0.0;
    l.x[3] = 0.0;
    l.x[4] = 0.0;

    u.x[0] = 1.0;
    u.x[1] = 1.0;
    u.x[2] = 1.0;
    u.x[3] = 1.0;
    u.x[4] = 1.0;

    for(i = 0; i<1000; i++){
        p.x[0] = 0.0;
        p.x[1] = 0.0;
        p.x[2] = 0.0;
        p.x[3] = 0.0;
        p.x[4] = 0.0;

        printf("%d,", 1<<atoi(argv[1]));
        printf("%d,", seed);
        hit_and_run_dist(p, 1<<atoi(argv[1]), 0.0001, sqrt(5), &square, &funky_dist);
        //rejection_sample_dist(p, &funky_dist, l, u, 1);
        print_point(p);
        printf("\n");
    }

    return 0;
}
*/