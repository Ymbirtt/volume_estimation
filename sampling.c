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

//Returns a freshly allocated point in space after performing a particular
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
        grid_walk(p, 1<<atoi(argv[1]), 0.001, &square);

        print_point(p);
        printf("\n");
    }

    return 0;
}
