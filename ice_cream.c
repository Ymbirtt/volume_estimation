#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<time.h>

#include"sampling.h"

float ICE_CREAM_HEIGHT;
float ICE_CREAM_FACTOR;


//Returns true iff p is inside the n-ball of radius r
int in_scoop(point p, float r){
    int i;
    float sq_sum = 0;
    r*= r;
    
    for(i=0; i<p.n; i++){
        sq_sum += p.x[i] * p.x[i];
    }
    return (sq_sum<r);
}

void print_radius(point p){
    float cone_radius = (ICE_CREAM_HEIGHT - p.x[0]) * ICE_CREAM_FACTOR;
    
    printf("Cone at distance %f has radius %f\n", p.x[0], cone_radius);
}

int in_cone(point p){
    int i;
    float cone_radius = (ICE_CREAM_HEIGHT - p.x[0]) * ICE_CREAM_FACTOR;
    float sq_sum;
    
    
    cone_radius *= cone_radius;
    
    for(i=1;i<p.n;i++){
        sq_sum += p.x[i]*p.x[i];
    }
    
    return(sq_sum<cone_radius);
}

void init_ice_cream(float h){
    ICE_CREAM_HEIGHT = h;
    ICE_CREAM_FACTOR = sqrt(1-(1/(h*h)))/(h-(1/h));
}

int in_ice_cream(point p){
    if (p.x[0] <= 1/ICE_CREAM_HEIGHT){
        return in_scoop(p,1);
    } else {
        return in_cone(p);
    }
}