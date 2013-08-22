#include"sampling.h"
#include"mt.h"
#include"glpk.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<time.h>
#include <float.h>
#include <limits.h>

//Functions for dealing with the convex hull of a set of points. Requires glpk
//compile with gcc (flags) (file) -lglpk

glp_prob* LP;
glp_smcp LP_PARAMS;

int query_hull(point query){
    int i;

    for (i=0; i<query.n; i++){
        glp_set_row_bnds(LP, i+1, GLP_FX, query.x[i], query.x[i]);
    }
    i = glp_simplex(LP, &LP_PARAMS);
    
    if (i!=0){
        return 0;
    } else {
        return 1;
    }
}

void prepare_point(double* dest, float* src, int n){
    int i;
    
    for(i=0; i<n; i++){
        dest[i+1] = src[i];
    }
    dest[i+1] = 1.0;
}

//Assigns the private LP problem based on the given set of points
void init_hull(point* ps, int n){
    int i;
    int* inds;
    int* range;
    double* ones;
    double* temp;

    range = malloc(n*sizeof(int));
    inds = malloc(ps[0].n*sizeof(int)) - sizeof(int);
    ones = malloc(n*sizeof(double)) - sizeof(double);
    temp = malloc((ps[0].n+1)*sizeof(double)) - sizeof(double);
    
    
    for (i=1; i<=n; i++){
        range[i] = i;
        ones[i]  = 1;
    }
    for (i=1; i<=ps[0].n+1; i++){
        inds[i]  = i;
    }

    LP_PARAMS.presolve = GLP_ON;

    LP_PARAMS.meth = GLP_PRIMAL;
    LP_PARAMS.r_test = GLP_RT_HAR;
    LP_PARAMS.msg_lev = GLP_MSG_OFF;
    LP_PARAMS.pricing = GLP_PT_PSE;
    LP_PARAMS.tol_bnd = 1e-7;
    LP_PARAMS.tol_dj = 1e-7;
    LP_PARAMS.tol_piv = 1e-10;
    LP_PARAMS.obj_ll = -DBL_MAX;
    LP_PARAMS.obj_ul = +DBL_MAX;
    LP_PARAMS.it_lim = INT_MAX;
    LP_PARAMS.tm_lim = INT_MAX;
    LP_PARAMS.out_frq = 500;
    LP_PARAMS.out_dly = 0;


    //glp_delete_prob(LP);
    LP = glp_create_prob();

    //Set objective direction - max or min
    glp_set_obj_dir(LP, GLP_MAX);

    //Set the number of free variables - there should be one for each point on
    //this convex hull. We aim to find coefficients such that a queried point is
    //the convex combination of these coefficients and the given set of points
    glp_add_cols(LP, n);

    //Set the number of constraints - there should be one for each dimension,
    //plus an extra to ensure that the coefficients all sum to 1
    glp_add_rows(LP,ps[0].n+1);

    
    print_problem();
    
    glp_set_row_bnds(LP, ps[0].n+1, GLP_FX, 1,1);
    
    //Each column contains the coordinates of the points on the convex hull
    //We will be taking convex combinations of these points
            
    //The first point should always be added to the hull
    prepare_point(temp, ps[0].x, ps[0].n);
    
    //Each coordinate takes its place on the coefficient matrix
    glp_set_mat_col(LP, 0+1, ps[0].n+1, inds, temp);

    //The target coefficients should be between 0 and 1
    glp_set_col_bnds(LP, 0+1, GLP_DB, 0,1);
    printf("===========\n");
    for(i=1; i<n; i++){
        assert(ps[i].n == ps[0].n);
        printf("===========\n");
        //print_problem();
        //If the point we're trying to add is already within the convex hull of
        //the previous points, we've no need to add it and can move on to the next 
        if(query_hull(ps[i])){
            printf("Point is already within convex hull, no need to add\n");
        } else {
            printf("New point on hull, adding\n");
        
            prepare_point(temp, ps[i].x, ps[i].n);
            
            //Each coordinate takes its place on the coefficient matrix
            glp_set_mat_col(LP, i+1, ps[0].n+1, inds, temp);

            //The target coefficients should be between 0 and 1
            glp_set_col_bnds(LP, i+1, GLP_DB, 0,1);
        }
        printf("===========\n");
    }
    printf("===========\n");
    print_problem();
    //query_hull(ps[i]);
    printf("===========\n");


    //The first row is the sum of all the coefficients. We constrain this such
    //that it sums to 1
    glp_set_mat_row(LP, ps[0].n+1, n, range, ones);
    glp_set_row_bnds(LP, ps[0].n+1, GLP_FX, 1,1);
    //The other row sums will equal the relevant coordiante of the queried point
    
}

void print_problem(){
    int num_cols = glp_get_num_cols(LP);
    int num_rows = glp_get_num_rows(LP);
    int i,j;
    int* range = malloc(glp_get_num_cols(LP)*sizeof(int)) - sizeof(int);
    double* vec = malloc(glp_get_num_cols(LP)*sizeof(double)) - sizeof(double);
    
    for(i=1; i<=num_cols; i++){
        range[i] = i;
        vec[i] = 0;
    }
    
    for(j=1; j<=num_rows; j++){
        glp_get_mat_row(LP, j, range, vec);
        for(i=1; i<=num_cols; i++){
            printf("% 1.2lf ", vec[i]);
            vec[i] = 0;
        }
        printf("\n");
    }
    
    printf("\n");
    
}

int main (void){
    point ps[8];
    point query;
    int i;

    for(i=0;i<8;i++){
        ps[i].n = 3;
        ps[i].x = malloc(3*sizeof(float));
    }

    query.n = 3;
    query.x = malloc(3*sizeof(float));
    query.x[0] = -0.9, query.x[1] = -1.0000001, query.x[2] = -0.9;

    
    ps[0].x[0] = -1, ps[0].x[1] = -1, ps[0].x[2] = -1;
    ps[1].x[0] = -1, ps[1].x[1] = -1, ps[1].x[2] = +1;
    ps[2].x[0] = -1, ps[2].x[1] = +1, ps[2].x[2] = -1;
    ps[3].x[0] = -1, ps[3].x[1] = +1, ps[3].x[2] = +1;
    ps[4].x[0] = +1, ps[4].x[1] = -1, ps[4].x[2] = -1;
    ps[5].x[0] = +1, ps[5].x[1] = -1, ps[5].x[2] = +1;
    ps[6].x[0] = +1, ps[6].x[1] = +1, ps[6].x[2] = -1;
    ps[7].x[0] = +1, ps[7].x[1] = +1, ps[7].x[2] = +1;

    init_hull(ps, 8);

    print_problem();
    
    if(query_hull(query)){
        printf("Queried point was within hull\n");
    } else {
        printf("Queried point was not within hull\n");
    }

}