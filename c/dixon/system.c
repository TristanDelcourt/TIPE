#include "system.h"
#include "vector.h"
#include "list_matrix_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

void swap_lines_horz(system_t s, int i, int j){
    int* temp = s->m[i];
    s->m[i] = s->m[j];
    s->m[j] = temp;
}

void swap_lines_vert(system_t s, int i, int j){
    int temp = s->perm[i];
    s->perm[i] = s->perm[j];
    s->perm[j] = temp;

    for(int k = 0; k<s->n1; k++){
        int temp = s->m[k][i];
        s->m[k][i] = s->m[k][j];
        s->m[k][j] = temp;
    }
}

int find_index(system_t s, int from, int look){
    for(int i = from; i < s->n1; i++){
        if(s->m[i][look]){
            return i;
        }
    }
    return -1;
}

system_t transpose(int** v, int n1, int n2){
    system_t s = malloc(sizeof(system_s));

    s->m = malloc(n2*sizeof(int*));
    for(int i = 0; i<n2; i++){
        s->m[i] = malloc(n1*sizeof(int));
        for(int j = 0; j<n1; j++){
            s->m[i][j] = v[j][i];
        }
    }

    s->n1 = n2;
    s->n2 = n1;
    return s;
}

void triangulate(system_t s){
    s->perm = malloc(s->n2*sizeof(int));
    for(int i = 0; i<s->n2; i++){
        s->perm[i] = i;
    }

    int i = 0;
    int j = 0;
    while(i<s->n1 && j<s->n2){
        int k = find_index(s, i, j);
        if(k != -1){
            if(i != j){
                swap_lines_vert(s, i, j);
            }

            swap_lines_horz(s, i, k);

            for(int l = i + 1; l < s->n1; l++){
                if(s->m[l][i] == 1){
                    sub_vect(s->m, l, i, s->n2);
                    mod_vect(s->m[l], 2, s->n2);
                }
            }
            i++;
            j = i;
        }
        else{
            j++;
        }
    }
}

void get_arbitary(system_t triangulated){
    for(int i = triangulated->n1-1; i>=0; i--){
        int j = 0;
        while(j < triangulated->n2 && !triangulated->m[i][j]){
            j++;
        }
        if(j<triangulated->n2){
            triangulated->arb = j+1;
            return;
        }
    }

    fprintf(stderr, "ERROR: All vectors are zero in system\n");
    exit(1);
}

void init_sol(system_t s){
    s->sol = malloc(s->n2*sizeof(int));
    for(int i = s->arb; i<s->n2; i++){
        s->sol[i] = 0;
    }
}

void iter_sol(system_t s){
    int i = s->arb;
    while(i<s->n2 && (s->sol[i] == 1)){
        s->sol[i] = 0;
        i++;   
    }
    if(i >= s->n2){
        s->done = true;
        return;
    }
    s->sol[i] = 1;
}

system_t init_gauss(int** v, int n1, int n2){
    //printf("Initial vectors\n");
    //print_ll(v, n1, n2);
    
    system_t s = transpose(v, n1, n2);
    s->done = false;

    //printf("Transposed\n");
    //print_ll(s->m, s->n1, s->n2);
    
    for(int i = 0; i<s->n1; i++){
        mod_vect(s->m[i], 2, s->n2);
    }

    //printf("Modded\n");
    //print_ll(s->m, s->n1, s->n2);

    triangulate(s);
    
    printf("Triangulated\n");
    print_ll(s->m, s->n1, s->n2);

    get_arbitary(s);
    init_sol(s);

    return s;
}

void gaussian_step(system_t s){
    iter_sol(s);

    for(int i = s->n1-1; i>=0; i--){
        int j = 0;
        while(j < s->n2 && !s->m[i][j]){
            j++;
        }

        if(j<s->n2){
            s->sol[j] = 0;

            for(int k = s->n2-1; k>j; k--){
                s->sol[j] -= s->m[i][k] * s->sol[k];
            }
            s->sol[j] = abs(s->sol[j]) % 2;
        }
    }
}

void free_system(system_t s){
    for(int i = 0; i<s->n1; i++){
        free(s->m[i]);
    }
    free(s->m);
    free(s->sol);
    free(s->perm);
    free(s);
}