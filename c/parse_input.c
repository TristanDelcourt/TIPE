#include "parse_input.h"
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <stdbool.h>

input_t* init_input(void){
    input_t* input = malloc(sizeof(input_t));
    input->bound = -1;
    input->output_file = NULL;
    input->sieving_interval = -1;
    input->extra = -1;
    input->quiet = false;
    input->algorithm = QSIEVE;
    mpz_init_set_ui(input->N, 0);
    return input;
}

bool valid_int(char* str){
    int i = 0;
    char c = str[i];
    while(c != '\0'){
        if(c<48 || c>57) return false;
        c = str[++i];
    }

    return true;
}

void free_input(input_t* input){
    if(input->output_file) free(input->output_file);
    mpz_clear(input->N);
    free(input);
}

input_t* parse_input(int argc, char** argv){
    input_t* input = init_input();

    int i = 1;
    while(i<argc){
        if(strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--bound") == 0){
            i++;
            if(i<argc){
                if(valid_int(argv[i])) input->bound = atoi(argv[i]);
                else return NULL;}
            else return NULL;
        }

        else if(strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--sieving_interval") == 0){
            i++;
            if(i<argc){
                if(valid_int(argv[i])) input->sieving_interval = atoi(argv[i]);
                else return NULL;}
            else return NULL;
        }
        
        else if(strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--extra") == 0){
            i++;
            if(i<argc){
                if(valid_int(argv[i])) input->extra = atoi(argv[i]);
                else return NULL;}
            else return NULL;
        }

        else if(strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--number") == 0){
            i++;
            if(i<argc){
                if(valid_int(argv[i])) mpz_set_str(input->N, argv[i], 10);
                else return NULL;}
            else return NULL;
        }

        else if(strcmp(argv[i], "-o") == 0){
            i++;
            if(i<argc) input->output_file = argv[i];
            else return NULL;
 }

        else if(strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--type") == 0){
            i++;
            if(i<argc) {
                if(strcmp(argv[i], "dixon") == 0) input->algorithm = DIXON;
                else if(strcmp(argv[i], "qsieve") == 0) input->algorithm = QSIEVE;
                else if(strcmp(argv[i], "mpqs") == 0) input->algorithm = MPQS;
                else if(strcmp(argv[i], "pmpqs") == 0) input->algorithm = PMPQS;
                else return NULL;}
            else return NULL;
        }

        else if(strcmp(argv[i], "-q") == 0 ||
                strcmp(argv[i], "-stfu") == 0 /*easter egg*/ ||
                strcmp(argv[i], "--quiet") == 0){
            input->quiet = true;
        }

        else return NULL;

        i++;
    }

    return input;
}