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
    input->algorithm = DEFAULT;
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

input_t* parse_input(int argc, char** argv){
    input_t* input = init_input();
    for(int i = 1; i<argc-1; i++){
        if(strcmp(argv[i], "-b") == 0){
            if(i+1<argc){
                if(valid_int(argv[i+1])) input->bound = atoi(argv[i+1]);
                else return NULL;}
            else return NULL;
        }

        else if(strcmp(argv[i], "-s") == 0){
            if(i+1<argc){
                if(valid_int(argv[i+1])) input->sieving_interval = atoi(argv[i+1]);
                else return NULL;}
            else return NULL;
        }
        
        else if(strcmp(argv[i], "-e") == 0){
            if(i+1<argc){
                if(valid_int(argv[i+1])) input->extra = atoi(argv[i+1]);
                else return NULL;}
            else return NULL;
        }

        else if(strcmp(argv[i], "-n") == 0){
            if(i+1<argc){
                if(valid_int(argv[i+1])) mpz_set_str(input->N, argv[i+1], 10);
                else return NULL;}
            else return NULL;
        }

        else if(strcmp(argv[i], "-o") == 0){
            if(i+1<argc) input->output_file = argv[i+1];
            else return NULL;
        }

        else if(strcmp(argv[i], "-t") == 0){
            if(i+1<argc) {
                if(strcmp(argv[i+1], "dixon")) input->algorithm = DIXON;
                else if(strcmp(argv[i+1], "qsieve")) input->algorithm = QSIEVE;
                else return NULL;}
            else return NULL;
        }

        else if(strcmp(argv[i], "-q") == 0){
            input->quiet = true;
        }

        else return NULL;
    }

    return input;
}