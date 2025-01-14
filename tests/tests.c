#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>

FILE* open_file(char** filename_parts, int n, int max_size, char* type){
    char* file_name = malloc(max_size*sizeof(char));
    file_name[0] = '\0';
    for(int i = 0; i<n; i++){
        strcat(file_name, filename_parts[i]);
    }
    FILE* f = fopen(file_name, type);
    free(file_name);
    return f;
}

void test_dixon_b_values(int start, int end, int step){
    char* file = "/home/idle/Work/TIPE/c/dixon/bin/dixon";
    char* arg1 = malloc(50*sizeof(char));
    arg1[0] = '\0';

    char* file_name_pt2 = malloc(50*sizeof(char));
    file_name_pt2[0] = '6';
    file_name_pt2[1] = '\0';
    char* argv[3] = {"./products/", file_name_pt2, "_digit_products.txt"};
    FILE* f = open_file(argv, 3, 150, "r");
    if (f == NULL) {
        fprintf(stderr, "Failed to open file\n");
        return;
    }

    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    int status = 0;
    while ((read = getline(&line, &len, f)) != -1) {
        line[read-1] = '\0';
        arg1[0] = '\0';
        strcat(arg1, line);

        pid_t pid = fork();
        if(pid == 0){
            int error = execl(file, "dixon", "1", arg1, NULL);
            if(error == -1){
                fprintf(stderr, "ERROR: Could not execute command\n: %s %s\n", file, arg1);
                exit(EXIT_FAILURE);
            }
        }
        struct timeval begin, end;
        gettimeofday(&begin, 0);
        pid = wait(&status);
        gettimeofday(&end, 0);
        long seconds = end.tv_sec - begin.tv_sec;
        long microseconds = end.tv_usec - begin.tv_usec;
        double elapsed = seconds + microseconds*1e-6;
        printf("n = %s / Time measured: %f seconds. / %d\n", line, elapsed, WEXITSTATUS(status));
    }
}

void create_products(int ndigit, int nb_tests){
    char* file_name_pt2 = malloc(50*sizeof(char));
    sprintf(file_name_pt2, "%d", ndigit);
    char* argv[3] = {"./primes/", file_name_pt2, "_digit_primes.txt"};
    FILE* f = open_file(argv, 3, 150, "r");

    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    uint32_t p;
    uint32_t allocsize = 50;
    uint32_t* primes = malloc(allocsize*sizeof(uint32_t));
    if (f == NULL) {
        fprintf(stderr, "Failed to open file\n");
        return;
    }
    int number_of_primes = 0;
    while ((read = getline(&line, &len, f)) != -1) {
        sscanf(line, "%"SCNu32"\n", &primes[number_of_primes]);
        number_of_primes++;
        if(number_of_primes >= allocsize){
            allocsize += 50;
            primes = realloc(primes, allocsize*sizeof(uint32_t));
        }
    }
    primes = realloc(primes, number_of_primes*sizeof(uint32_t));

    fclose(f);
    if(line) free(line);

    char* argv2[3] = {"./products/", file_name_pt2, "_digit_products.txt"};
    FILE* pf = open_file(argv2, 3, 150, "w");
    if (pf == NULL) {
        fprintf(stderr, "Failed to open file\n");
        return;
    }

    uintmax_t prod;
    for(int i = 0; i<nb_tests; i++){
        int k = rand()%number_of_primes;
        int l = rand()%number_of_primes;
        prod = (uintmax_t) primes[k] * primes[l];
        fprintf(pf, "%"PRIuMAX"\n", prod);
    }
    fclose(pf);
    free(primes);
    free(file_name_pt2);
}

void parse_primelist(char* source){
    FILE* f = fopen("someprimes.txt", "r");
    if (f == NULL) {
        fprintf(stderr, "Failed to open file\n");
        exit(1);
    }
    char* line = NULL;
    size_t len = 0;
    ssize_t read;

    char* file_name_pt2 = malloc(50*sizeof(char));
    char* argv[3] = {"./primes/", "2", "_digit_primes.txt"};
    FILE* pf = open_file(argv, 3, 150, "r");
    if (pf == NULL) {
        fprintf(stderr, "Failed to open file\n");
        exit(1);
    }

    int current_len = 1;
    while ((read = getline(&line, &len, f)) != -1) {
        int size = (int)read -1;
        if(size!=current_len){
            fclose(pf);
            current_len = (int)size;
            sprintf(file_name_pt2, "%d", size);
            char* argv[3] = {"./primes/", file_name_pt2, "_digit_primes.txt"};
            FILE* pf = open_file(argv, 3, 150, "r");
        }

        fprintf(pf, "%s", line);
    }

    fclose(f);
    fclose(pf);
    if(line) free(line);
    free(file_name_pt2);
}

int main(int argc, char** argv){
    srand(time(NULL));
    if(argc>=2){
        bool recompute = atoi(argv[1]);
        if(recompute) parse_primelist("someprimes.txt");
    }

    if(argc>=2){
        bool recompute = atoi(argv[2]);
        if(recompute){
            int nb_tests = atoi(argv[3]);
            for(int i = 2; i<9; i++){
                create_products(i, nb_tests);
            }
        }
    }

    test_dixon_b_values(0, 0, 0);

    return 0;
}