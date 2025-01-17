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

void test_dixon_b_values(int start, int end, int step, int digits){
    char* file = "/home/idle/Work/TIPE/c/dixon/bin/dixon";
    char* number = malloc(50*sizeof(char));
    char* b_value = malloc(50*sizeof(char));
    number[0] = '\0';

    char* file_name_pt2 = malloc(50*sizeof(char));
    sprintf(file_name_pt2, "%d", digits);
    char* argv[3] = {"./products/", file_name_pt2, "_digit_products.txt"};
    
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    int status;

    int n_tests = (end-start)/step + 1;
    double* averages = malloc(n_tests*sizeof(double));
    float* success_rates = malloc(n_tests*sizeof(double));
    int i = 0;
    for(int b = start; b<=end; b+=step){
        printf("=====================\nB = %d\n=====================\n", b);

        double average = 0;
        float success_rate = 0;
        FILE* f = open_file(argv, 3, 150, "r");
        if (f == NULL) {
            fprintf(stderr, "Failed to open file\n");
            return;
        }

        int n_lines = 0;
        while ((read = getline(&line, &len, f)) != -1) {
            n_lines++;
            line[read-1] = '\0';
            number[0] = '\0';
            sprintf(b_value, "%d", b);
            strcat(number, line);

            printf("n = %s", line);
            fflush(stdout);
            pid_t pid = fork();
            if(pid == 0){
                int error = execl(file, "dixon", "1", number, b_value, NULL);
                if(error == -1){
                    fprintf(stderr, "ERROR: Could not execute command\n: '%s %s %s'\n", file, number, b_value);
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
            average += elapsed;
            if(WEXITSTATUS(status) == 0){
                success_rate += 1.0;
            }
            printf(" / %fs / %d\n", elapsed, WEXITSTATUS(status));
        }
        fclose(f);

        averages[i] = average/n_lines;
        success_rates[i] =  success_rate/n_lines;
        i++;
    }

    printf("=====================\n=====================\n");

    int b = start;
    for(int i = 0; i<n_tests; i++ ){
        printf("Avg: b=%d -> %fs / %.1f%%\n", b, averages[i], success_rates[i]*100);
        b+=step;
    }
    
    free(averages);
    free(file_name_pt2);
    free(number);
    free(b_value);
    free(line);
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

    if(argc>=4){
        bool recompute = atoi(argv[2]);
        if(recompute){
            int nb_tests = atoi(argv[3]);
            for(int i = 2; i<9; i++){
                create_products(i, nb_tests);
            }
        }
    }
    /*
    6 digit products
    Avg: b=250 -> 1.409920s / 100.0%
    Avg: b=260 -> 1.687535s / 100.0%
    Avg: b=270 -> 1.706647s / 100.0%
    Avg: b=280 -> 1.876258s / 100.0%
    */

    test_dixon_b_values(600, 600, 10, 6);

    return 0;
}