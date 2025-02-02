CC := gcc
CFLAGS := -Wall -Wextra -g
CLINKS := -lgmp -lm
DIR := ./bin

all: factor
dir:
	@mkdir -p $(DIR)

list_matrix_utils.o: dir list_matrix_utils.c list_matrix_utils.h
	$(CC) $(CFLAGS) -c list_matrix_utils.c -o $(DIR)/list_matrix_utils.o $(CLINKS)

vector.o: dir vector.c vector.h
	$(CC) $(CFLAGS) -c vector.c -o $(DIR)/vector.o $(CLINKS)

system.o: dir system.c system.h
	$(CC) $(CFLAGS) -c system.c -o $(DIR)/system.o $(CLINKS)

factorbase.o: dir factorbase.c factorbase.h
	$(CC) $(CFLAGS) -c factorbase.c -o $(DIR)/factorbase.o $(CLINKS)

tonellishanks.o: dir tonellishanks.c tonellishanks.h
	$(CC) $(CFLAGS) -c tonellishanks.c -o $(DIR)/tonellishanks.o $(CLINKS)

parse_input.o: dir parse_input.c parse_input.h
	$(CC) $(CFLAGS) -c parse_input.c -o $(DIR)/parse_input.o $(CLINKS)

polynomial.o: dir ./mpqs/polynomial.c ./mpqs/polynomial.h
	$(CC) $(CFLAGS) -c ./mpqs/polynomial.c -o $(DIR)/polynomial.o $(CLINKS)

common_mpqs.o: dir ./mpqs/common_mpqs.c ./mpqs/common_mpqs.h
	$(CC) $(CFLAGS) -c ./mpqs/common_mpqs.c -o $(DIR)/common_mpqs.o $(CLINKS)

parallel_mpqs.o: dir ./mpqs/parallel_mpqs.c ./mpqs/parallel_mpqs.h
	$(CC) $(CFLAGS) -c ./mpqs/parallel_mpqs.c -o $(DIR)/parallel_mpqs.o $(CLINKS)

# Algorithms

qsieve.o: ./qsieve/qsieve.c ./qsieve/qsieve.h
	$(CC) $(CFLAGS) -c ./qsieve/qsieve.c -o $(DIR)/qsieve.o $(CLINKS)

dixon.o: ./dixon/dixon.c ./dixon/dixon.h
	$(CC) $(CFLAGS) -c ./dixon/dixon.c -o $(DIR)/dixon.o $(CLINKS)

mpqs.o: ./mpqs/mpqs.c ./mpqs/mpqs.h
	$(CC) $(CFLAGS) -c ./mpqs/mpqs.c -o $(DIR)/mpqs.o $(CLINKS)

factor: dir list_matrix_utils.o vector.o system.o factorbase.o tonellishanks.o polynomial.o parse_input.o qsieve.o dixon.o mpqs.o common_mpqs.o parallel_mpqs.o
	$(CC) $(CFLAGS) main.c \
		$(DIR)/list_matrix_utils.o \
		$(DIR)/vector.o \
		$(DIR)/system.o \
		$(DIR)/factorbase.o \
		$(DIR)/tonellishanks.o \
		$(DIR)/polynomial.o \
		$(DIR)/parse_input.o \
		$(DIR)/qsieve.o \
		$(DIR)/dixon.o \
		$(DIR)/mpqs.o \
		$(DIR)/common_mpqs.o \
		$(DIR)/parallel_mpqs.o \
		-o $(DIR)/factor $(CLINKS)

# clean
clean:
	@rm -fr $(DIR)