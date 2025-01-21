CC := gcc
CFLAGS := -Wall -Wextra -g
DIR := ./bin

all: factor
dir:
	@mkdir -p $(DIR)

list_matrix_utils.o: dir list_matrix_utils.c list_matrix_utils.h
	$(CC) $(CFLAGS) -c list_matrix_utils.c -o $(DIR)/list_matrix_utils.o

vector.o: dir vector.c vector.h
	$(CC) $(CFLAGS) -c vector.c -o $(DIR)/vector.o

system.o: dir system.c system.h
	$(CC) $(CFLAGS) -c system.c -o $(DIR)/system.o

factorbase.o: dir factorbase.c factorbase.h
	$(CC) $(CFLAGS) -c factorbase.c -o $(DIR)/factorbase.o

parse_input.o: dir parse_input.c parse_input.h
	$(CC) $(CFLAGS) -c parse_input.c -o $(DIR)/parse_input.o

# Algorithms

qsieve.o:
	$(CC) $(CFLAGS) -c ./qsieve/qsieve.c -o $(DIR)/qsieve.o

dixon.o:
	$(CC) $(CFLAGS) -c ./dixon/dixon.c -o $(DIR)/dixon.o

factor: dir list_matrix_utils.o vector.o system.o factorbase.o parse_input.o qsieve.o dixon.o
	$(CC) $(CFLAGS) main.c \
		$(DIR)/list_matrix_utils.o \
		$(DIR)/vector.o \
		$(DIR)/system.o \
		$(DIR)/factorbase.o \
		$(DIR)/parse_input.o \
		$(DIR)/qsieve.o \
		$(DIR)/dixon.o \
		 -o $(DIR)/factor -lgmp

# clean
clean:
	@rm -fr $(DIR)