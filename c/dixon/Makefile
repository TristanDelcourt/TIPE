CC := gcc
CFLAGS := -Wall -Wextra -g
DIR := ./bin

all: dixon
dir:
	@mkdir -p $(DIR)

list_matrix_utils.o: dir list_matrix_utils.c list_matrix_utils.h
	$(CC) $(CFLAGS) -c list_matrix_utils.c -o $(DIR)/list_matrix_utils.o

vector.o: dir vector.c vector.h
	$(CC) $(CFLAGS) -c vector.c -o $(DIR)/vector.o

system.o: dir system.c system.h
	$(CC) $(CFLAGS) -c system.c -o $(DIR)/system.o

dixon: dir dixon.c system.o vector.o list_matrix_utils.o
	$(CC) $(CFLAGS) dixon.c $(DIR)/list_matrix_utils.o $(DIR)/vector.o $(DIR)/system.o -o $(DIR)/dixon -lgmp


# clean
clean:
	@rm -fr $(DIR)