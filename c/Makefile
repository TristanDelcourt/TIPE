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

factor_base.o: dir factor_base.c factor_base.h
	$(CC) $(CFLAGS) -c factor_base.c -o $(DIR)/factor_base.o

parse_input.o: dir parse_input.c parse_input.h
	$(CC) $(CFLAGS) -c parse_input.c -o $(DIR)/parse_input.o

factor: dir %.o
	$(CC) $(CFLAGS) factor.c $(DIR)/*.o -o $(DIR)/factor -lgmp

# clean
clean:
	@rm -fr $(DIR)