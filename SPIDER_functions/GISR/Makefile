CC=gcc
CFLAGS=-lm -O3

clean:
	rm add_column remove_column

all: add_column remove_column

add_column: src/add_column.c
	$(CC) $^ -o $@ $(CFLAGS)

remove_column: src/remove_column.c
	$(CC) $^ -o $@ $(CFLAGS)
