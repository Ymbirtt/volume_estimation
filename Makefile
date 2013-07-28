FLAGS = -Wall -g -O3 --std=c99
CC = gcc

all:
	make sampling

mt: mt.c
	$(CC) $(FLAGS) -o ./mt.o ./mt.c

sampling: sampling.c mt.o
	$(CC) $(FLAGS) -o ./sampling.exe ./mt.o ./sampling.c

clean:
	rm -f ./*.exe
	rm -f ./*.o
