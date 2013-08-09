FLAGS = -Wall -g -O3 --std=c99 -fopenmp
CC = gcc

all:
	make volume

mt: mt.c
	$(CC) $(FLAGS) -o ./mt.o ./mt.c

sampling: sampling.c mt.o
	$(CC) $(FLAGS) -o ./sampling.o ./mt.o ./sampling.c

volume: sampling.o mt.o volume_estimation.c
	$(CC) $(FLAGS) -o ./volume.exe ./sampling.o ./mt.o ./volume_estimation.c
    
clean:
	rm -f ./*.exe
	rm -f ./*.o
