FLAGS = -Wall -g -O3 --std=c99 -pg
CC = gcc

all:
	make volume

mt: mt.c
	$(CC) $(FLAGS) -o ./mt.o ./mt.c

sampling: sampling.c mt.o ziggurat.o
	$(CC) $(FLAGS) -o ./sampling.o ./mt.o ./ziggurat.o ./sampling.c

ziggurat: ziggurat.c
	$(CC) $(FLAGS) -o ./ziggurat.o ./ziggurat.c

ice_cream: ice_cream.c sampling.c
	$(CC) $(FLAGS) -o ./ice_cream.o ./sampling.o ./ice_cream.c

volume: ziggurat.o sampling.o mt.o ice_cream.o volume_estimation.c
	$(CC) $(FLAGS) -o ./volume.exe ./sampling.o ./ziggurat.o ./mt.o ./ice_cream.o ./volume_estimation.c
    
clean:
	rm -f ./*.exe
	rm -f ./*.o
