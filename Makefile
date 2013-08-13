FLAGS = -Wall -g  --std=c99 -pg
CC = gcc

all:
	make volume

mt: mt.c
	$(CC) $(FLAGS) -o ./mt.o ./mt.c

sampling: sampling.c mt.o ziggurat.o
	$(CC) $(FLAGS) -o ./sampling.o ./mt.o ./ziggurat.o ./sampling.c

ziggurat: ziggurat.c
	$(CC) $(FLAGS) -o ./ziggurat.o ./ziggurat.c
    
volume: ziggurat.o sampling.o mt.o lovasz-vempala.c
	$(CC) $(FLAGS) -o ./volume.exe ./sampling.o ./ziggurat.o ./mt.o ./lovasz-vempala.c
    
clean:
	rm -f ./*.exe
	rm -f ./*.o
