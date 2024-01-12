CC = gcc
CFLAGS  = -Wall -fopenmp -g -lm

default: scene

scene:
	$(CC) $(CFLAGS) -o scene lensing.c linalg.c object.c scene.c tracing.c

clean: 
	$(RM) scene *.o *~