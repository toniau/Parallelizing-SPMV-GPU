#/*
# *********************************************
# *  314 Principles of Programming Languages  *
# *  Fall 2016                                *
# *********************************************
# */

CCFLAGS = -Wall -fopenmp -std=c11

all: verify

verify: verify.c mmio.c mmio.h utils.c utils.h 
	gcc $(CCFLAGS) verify.c mmio.c utils.c utils.h -o verify


clean:
	rm -f  verify




