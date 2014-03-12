/* Generate binomial random fields for Monte Carlo simulations using the
 * Mersenne Twister algorithm. 
 *
 * Written in 2006 for the benefit of Tuorla Observatory Cosmological
 * Simulations research group by Pauli Pihajoki.
 * 
 * Copyright (c) 2006 Pauli Pihajoki. All rights reserved. 
 */

#include <stdlib.h>
#include <stdio.h>

/* Use Matsumoto's Mersenne twister algorithm. */
#include "mersenne/mt19937ar.h"

/* Print usage. */
void usage(void) {
	printf("fieldgen: arguments: outputfile, seed, box size, "
	       "number of points\n");
}

/* Controls how big blocks are generated and written in turns. */
#define BLOCK_SIZE	1000


int main(int argc, char **argv) {
	FILE *OUT;
	int box_size;
	unsigned long seed, points;
 
	if (argc != 5) {
		usage();
		exit(EXIT_FAILURE);
	}

	/* Open output file. */
	OUT = fopen(argv[1], "w");
	if (OUT == NULL) {
		perror("fopen");
		exit(EXIT_FAILURE);
	}

	/* Seed and number of points should be unsigned long. */
	seed = atol(argv[2]);
	points = atol(argv[4]);
	if (points <= 0) {
		fprintf(stderr, "Invalid number of points: %ld\n", points);
		exit(EXIT_FAILURE);
	}

	/* Box size should be an integer. */
	box_size = atoi(argv[3]);
	if (box_size <= 0) {
		fprintf(stderr, "Invalid box size: %d\n", box_size);
		exit(EXIT_FAILURE);
	}

	/* Print HaloAnalyzer-compatible field preamble, constisting on the
	 * number of points in the field and the lower and upper bounds. */
	fprintf(OUT, "%ld\n", points);
	fprintf(OUT, "0 0 0\n%d %d %d\n", box_size, box_size, box_size);

	/* Initialize the mersenne twister algorithm. */
	init_genrand(seed);

	/* Generate a block of random numbers, and print them to output file.
	 * Rinse, lather and repeat. */
	{
		unsigned long i,j, old=0; 
		unsigned long chg=points/100;
		double p[BLOCK_SIZE][3];
		for (i=0; i < points; i+=BLOCK_SIZE) {
			/* Print some progress info. */
			if (i-old > chg) {
				old = i;
				printf("%ld%% done\n", i*100/points);
			}

			/* Generate a block of mt random numbers from [0,1]
			 * and multiply to [0,box_size]. Then write them. */
			for (j=0; j < BLOCK_SIZE; j++) {
				p[j][0] = box_size * genrand_real1();
				p[j][1] = box_size * genrand_real1();
				p[j][2] = box_size * genrand_real1();
			}
			for (j=0; j < BLOCK_SIZE; j++) {
				fprintf(OUT, "%16.16lg %16.16lg %16.16lg\n",
					p[j][0], p[j][1], p[j][2]);
			}
		}
	}

	/* Cleanup. */
	fclose(OUT);

	exit(EXIT_SUCCESS);
}
