/* Calculate two-point correlation function for parallel and transverse
 * components. 
 *
 * Usage: compcorr OUTFILE DATAFILE CMPFILE method
 *
 *
 * File CMPFILE should have the following format:
 *
 * number of points (N)
 * data low bound (3 components, separated by whitespace)
 * data high bound (3 components, separated by whitespace)
 * point 1 (3 components, separated by whitespace)
 * .
 * .
 * .
 * point N (3 components, separated by whitespace)
 *
 *
 * File DATAFILE should be otherwise similar to CMPFILE, with the following
 * added to the very beginning of the file:
 *
 * number of parallel scale bins (N_p)
 * bin 1 (2 components, low and high bounds, separated by whitespace)
 * .
 * .
 * bin N_p 
 * number of transverse scale bins (N_t)
 * bin 1 (2 components, low and high bounds, separated by whitespace)
 * .
 * .
 * bin N_t 
 *
 *
 * Result in OUTFILE will have one row for each bin pair, with the following
 * values:
 *
 * parallel bin low bound, high bound, transverse bin low bound, high bound, 
 * value of ksi, Poisson error
 *
 *
 * Written in 2006 for the benefit of Tuorla Observatory Cosmological
 * Simulations research group by Pauli Pihajoki.
 * 
 * Copyright (c) 2006 Pauli Pihajoki. All rights reserved. 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Datatypes. */
typedef unsigned long long ULONG;
typedef long double LDOUBLE;


/* Struct defs and functionality using them. */
typedef struct _p3 {
	LDOUBLE x, y, z;
} p3;

#define IPROD(a, b) ((a).x*(b).x + (a).y*(b).y + (a).z*(b).z)
#define NORM(a) (sqrt(IPROD((a), (a))))

typedef struct _field {
	ULONG num_pts;
	p3 low_bound, high_bound;
	p3 *pts;
} field;

typedef struct _dataset {
	int num_p_bins;
	LDOUBLE **p_bins;
	int num_t_bins;
	LDOUBLE **t_bins;
	field data;
} dataset;

/* Func defs. */
void read_bins(FILE *f, dataset *storage);
void read_data(FILE *f, field *storage);
void correlate(int correlation_function, dataset *data, field *cmp, FILE *o);

/* A prototype for a correlation function. 
 * Parameters: ratio of N_r/N., DD, DR, RR. */
typedef LDOUBLE (*correlation_function)(LDOUBLE, ULONG, ULONG, ULONG);

/* Supported correlation methods. */
#define PEEBLES_HAUSER		1
#define DAVIS_PEEBLES		2
#define HAMILTON		3
#define LANDY_SZALAY		4
#define NUM_METHODS		4

/* Corresponding functions. */
inline LDOUBLE ph_corr(LDOUBLE rat, ULONG dd, ULONG dr, ULONG rr) {
	return rat*rat*dd/rr - 1;
}
inline LDOUBLE dp_corr(LDOUBLE rat, ULONG dd, ULONG dr, ULONG rr) {
	return rat*dd/dr - 1;
}
inline LDOUBLE ham_corr(LDOUBLE rat, ULONG dd, ULONG dr, ULONG rr) {
	return (LDOUBLE)dd*rr/(dr*dr) - 1;
}
inline LDOUBLE ls_corr(LDOUBLE rat, ULONG dd, ULONG dr, ULONG rr) {
	return 1 + rat*rat*dd/rr - 2*rat*dr/rr;
}

int main(int argc, char **argv) {
	FILE *IN, *OUT, *CMP;
	dataset data;
	field cmp;
	int i, corr_func;

	if (argc != 4 && argc != 5) {
		printf( "correlator: arguments: outputfile datafile "
			"comparison sample file [correlation method]\n"
			"Supported correlation methods:\n"
			"1 => Peebles-Hauser\n"
			"2 => Davis-Peebles\n"
			"3 => Hamilton\n"
			"4 => Landy-Szalay (default)\n");
		exit(EXIT_FAILURE);
	}

	/* Verify the correlation selection. */
	if (argc == 5) {
		corr_func = atoi(argv[4]);
		if (corr_func <= 0 || corr_func > NUM_METHODS) {
			fprintf(stderr, "Invalid correlation method %d\n",
					corr_func);
			exit(EXIT_FAILURE);
		}
	}
	else
		/* Prefer Landy-Szalay correlation function by default. */
		corr_func = 4;
			


	/* Open the file handles, exit on errors.*/
	OUT = fopen(argv[1], "w");
	if (OUT == NULL) {
		perror("fopen");
		exit(EXIT_FAILURE);
	}

	IN = fopen(argv[2], "r");
	if (IN == NULL) {
		perror("fopen");
		exit(EXIT_FAILURE);
	}

	CMP = fopen(argv[3], "r");
	if (CMP == NULL) {
		perror("fopen");
		exit(EXIT_FAILURE);
	}

	/* Read in the data. */
	read_bins(IN, &data);
	read_data(IN, &data.data);
	read_data(CMP, &cmp);

	/* Close input files. */
	fclose(IN);
	fclose(CMP);

	/* Calculate and store the correlation. */
	correlate(corr_func, &data, &cmp, OUT);

	/* Final cleanup. */
	fclose(OUT);
	free(cmp.pts);
	free(data.data.pts);
	for (i=0; i < data.num_p_bins; i++) {
		free(data.p_bins[i]);
	}
	free(data.p_bins);
	for (i=0; i < data.num_t_bins; i++) {
		free(data.t_bins[i]);
	}
	free(data.t_bins);

	return 0;
}

inline void check_malloc(void *p, int line) {
	if (p == NULL) {
		fprintf(stderr, "malloc returned null at line %d", line);
		exit(EXIT_FAILURE);
	}
}

void read_bins(FILE *f, dataset *d) {
	int i;

	/* Process parallel bins first, then transverse bins. */
	
	/* First the number of bins. */
	fscanf(f, "%d", &d->num_p_bins);

	if (d->num_p_bins <= 0) {
		fprintf(stderr, "Invalid number of p_bins: %d\n",
				d->num_p_bins);
		exit(EXIT_FAILURE);
	}

	/* Allocate space. */
	printf("Allocating space for %d p_bins...\n", d->num_p_bins);
	d->p_bins = (LDOUBLE**)malloc(d->num_p_bins * sizeof(LDOUBLE*));
	check_malloc(d->p_bins, __LINE__);

	for (i=0; i < d->num_p_bins; i++) {
		d->p_bins[i] = (LDOUBLE*)malloc(2 * sizeof(LDOUBLE));
		check_malloc(d->p_bins[i], __LINE__); 
	}

	/* Read the p_bins. */
	for (i=0; i < d->num_p_bins; i++) {
		printf("Reading in bin %d... ", i);
		fscanf(f, "%Lf %Lf", &d->p_bins[i][0], &d->p_bins[i][1]);
		printf("%Lf %Lf\n", d->p_bins[i][0], d->p_bins[i][1]);
	}

	/* First the number of bins. */
	fscanf(f, "%d", &d->num_t_bins);

	if (d->num_t_bins <= 0) {
		fprintf(stderr, "Invalid number of t_bins: %d\n",
				d->num_t_bins);
		exit(EXIT_FAILURE);
	}

	/* Allocate space. */
	printf("Allocating space for %d t_bins...\n", d->num_t_bins);
	d->t_bins = (LDOUBLE**)malloc(d->num_t_bins * sizeof(LDOUBLE*));
	check_malloc(d->t_bins, __LINE__);

	for (i=0; i < d->num_t_bins; i++) {
		d->t_bins[i] = (LDOUBLE*)malloc(2 * sizeof(LDOUBLE));
		check_malloc(d->t_bins[i], __LINE__); 
	}

	/* Read the t_bins. */
	for (i=0; i < d->num_t_bins; i++) {
		printf("Reading in bin %d... ", i);
		fscanf(f, "%Lf %Lf", &d->t_bins[i][0], &d->t_bins[i][1]);
		printf("%Lf %Lf\n", d->t_bins[i][0], d->t_bins[i][1]);
	}
}

void read_data(FILE *f, field *d) {
	ULONG i;

	/* First, get the number of points. */
	fscanf(f, "%Lu", &d->num_pts);

	if (d->num_pts <= 0) {
		fprintf(stderr, "Invalid number of data points: %Lu\n",
				d->num_pts);
		exit(EXIT_FAILURE);
	}
	
	/* Allocate space. */
	printf("Allocating space for %Lu points...\n", d->num_pts);
	d->pts = (p3*)malloc((int)d->num_pts * sizeof(p3));
	check_malloc(d->pts, __LINE__);
	
	/* Low and high bounds. */
	printf("Reading high and low bounds...\n");
	fscanf(f, "%Lf %Lf %Lf", 
		&d->low_bound.x, &d->low_bound.y, &d->low_bound.z);
	fscanf(f, "%Lf %Lf %Lf", 
		&d->high_bound.x, &d->high_bound.y, &d->high_bound.z);
	printf("low: %Lf %Lf %Lf high: %Lf %Lf %Lf\n",
		d->low_bound.x, d->low_bound.y, d->low_bound.z,
		d->high_bound.x, d->high_bound.y, d->high_bound.z);

	/* Then get the data. */
	for (i=0; i < d->num_pts; i++) {
		fscanf(f, "%Lf %Lf %Lf", 
			&d->pts[i].x, &d->pts[i].y, &d->pts[i].z);
	}
}

void correlate(int corr_func, dataset *d, field *c, FILE *out) {
	ULONG **dd; 
	ULONG **dr;
	ULONG **rr;
	ULONG i, j;
	int k, l;
	p3 sep, los;
	LDOUBLE para, trans;
	correlation_function cfunc;

	/* Get the right function. */
	switch (corr_func) {
		case PEEBLES_HAUSER:
			cfunc = &ph_corr;
			break;
		case DAVIS_PEEBLES:
			cfunc = &dp_corr;
			break;
		case HAMILTON:
			cfunc = &ham_corr;
			break;
		case LANDY_SZALAY:
		default:
			cfunc = &ls_corr;
			break;
	}

	/* Allocate. */
	dd = (ULONG**)malloc(d->num_p_bins * sizeof(ULONG*));
	check_malloc(dd, __LINE__);
	dr = (ULONG**)malloc(d->num_p_bins * sizeof(ULONG*));
	check_malloc(dr, __LINE__);
	rr = (ULONG**)malloc(d->num_p_bins * sizeof(ULONG*));
	check_malloc(rr, __LINE__);
	for (k=0; k < d->num_p_bins; k++) {
		dd[k] = (ULONG*)malloc(d->num_t_bins * sizeof(ULONG));
		check_malloc(dd[k], __LINE__);
		dr[k] = (ULONG*)malloc(d->num_t_bins * sizeof(ULONG));
		check_malloc(dr[k], __LINE__);
		rr[k] = (ULONG*)malloc(d->num_t_bins * sizeof(ULONG));
		check_malloc(rr[k], __LINE__);
	}

	/* Zero the storages. */
	for (k=0; k < d->num_p_bins; k++) {
		memset(dd[k], 0, d->num_t_bins * sizeof(ULONG));
		memset(dr[k], 0, d->num_t_bins * sizeof(ULONG));
		memset(rr[k], 0, d->num_t_bins * sizeof(ULONG));
	}


	/* Calculate DD and DR pairs, and add for each bin. */
	/* XXX: We require that one pair can only end up in one bin, so we can
	 * stop going through the bins as soon as we find a match. This means
	 * that overlapping bins to achieve some effect will not result in
	 * expected behaviour. */
	printf("DD and DR pairs\n");
	for (i=0; i < d->data.num_pts; i++) {
		/* Print some progress info. */
		printf("%Lu ", i);
		if (i != 0 && (i % 10 == 0))
			printf("\n");

		for (j=i+1; j < d->data.num_pts; j++) {
			sep.x = d->data.pts[i].x - d->data.pts[j].x;
			sep.y = d->data.pts[i].y - d->data.pts[j].y;
			sep.z = d->data.pts[i].z - d->data.pts[j].z;
			los.x = 0.5*(d->data.pts[i].x + d->data.pts[j].x);
			los.y = 0.5*(d->data.pts[i].y + d->data.pts[j].y);
			los.z = 0.5*(d->data.pts[i].z + d->data.pts[j].z);
			para = IPROD(sep, los)/NORM(los);
			trans = IPROD(sep, sep) - para*para;
			for (k=0; k < d->num_p_bins; k++) {
				for (l=0; l < d->num_t_bins; l++) {
					if (	d->p_bins[k][0] <= para &&
						para <= d->p_bins[k][1] &&
						d->t_bins[l][0] <= trans
						&& trans <= d->t_bins[l][1]) {
						dd[k][l] += 2;
						break;
					}
				}
			}
		}

		for (j=0; j < c->num_pts; j++) {
			sep.x = d->data.pts[i].x - c->pts[j].x;
			sep.y = d->data.pts[i].y - c->pts[j].y;
			sep.z = d->data.pts[i].z - c->pts[j].z;
			los.x = 0.5*(d->data.pts[i].x + c->pts[j].x);
			los.y = 0.5*(d->data.pts[i].y + c->pts[j].y);
			los.z = 0.5*(d->data.pts[i].z + c->pts[j].z);
			para = IPROD(sep, los)/NORM(los);
			trans = IPROD(sep, sep) - para*para;
			for (k=0; k < d->num_p_bins; k++) {
				for (l=0; l < d->num_t_bins; l++) {
					if (	d->p_bins[k][0] <= para &&
						para <= d->p_bins[k][1] &&
						d->t_bins[l][0] <= trans
						&& trans <= d->t_bins[l][1]) {
						dr[k][l]++;
						break;
					}
				}
			}
		}
	}
	printf("\n");

	/* Calculate RR pairs. */
	printf("RR pairs\n");
	for (i=0; i < c->num_pts; i++) {
		/* Print some progress info. */
		printf("%Lu ", i);
		if (i != 0 && (i % 10 == 0))
			printf("\n");

		for (j=i+1; j < c->num_pts; j++) {
			sep.x = c->pts[i].x - c->pts[j].x;
			sep.y = c->pts[i].y - c->pts[j].y;
			sep.z = c->pts[i].z - c->pts[j].z;
			los.x = 0.5*(c->pts[i].x + c->pts[j].x);
			los.y = 0.5*(c->pts[i].y + c->pts[j].y);
			los.z = 0.5*(c->pts[i].z + c->pts[j].z);
			para = IPROD(sep, los)/NORM(los);
			trans = IPROD(sep, sep) - para*para;
			for (k=0; k < d->num_p_bins; k++) {
				for (l=0; l < d->num_t_bins; l++) {
					if (	d->p_bins[k][0] <= para &&
						para <= d->p_bins[k][1] &&
						d->t_bins[l][0] <= trans
						&& trans <= d->t_bins[l][1]) {
						rr[k][l] += 2;
						break;
					}
				}
			}
		}
	}
	printf("\n");

	/* Calculate the correlation and Poisson error
	 * for each bin. Then print bin low bound, high bound, ksi and error
	 * on each row. */
	{
		LDOUBLE ksi, error;
		LDOUBLE rat = c->num_pts/(LDOUBLE)d->data.num_pts;
		for (k=0; k < d->num_p_bins; k++) {
			for (l=0; l < d->num_t_bins; l++) {

			printf("Bin %d-%d [%Lg, %Lg] [%Lg, %Lg]\n", 
				k, l, d->p_bins[k][0], d->p_bins[k][1],
				d->t_bins[l][0], d->t_bins[l][1]);

			printf("dd %Lu ", dd[k][l]);
			printf("dr %Lu ", dr[k][l]);
			printf("rr %Lu\n", rr[k][l]);

			if (rr[k][l] == 0) {
				/* XXX: Something like this is what _should_
				 * be done, however, it would require
				 * enabling GNU extensions or somesuch. */
				/*
				ksi = error = NAN;
				*/
				fprintf(out, "%Lg %Lg %Lg %Lg %s %s\n",
					d->p_bins[k][0], d->p_bins[k][1],
					d->t_bins[l][0], d->t_bins[l][1],
					"NaN", "NaN");
			}
			else {
				/* Use the given estimator. */
				ksi = cfunc(rat, dd[k][l], dr[k][l], rr[k][l]);

				/* For the Poisson errors, we have:
				 * (1+ksi)/sqrt(x), where x can be DD,
				 * (N/N_r)*DR or (N/N_r)^2*RR. Of these,
				 * Martinez & Saar suggest using one of the
				 * latter two. Of these, again the latter
				 * is more convenient, as the RR counts should
				 * be consistently large. */
				/*
				error = (1 + ksi)/sqrt(dd[k]);
				error = (1 + ksi)/sqrt(dr[k]/rat);
				*/
				error = (1 + ksi)/sqrt(rr[k][l]/(rat*rat));

				fprintf(out, "%Lg %Lg %Lg %Lg %Lg %Lg\n",
					d->p_bins[k][0], d->p_bins[k][1],
					d->t_bins[l][0], d->t_bins[l][1],
					ksi, error);
			}

			}
		}
	}

	/* Cleanup. */
	for (k=0; k < d->num_t_bins; k++) {
		free(dd[k]);
		free(dr[k]);
		free(rr[k]);
	}
	free(dd);
	free(dr);
	free(rr);
}


