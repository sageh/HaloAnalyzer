/* Calculate different two-point correlation functions using parallelized
 * computing.
 *
 * Usage: correlator OUTFILE DATAFILE CMPFILE method
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
 * number of scale bins (N_s)
 * bin 1 (2 components, low and high bounds, separated by whitespace)
 * .
 * .
 * bin N_s 
 *
 *
 * Result in OUTFILE will have one row for each bin, with the following
 * values:
 *
 * bin low bound, high bound, value of ksi, Poisson error
 *
 *
 * Written in 2006 for the benefit of Tuorla Observatory Cosmological
 * Simulations research group by Pauli Pihajoki.
 * 
 * Copyright (c) 2006 Pauli Pihajoki. All rights reserved. 
 */

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Datatypes. */
typedef unsigned long ULONG;
typedef double LDOUBLE;

/* Struct defs and inlines using them. */
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
	int num_bins;
	LDOUBLE **bins;
	field data;
} dataset;

/* Func defs. */
int read_bins(FILE *f, dataset *storage);
ULONG read_data(FILE *f, field *storage);
void correlate(int correlation_function, dataset *data, field *cmp, FILE *o,
		int processor_amount, int this_processor_id);

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

/* Malloc checker. */
inline void check_malloc(void *p, int line) {
	if (p == NULL) {
		fprintf(stderr, "malloc returned null at line %d", line);
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char **argv) {
	FILE *IN, *OUT, *CMP;
	dataset data;
	field cmp;
	int i; 
	/* Use Landy-Szalay correlation function by default. */
	int corr_func = 4;
	ULONG j;
	/* MPI specific variables */
	int mpi_rank, mpi_size;

	/* Initialize MPI library. */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	if (mpi_rank == 0) 
		printf("Initialized MPI. Using %d nodes.\n", mpi_size);

	/* If we are the main process, then handle file I/O. */
	if (mpi_rank == 0) {
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
				fprintf(stderr, "Invalid correlation "
						"method %d\n", corr_func);
				exit(EXIT_FAILURE);
			}
		}

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
	}
	else
		IN = OUT = CMP = NULL;

	/* Send the data to every other process.
	 * XXX: For now, we don't have to worry about sending the
	 * low and high bounds of the data, because that information is
	 * not being used _at the moment_. */
#ifdef DEBUG
	printf("p%d before bcast: num_bins %d\n", mpi_rank, data.num_bins);
	printf("p%d before bcast: data.num_pts %lu\n", 
			mpi_rank, data.data.num_pts);
	printf("p%d before bcast: cmp.num_pts: %lu\n", mpi_rank, cmp.num_pts);
#endif

	MPI_Bcast(&data.num_bins, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&data.data.num_pts, 1, MPI_UNSIGNED_LONG, 0, 
			MPI_COMM_WORLD);
	MPI_Bcast(&cmp.num_pts, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

	/* Print some some debug info if we are a slave process. */
#ifdef DEBUG
	printf("p%d after bcast: num_bins %d\n", mpi_rank, data.num_bins);
	printf("p%d after bcast: data.num_pts %lu\n", 
			mpi_rank, data.data.num_pts);
	printf("p%d after bcast: cmp.num_pts: %lu\n", mpi_rank, cmp.num_pts);
#endif

	/* If we are not root, then then malloc some space before receving
	 * the actual data. */
	if (mpi_rank != 0) {
		/* Bins. */
		printf("p%d: Allocating space for %d bins...\n", 
				mpi_rank, data.num_bins);
		data.bins = (LDOUBLE**)malloc(data.num_bins * sizeof(LDOUBLE*));
		check_malloc(data.bins, __LINE__);

		for (i=0; i < data.num_bins; i++) {
			data.bins[i] = (LDOUBLE*)malloc(2 * sizeof(LDOUBLE));
			check_malloc(data.bins[i], __LINE__); 
		}

		/* Data points. */
		printf("p%d: Allocating space for %lu points...\n", 
				mpi_rank, data.data.num_pts);
		data.data.pts = (p3*)malloc(
				(int)data.data.num_pts * sizeof(p3));
		check_malloc(data.data.pts, __LINE__);

		/* Comparison field. */
		printf("p%d: Allocating space for %lu points...\n", 
				mpi_rank, cmp.num_pts);
		cmp.pts = (p3*)malloc((int)cmp.num_pts * sizeof(p3));
		check_malloc(cmp.pts, __LINE__);
	}

	/* Broadcast bins, data set and the comparison field. */
	for (i=0; i < data.num_bins; i++) {
		MPI_Bcast(&data.bins[i][0], 1, MPI_DOUBLE, 0, 
				MPI_COMM_WORLD);
		MPI_Bcast(&data.bins[i][1], 1, MPI_DOUBLE, 0, 
				MPI_COMM_WORLD);

		/* Debug info if we're not root. */
#ifdef DEBUG
		if (mpi_rank != 0) 
			printf("p%d received: bin %d [%lg, %lg]\n",
				mpi_rank, i, data.bins[i][0], data.bins[i][1]);
#endif
	}

	for (j=0; j < data.data.num_pts; j++) {
		MPI_Bcast(&data.data.pts[j].x, 1, MPI_DOUBLE, 0,
				MPI_COMM_WORLD);
		MPI_Bcast(&data.data.pts[j].y, 1, MPI_DOUBLE, 0,
				MPI_COMM_WORLD);
		MPI_Bcast(&data.data.pts[j].z, 1, MPI_DOUBLE, 0,
				MPI_COMM_WORLD);
	}

	for (j=0; j < cmp.num_pts; j++) {
		MPI_Bcast(&cmp.pts[j].x, 1, MPI_DOUBLE, 0,
				MPI_COMM_WORLD);
		MPI_Bcast(&cmp.pts[j].y, 1, MPI_DOUBLE, 0,
				MPI_COMM_WORLD);
		MPI_Bcast(&cmp.pts[j].z, 1, MPI_DOUBLE, 0,
				MPI_COMM_WORLD);
	}

	/* Calculate and store the correlation. */
	correlate(corr_func, &data, &cmp, OUT, mpi_size, mpi_rank);

	/* Final cleanup. */
	if (mpi_rank == 0)
		fclose(OUT);

	free(cmp.pts);
	free(data.data.pts);
	for (i=0; i < data.num_bins; i++) {
		free(data.bins[i]);
	}
	free(data.bins);

	MPI_Finalize();

	return 0;
}

int read_bins(FILE *f, dataset *d) {
	int i;
	
	/* First the number of bins. */
	fscanf(f, "%d", &d->num_bins);

	if (d->num_bins <= 0) {
		fprintf(stderr, "Invalid number of bins: %d\n",
				d->num_bins);
		exit(EXIT_FAILURE);
	}

	/* Allocate space. */
	printf("Allocating space for %d bins...\n", d->num_bins);
	d->bins = (LDOUBLE**)malloc(d->num_bins * sizeof(LDOUBLE*));
	check_malloc(d->bins, __LINE__);

	for (i=0; i < d->num_bins; i++) {
		d->bins[i] = (LDOUBLE*)malloc(2 * sizeof(LDOUBLE));
		check_malloc(d->bins[i], __LINE__); 
	}

	/* Read the bins. */
	for (i=0; i < d->num_bins; i++) {
		printf("Reading in bin %d... ", i);
		fscanf(f, "%lg %lg", &d->bins[i][0], &d->bins[i][1]);
		printf("%lg %lg\n", d->bins[i][0], d->bins[i][1]);
	}

	/* Return the number of bins read. */
	return d->num_bins;
}

ULONG read_data(FILE *f, field *d) {
	ULONG i;

	/* First, get the number of points. */
	fscanf(f, "%lu", &d->num_pts);

	if (d->num_pts <= 0) {
		fprintf(stderr, "Invalid number of data points: %lu\n",
				d->num_pts);
		exit(EXIT_FAILURE);
	}
	
	/* Allocate space. */
	printf("Allocating space for %lu points...\n", d->num_pts);
	d->pts = (p3*)malloc((int)d->num_pts * sizeof(p3));
	check_malloc(d->pts, __LINE__);
	
	/* Low and high bounds. */
	printf("Reading high and low bounds...\n");
	fscanf(f, "%lg %lg %lg", 
		&d->low_bound.x, &d->low_bound.y, &d->low_bound.z);
	fscanf(f, "%lg %lg %lg", 
		&d->high_bound.x, &d->high_bound.y, &d->high_bound.z);
	printf("low: %lg %lg %lg high: %lg %lg %lg\n",
		d->low_bound.x, d->low_bound.y, d->low_bound.z,
		d->high_bound.x, d->high_bound.y, d->high_bound.z);

	/* Then get the data. */
	for (i=0; i < d->num_pts; i++) {
		fscanf(f, "%lg %lg %lg", 
			&d->pts[i].x, &d->pts[i].y, &d->pts[i].z);
	}

	/* Return the number of data points read. */
	return d->num_pts;
}

void correlate(int corr_func, dataset *d, field *c, FILE *out,
		int mpi_size, int mpi_rank) {
	ULONG *dd, *dd_; 
	ULONG *dr, *dr_;
	ULONG *rr, *rr_;
	ULONG i, j;
	int k;
	p3 sep;
	LDOUBLE l;
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
	dd = (ULONG*)malloc(d->num_bins * sizeof(ULONG));
	check_malloc(dd, __LINE__);
	dr = (ULONG*)malloc(d->num_bins * sizeof(ULONG));
	check_malloc(dr, __LINE__);
	rr = (ULONG*)malloc(d->num_bins * sizeof(ULONG));
	check_malloc(rr, __LINE__);

	/* Zero the storages. */
	memset(dd, 0, d->num_bins * sizeof(ULONG));
	memset(dr, 0, d->num_bins * sizeof(ULONG));
	memset(rr, 0, d->num_bins * sizeof(ULONG));

	/* Only allocate the final sum storages for the root process. */
	if (mpi_rank == 0) {
		dd_ = (ULONG*)malloc(d->num_bins * sizeof(ULONG));
		check_malloc(dd_, __LINE__);
		dr_ = (ULONG*)malloc(d->num_bins * sizeof(ULONG));
		check_malloc(dr_, __LINE__);
		rr_ = (ULONG*)malloc(d->num_bins * sizeof(ULONG));
		check_malloc(rr_, __LINE__);
		memset(dd_, 0, d->num_bins * sizeof(ULONG));
		memset(dr_, 0, d->num_bins * sizeof(ULONG));
		memset(rr_, 0, d->num_bins * sizeof(ULONG));
	}
	else 
		dd_ = dr_ = rr_ = NULL;


	/* Calculate DD and DR pairs, and add for each bin. */
	/* XXX: We require that one pair can only end up in one bin, so we can
	 * stop going through the bins as soon as we find a match. This means
	 * that overlapping bins to achieve some effect will not result in
	 * expected behaviour. */
#ifdef PROGRESS
	if (mpi_rank == 0)
		printf("DD and DR pairs\n");
#endif
	for (i=mpi_rank; i < d->data.num_pts; i+=mpi_size) {
		/* Print some progress info, but only if we are the 
		 * root process. */
#ifdef PROGRESS
		if (mpi_rank == 0) {
			printf("%lu ", i);
			if (i != 0 && ((i/mpi_size) % 5 == 0)) {
				printf("\n");
				fflush(stdout);
			}
		}
#endif

		for (j=i+1; j < d->data.num_pts; j++) {
			sep.x = d->data.pts[i].x - d->data.pts[j].x;
			sep.y = d->data.pts[i].y - d->data.pts[j].y;
			sep.z = d->data.pts[i].z - d->data.pts[j].z;
			l = NORM(sep);
			for (k=0; k < d->num_bins; k++) {
				if (d->bins[k][0] <= l
					&& l <= d->bins[k][1]) {
					dd[k]+=2;
					break;
				}
			}
		}

		for (j=0; j < c->num_pts; j++) {
			sep.x = d->data.pts[i].x - c->pts[j].x;
			sep.y = d->data.pts[i].y - c->pts[j].y;
			sep.z = d->data.pts[i].z - c->pts[j].z;
			l = NORM(sep);
			for (k=0; k < d->num_bins; k++) {
				if (d->bins[k][0] <= l
					&& l <= d->bins[k][1]) {
					dr[k]++;
					break;
				}
			}
		}
	}
#ifdef PROGRESS
	if (mpi_rank == 0)
		printf("\n");
#endif

	/* Calculate RR pairs. */
#ifdef PROGRESS
	if (mpi_rank == 0)
		printf("RR pairs\n");
#endif
	for (i=mpi_rank; i < c->num_pts; i+=mpi_size) {
		/* Print some progress info, but only if we are the 
		 * root process. */
#ifdef PROGRESS
		if (mpi_rank == 0) {
			printf("%lu ", i);
			if (i != 0 && ((i/mpi_size) % 5 == 0)) {
				printf("\n");
				fflush(stdout);
			}
		}
#endif

		for (j=i+1; j < c->num_pts; j++) {
			sep.x = c->pts[i].x - c->pts[j].x;
			sep.y = c->pts[i].y - c->pts[j].y;
			sep.z = c->pts[i].z - c->pts[j].z;
			l = NORM(sep);
			for (k=0; k < d->num_bins; k++) {
				if (d->bins[k][0] <= l
					&& l <= d->bins[k][1]) {
					rr[k]+=2;
					break;
				}
			}
		}
	}
#ifdef PROGRESS
	if (mpi_rank == 0)
		printf("\n");
#endif


	/* Gather all the pair counts from the slave processes and sum
	 * them at the root process. */
#ifdef DEBUG
	for (k=0; k < d->num_bins; k++) {
		/* DEBUG: Check what we have before reducing. */
		printf("p%d: Bin %d [%lg, %lg] ", mpi_rank,
			k, d->bins[k][0], d->bins[k][1]);
		printf("dd %lu ", dd[k]);
		printf("dr %lu ", dr[k]);
		printf("rr %lu\n", rr[k]);
	}
#endif

	MPI_Reduce(dd, dd_, d->num_bins, MPI_UNSIGNED_LONG, 
			MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(dr, dr_, d->num_bins, MPI_UNSIGNED_LONG, 
			MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(rr, rr_, d->num_bins, MPI_UNSIGNED_LONG, 
			MPI_SUM, 0, MPI_COMM_WORLD);

	/* Calculate the correlation and Poisson error
	 * for each bin. Then print bin low bound, high bound, ksi and error
	 * on each row. */
	if (mpi_rank == 0) {
		LDOUBLE ksi, error;
		LDOUBLE rat = c->num_pts/(LDOUBLE)d->data.num_pts;
		for (k=0; k < d->num_bins; k++) {
			printf("Bin %d [%lg, %lg] ", 
				k, d->bins[k][0], d->bins[k][1]);
			printf("dd %lu ", dd_[k]);
			printf("dr %lu ", dr_[k]);
			printf("rr %lu\n", rr_[k]);
			if (rr_[k] == 0) {
				/* XXX: Something like this is what _should_
				 * be done, however, it would require
				 * enabling GNU extensions or somesuch. */
				/*
				ksi = error = NAN;
				*/
				fprintf(out, "%lg %lg %s %s\n",
					d->bins[k][0], d->bins[k][1],
					"nan", "nan");
			}
			else {
				/* Use the given estimator. */
				ksi = cfunc(rat, dd_[k], dr_[k], rr_[k]);

				/* For the Poisson errors, we have:
				 * (1+ksi)/sqrt(x), where x can be DD,
				 * (N/N_r)*DR or (N/N_r)^2*RR. Of these,
				 * Martinez & Saar suggest using one of the
				 * latter two. Of these, again the latter
				 * is more convenient, as the RR counts should
				 * be consistently large. */
				/*
				error = (1 + ksi)/sqrt(dd_[k]);
				error = (1 + ksi)/sqrt(dr_[k]/rat);
				*/
				error = (1 + ksi)/sqrt(rr_[k]/(rat*rat));

				fprintf(out, "%lg %lg %lg %lg\n",
					d->bins[k][0], d->bins[k][1],
					ksi, error);
			}
		}
	}

	/* Cleanup. */
	free(dd);
	free(dr);
	free(rr);

	if (mpi_rank == 0) {
		free(dd_);
		free(dr_);
		free(rr_);
	}
}


