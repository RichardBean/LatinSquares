/* This is a program to check through n x n (n <= 8) Latin squares in a file to see if any of them have critical sets of 
 * n^2 / 4 - 1 (which can be modified to n^2 / 4 in the line with GRB_LESS_EQUAL)
 *
 * This program can be used to verify the results in the paper
 * "The size of the smallest uniquely completable set in order 8 Latin squares" by Richard Bean
 * in Journal of Combinatorial Mathematics and Combinatorial Computing, v. 52 (2005), pp. 159–168.
 * For instance "tradegu 8list 1 283657 8 3 24" followed by three runs for the unsolved squares would verify
 * the results in the paper.
 *
 * This version is from 9 July 2011. The next improvement is to use Knuth's Dancing Links algorithm for 
 * Latin square completions.
 * This might be useful if anyone tries to use 4-row/col/elt trades to solve the conjecture about 8x8 squares
 * not having any critical sets of size 16 (except for the square based on Z_8) in the paper,
 * because it would make finding the trades much quicker.
 *
 * The program requires the gurobi library and header files installed.
 * Free academic licenses for gurobi are available from http://www.gurobi.com/html/academic.html
 * Compile with: gcc -O3 -o tradegu tradegu.c -lgurobi45
 * Usage: tradegu filename linestart lineend size k limit
 * where: linestart = line to start at, lineend = line to end at (first line is 1)
 * size = order of Latin squares in file
 * k = maximum number of rows / columns / elements in trades from Latin square to consider
 * limit = limit of maximum size of trade to use in MIP
 *
 * The parameters ... 4 9 produce the same results as parameters ... 3 9 (finds the same trades - just takes longer)
 * The same applies to ... 3 6 and ... 2 6.
 *
 * Run times on a 2011 Macbook Pro with 2.3 GHz Intel Core i7 (2820QM - 4 cores, 8 threads)
 *
 * "tradegu 6list 1 12 6 3 10" takes 1 second (solves 12 out of 12)
 *
 * (in 7list, the 6th square has no trades smaller than size 9)
 * "tradegu 7list 1 147 7 2 4" takes <1 second (solves 1 out of 147 - the main class with 42 intercalates ie STS(7))
 * "tradegu 7list 1 147 7 2 6" takes 2 seconds (solves 5 out of 147)
 * "tradegu 7list 1 147 7 2 8" takes 3 seconds (solves 6 out of 147)
 * "tradegu 7list 1 147 7 2 10" takes 5 seconds (solves 7 out of 147)
 * "tradegu 7list 1 147 7 2 14" takes 8 seconds (solves 7 out of 147)
 * "tradegu 7list 1 147 7 3 6" takes 3 seconds (solves 5 out of 147)
 * "tradegu 7list 1 147 7 3 8" takes 13 seconds (solves 18 out of 147)
 * "tradegu 7list 1 147 7 3 9" takes 23 seconds (solves 59 out of 147)
 * "tradegu 7list 1 147 7 3 10" takes 59 seconds (solves 146 out of 147, misses out #90 which takes <1s with "... 7 3 11")
 * "tradegu 7list 1 147 7 3 11" takes 71 seconds (solves 147 out of 147)
 * "tradegu 7list 1 147 7 3 12" takes 100 seconds
 * "tradegu 7list 1 147 7 3 21" takes 197 seconds
 * "tradegu 7list 1 147 7 4 9" takes 234 seconds (solves 59 out of 147)
 * "tradegu 7list 1 147 7 4 10" takes 288 seconds (solves 147 out of 147)
 * 
 * "tradegu 8list 1 283657 8 2 4" takes 307 seconds (solves 528 of 283657)
 *
 * (in 8list, the last three squares have no intercalates, but all have trades of size 6)
 * "tradegu 8list 283655 283657 8 3 8" takes 1 second (solves 0 of 3)
 * "tradegu 8list 283655 283657 8 3 9" takes 10 seconds (solves 0 of 3)
 * "tradegu 8list 283655 283657 8 3 10" takes 275 seconds (solves 3 of 3)
 * "tradegu 8list 283655 283657 8 4 6" takes 33 seconds (solves 0 of 3)
 * "tradegu 8list 283655 283657 8 4 8" takes 35 seconds (solves 0 of 3)
 * "tradegu 8list 283655 283657 8 4 9" takes 43 seconds (solves 0 of 3)
 *
 * "tradegu 8list 1 1000 8 2 8" takes 68 seconds (solves 467 of 1000)
 * "tradegu 8list 1 1000 8 2 16" takes 108 seconds (solves 516 of 1000)
 * "tradegu 8list 1 1000 8 3 8" takes 141 seconds (solves 808 of 1000)
 * "tradegu 8list 1 1000 8 3 9" takes 147 seconds (solves 858 of 1000)
 * "tradegu 8list 1 1000 8 3 10" takes 228 seconds (solves 998 of 1000)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gurobi_c.h"

#define SIZE 8
int             s1[SIZE][SIZE];
int             limit;

int             trades = 0;

struct trade {
	int             on;
	int             filled;
	unsigned long   sq;
};

struct trade   *tlist;

struct bitmap {
	int             count;
	int             values[SIZE + 1];
};

int             fill(int s[SIZE][SIZE], int level, int pos, int size);

void
add(unsigned long d, int size, int filled)
{
	int             i;
	struct trade   *more_trades;

	if (trades == 0) {
		tlist = (struct trade *) malloc(sizeof(struct trade));
		tlist[0].sq = d;
		tlist[0].on = 1;
		tlist[0].filled = filled;
		trades++;
		return;
	}
	for (i = 0; i < trades; i++) {
		if ((d & tlist[i].sq) == tlist[i].sq)
			return;
		if ((d & tlist[i].sq) == d)
			tlist[i].on = 0;
	}
	trades++;
	more_trades = (struct trade *) realloc(tlist, sizeof(struct trade) * trades);
	if (!more_trades) {
		printf("realloc failed\n");
		exit(1);
	} else
		tlist = more_trades;
	tlist[trades - 1].sq = d;
	tlist[trades - 1].on = 1;
	tlist[trades - 1].filled = filled;
	return;
}

struct bitmap
testone(int s[SIZE][SIZE], int a, int b, int size)
{
	/*
	 * A function for testing if the empty cell in row a, column b of the
	 * Latin square in the array s has a forced completion.
	 */

	struct bitmap   bm;
	int             t[SIZE + 1];
	int             that;
	int             i;

	/* clear the test and values array */
	for (i = 1; i <= size; i++) {
		t[i] = 0;
		bm.values[i] = 0;
	}
	bm.count = 0;

	/*
	 * for every element in the same row and column as the element being
	 * tested, set the corresponding test array element to be 1
	 */
	for (i = 0; i < size; i++)
		t[s[a][i]] = t[s[i][b]] = 1;

	/*
	 * keep track of all the different elements used in that row and
	 * column
	 */
	for (i = 1; i <= size; i++)
		if (!t[i]) {
			bm.values[i] = 1;
			that = i;
			bm.count++;
		}
	/* return those values */
	return bm;
}

int
printt(unsigned long t, int size, int ind[SIZE * SIZE], double val[SIZE * SIZE])
{
	int             x, y, z = 0;
	unsigned long   o = 1;

	for (y = 0; y < size; y++)
		for (x = 0; x < size; x++)
			if (t & (o << (y * size + x))) {
				ind[z++] = y * size + x;
				val[y * size + x] = 1;
			}
}

int
vfill(int *v, int n, int k)
{
	int             s[SIZE][SIZE], i, j, a;
	memcpy(s, s1, sizeof(s));
	for (i = 0; i < n; i++)
		for (a = 0; a < k; a++)
			s[v[a]][i] = 0;
	fill(s, n * (n - k), 0, n);

	memcpy(s, s1, sizeof(s));
	for (i = 0; i < n; i++)
		for (a = 0; a < k; a++)
			s[i][v[a]] = 0;
	fill(s, n * (n - k), 0, n);

	memcpy(s, s1, sizeof(s));
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (a = 0; a < k; a++)
				if (s[i][j] - 1 == v[a])
					s[i][j] = 0;
	fill(s, n * (n - k), 0, n);
}

int
fill(int s[SIZE][SIZE], int level, int pos, int size)
{
	struct bitmap   ret;
	int             a, b, c, poss = 0;

	if (level == size * size) {
		unsigned long   d = 0;
		unsigned long   o = 1;
		int             t = 0;
		for (a = 0; a < size; a++)
			for (b = 0; b < size; b++)
				if (s[a][b] != s1[a][b]) {
					d |= (unsigned long) (o << (a * size + b));
					t++;
				}
		if (d && t <= limit)
			add(d, size, t);
		return 1;
	}
	a = pos % size;
	b = pos / size;
	while (s[b][a]) {
		a++;
		if (a == size) {
			b++;
			a = 0;
		}
	}


	ret = testone(s, b, a, size);
	for (c = 1; c <= size; c++) {
		if (ret.values[c]) {
			s[b][a] = c;
			poss += fill(s, level + 1, b * size + a + 1, size);
			s[b][a] = 0;
		}
	}
	return poss;
}

int
main(int argc, char **argv)
{
	int             s[SIZE][SIZE], i, j, n, len, k, *v, a, x, linestart, lineend,
	                line;
	FILE           *file;
	char            str[100];

	GRBenv         *env = NULL;
	GRBmodel       *model = NULL;
	int             error = 0;
	double          sol[SIZE * SIZE];
	int             ind[SIZE * SIZE];
	double          val[SIZE * SIZE];
	double          obj[SIZE * SIZE];
	char            vtype[SIZE * SIZE];
	int             optimstatus;
	double          objval;

	if (argc != 7) {
		printf("usage: %s filename linestart lineend size k limit\n", argv[0]);
		exit(0);
	}
	if ((file = fopen(argv[1], "r")) == NULL) {
		printf("failed to open %s\n", argv[1]);
		exit(0);
	}
	linestart = atoi(argv[2]);
	lineend = atoi(argv[3]);
	n = atoi(argv[4]);
	k = atoi(argv[5]);
	limit = atoi(argv[6]);

	v = malloc((n + 2) * sizeof(int));	/* for doing n choose k soon */

	/* go to linestart */

	for (i = 1; i < linestart; i++)
		for (j = 0; j < n; j++)
			fscanf(file, "%s", str);

	/* Create environment */

	error = GRBloadenv(&env, NULL);
	if (error || env == NULL) {
		fprintf(stderr, "Error: could not create environment\n");
		exit(1);
	}
	error = GRBsetintparam(env, "OutputFlag", 0);
	if (error)
		goto QUIT;

	for (line = linestart; line <= lineend; line++) {
		printf("%d ", line);
		fflush(stdout);
		for (i = 0; i < n; i++) {
			fscanf(file, "%s", str);
			for (j = 0; j < n; j++)
				s1[i][j] = str[j] - '0';
		}

		/* do n choose k to find the trades */


		for (i = 0; i < k; i++)
			v[i] = i;
		vfill(v, n, k);
		v[k] = n;

		while (v[0] < n - k) {
			j = -1;
			do {
				j++;
			} while (v[j + 1] <= v[j] + 1);

			v[j]++;

			for (i = 0; i < j; i++)
				v[i] = i;

			vfill(v, n, k);
		}


		/* Create an empty model */

		error = GRBnewmodel(env, &model, "mip1", 0, NULL, NULL, NULL, NULL, NULL);
		if (error)
			goto QUIT;


		/* Add variables */

		for (i = 0; i < n * n; i++) {
			obj[i] = 1;
			vtype[i] = GRB_BINARY;
		}
		error = GRBaddvars(model, n * n, 0, NULL, NULL, NULL, obj, NULL, NULL, vtype,
				   NULL);
		if (error)
			goto QUIT;

		/* Integrate new variables */

		error = GRBupdatemodel(model);
		if (error)
			goto QUIT;

		/*
		 * First constraint: we're looking for a solution of size <=
		 * n*n/4 - 1
		 */

		for (i = 0; i < n * n; i++) {
			ind[i] = i;
			val[i] = 1;
		}

		error = GRBaddconstr(model, n * n, ind, val, GRB_LESS_EQUAL, n * n / 4 - 1, NULL);
		if (error)
			goto QUIT;

		/*
		 * other constraints: must have at least one entry in each
		 * trade
		 */

		for (i = 0; i < trades; i++) {
			if (tlist[i].on) {
				printt(tlist[i].sq, n, ind, val);
				error = GRBaddconstr(model, tlist[i].filled, ind, val, GRB_GREATER_EQUAL, 1.0, NULL);
				if (error)
					goto QUIT;
			}
		}

		/* Optimize model */

		error = GRBoptimize(model);
		if (error)
			goto QUIT;

		/* Write model to 'mip1.lp' */

		/*
		 * error = GRBwrite(model, "mip1.lp"); if (error) goto QUIT;
		 */

		error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
		if (error)
			goto QUIT;

		if (optimstatus == GRB_OPTIMAL) {
			error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL,
					      &objval);
			if (error)
				goto QUIT;
			printf("%d\n", (int) objval); /* solution found */
		} else if (optimstatus == GRB_INFEASIBLE) {
			printf("infeasible\n");
		} else {
			printf("stopped_early%d\n", optimstatus); 
		}
		fflush(stdout);

		/* Free model */
	
		GRBfreemodel(model);
		if (trades) free(tlist);
		trades = 0;
	}

QUIT:

	/* Error reporting */

	if (error) {
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(1);
	}
	/* Free environment */

	GRBfreeenv(env);
	fclose(file);

	return 0;
}
