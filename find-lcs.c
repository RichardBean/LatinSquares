/* find-lcs, a program for finding large critical sets
   by beginning with one of 3 types of PLS (partial Latin square), finding all 
   completions, and searching for a critical set within these LSs
   by removing elements according to one of 4 patterns. */

/* May 13 2002, rwb (pig's tail) eskimo.com */

/* complex fill2() when testing for unique completion, simple fill()
when generating all completions.  simple fill() used from PLS with 1..n
in first row and column */

/* example running times on an athlon 4 1200 mhz */
/* where x:y is given, x=number of intercalates, y=size of cs found */

/* find-lcs 6 18  0 0 0 - 111s */
/* find-lcs 6 18 27 0 0 - 1.7s */
/* find-lcs 6 18 27 1 0 - 0.14s - nothing found (9408->20 squares) */
/* find-lcs 6 18 27 1 1 - 0.3s */
/* find-lcs 6 18 27 1 2 - 0.34s - nothing found (9408->20 squares) */

/* find-lcs 7 25  0 0 0 - 3.5s - 26:25 */
/* find-lcs 7 25  0 0 1 - >10m         */
/* find-lcs 7 25  0 1 0 - >74m         */  /* paradoxical? - see 7 25 0 0 0 */
/* find-lcs 7 25  0 2 0 -  99s - nothing found (10752 squares) */
/* find-lcs 7 25  0 2 1 - 181s */
/* find-lcs 7 25  0 2 1 - 353s - nothing found (10752 squares) */
/* find-lcs 7 25 26 1 0 -  31s - 42:25 */
/* find-lcs 7 25 42 0 0 -  20s */
/* find-lcs 7 25 42 0 1 -  24s */
/* find-lcs 7 25 42 1 0 -   9s */
/* find-lcs 7 25 42 1 1 - 2.5s */

/* find-lcs 8 37  0 0 0 - 0.02s */

/* find-lcs 9  0  0 0 0 - 0.15s - 41:40 */
/* find-lcs 9  0  0 0 1 - 0.46s - 41:38 */
/* find-lcs 9  0  0 0 2 - 0.46s - 41:32 */
/* find-lcs 9 39  0 0 1 - >4m           */
/* find-lcs 9 41  0 0 0 - 1.4s  - 36:41 */
/* find-lcs 9 41  0 0 1 - >11m          */
/* find-lcs 9 42  0 0 0 - 20.5s - 37:42 */
/* find-lcs 9 42  0 2 0 - 3.5s  - 43:42 */
/* find-lcs 9 43  0 0 0 -               */
/* find-lcs 9 43  0 1 0 - >3m           */
/* find-lcs 9 43  0 2 0 - 13.5s - 46:43 */
/* find-lcs 9 44  0 0 0 -               */
/* find-lcs 9 44  0 2 0 - 964s          */
/* find-lcs 9 44 53 1 0 - ~1m           */
/* find-lcs 9 44 53 2 0 - ~3s   - 53:44 */
/* find-lcs 9 44 53 2 1 - >25m          */
/* find-lcs 9 44 53 2 2 -               */
/* find-lcs 9 45  0 0 0 -               */
/* find-lcs 9 45  0 0 1 -               */
/* find-lcs 9 45  0 0 2 -               */
/* find-lcs 9 45  0 1 0 - ~3h?  - 37:45 */
/* find-lcs 9 45  0 1 1 -               */
/* find-lcs 9 45  0 1 2 -               */

/* find-lcs 10 0   0 0 0 - 0.23s - 57:51 */
/* find-lcs 10 0   0 0 1 - 1.19s - 57:43 */
/* find-lcs 10 0   0 0 1 - 0.63s - 57:46 */
/* find-lcs 10 52  0 0 0 - 1.85s - 57:52 */
/* find-lcs 10 53  0 0 0 - >6m   -       */
/* find-lcs 10 53 84 2 0 - >9m   -       */

#define SIZE 26
#include <stdio.h>
#include <sys/time.h>		/* for srandom() */

void recurse (int s[SIZE][SIZE], int cs, int size, int ic);
void recursemax (int s[SIZE][SIZE], int cs, int size, int ic);
void recursemin (int s[SIZE][SIZE], int cs, int size, int ic);
void recursernd (int s[SIZE][SIZE], int cs, int size, int ic);
int fill2 (int s[SIZE][SIZE], int level, int pos, int size);
int fill (int s[SIZE][SIZE], int level, int pos, int size);
int icount (int s[SIZE][SIZE], int size);
struct bitmap testone (int s[SIZE][SIZE], int a, int b, int size);
void print (int s[SIZE][SIZE], int size);
void display (int array[200][200]);
int max = 0, minsize, mininter, recflag;
int count = 0;
int array[200][200];		/* array */

main (int argc, char **argv)
{
  int s[SIZE][SIZE], size, i, j, k;
  int use;
  struct timeval tp;
  struct timezone tzp;
  if (argc != 6)
    {
      printf
	("usage: %s order-of-LS minimum-size-wanted minimum-intercalates useflag recflag\n",
	 argv[0]);
      printf
	("useflag: 0 for empty square, 1 for 1..n in first row+col, 2 for 1..n in first row+col and 1s on main diagonal\n");
      printf
	("recflag: 0 for simple recursion, 1 for remove (i,j) where x_{ij} is max, 2 for where x_{ij} is min, 3 is random\n");
      exit (0);
    }
  size = atoi (argv[1]);
  minsize = atoi (argv[2]);
  mininter = atoi (argv[3]);
  use = atoi (argv[4]);
  recflag = atoi (argv[5]);
  memset (s, 0, sizeof (s));
  memset (array, 0, sizeof (array));

  /* use standard form */
  for (i = 0; i < SIZE; i++)
    {
      if (use == 1 || use == 2)
	{
	  s[0][i] = i + 1;
	  s[i][0] = i + 1;
	}
      if (use == 2)
	s[i][i] = 1;
    }
  if (recflag == 3)
    {
      gettimeofday (&tp, &tzp);
      srandom ((int) tp.tv_usec);
    }
  /* call fill() with different arguments depending on the invocation
     on the command line */
  if (use == 0)
    fill (s, 0, 0, size);
  if (use == 1)
    fill (s, size + size - 1, 0, size);
  if (use == 2)
    fill (s, size + size + size - 2, 0, size);
}

int
icount (int s[SIZE][SIZE], int size)
{
  /* O(size^3) intercalate counting */
  int r, c, f, ci[SIZE][SIZE+1], i;

  /* build a row and column index for each row and column */

  /* ri[x][y] records in which column element y occurs is in row x */
  /* ci[x][y] records in which row element y is in row x */

  i = 0;
  for (r = 0; r < size; r++)
    for (c = 0; c < size; c++)
      ci[c][s[r][c]] = r;
  for (r = 0; r < size - 1; r++)
    for (c = 0; c < size - 1; c++)
      for (f = c + 1; f < size; f++)
	{
	  int e = s[r][f], y = ci[c][e];
	  if (y < r)
	    continue;
	  i += s[r][c] == s[y][f];
	}
  return i;
}

struct bitmap
{
  int count;
  int values[SIZE + 1];
};

struct bitmap
testone (int s[SIZE][SIZE], int a, int b, int size)
{
  /* A function for testing if the empty cell in row a,
     column b of the Latin square in the array s
     has a forced completion. */

  struct bitmap bm;
  int t[SIZE + 1];
  int that;
  int i;

  /* clear the test and values array */
  for (i = 1; i <= size; i++)
    {
      t[i] = 0;
      bm.values[i] = 0;
    }
  bm.count = 0;

  /* for every element in the same row and column
     as the element being tested, set the
     corresponding test array element to be 1 */
  for (i = 0; i < size; i++)
    t[s[a][i]] = t[s[i][b]] = 1;

  /* keep track of all the different elements used
     in that row and column */
  for (i = 1; i <= size; i++)
    if (!t[i])
      {
	bm.values[i] = 1;
	that = i;
	bm.count++;
      }

  /* return those values */
  return bm;
}

void
print (int s[SIZE][SIZE], int size)
{
  /* print out the Latin square in the array s */
  int x, y;
  for (y = 0; y < size; y++)
    {
      for (x = 0; x < size; x++)
	{
	  printf ("%c", s[y][x] + 'a' - 1);
	}
      printf (" ");
    }
  printf ("\n");
}

void
recurse (int s[SIZE][SIZE], int cs, int size, int ic)
{
  /* recursively remove entries from a uniquely completable set,
     until only a critical set is left */
  int q, badflag = 0;
  if (cs < minsize)
    return;
  for (q = 0; q < size * size; q++)
    {
      if (s[q / size][q % size])
	{
	  int tmp = s[q / size][q % size];
	  s[q / size][q % size] = 0;
	  if (fill2 (s, cs - 1, 0, size) == 1)
	    {
	      badflag = 1;
	      break;
	    }
	  else
	    s[q / size][q % size] = tmp;
	}
    }
  if (badflag)
    recurse (s, cs - 1, size, ic);
  else
    {
      print (s, size);
      printf ("%d:%d\n", ic, cs);
      fflush (stdout);
    }
}

void
recursemin (int s[SIZE][SIZE], int cs, int size, int ic)
{
  /* same as recurse(), except remove entries where ret.count is lowest */
  int q, badflag = -1;
  int retmin = 99;
  if (cs < minsize)
    return;
  for (q = 0; q < size * size; q++)
    {
      if (s[q / size][q % size])
	{
	  int tmp = s[q / size][q % size], x;
	  struct bitmap ret;
	  s[q / size][q % size] = 0;
	  x = fill2 (s, cs - 1, 0, size);
	  ret = testone (s, q / size, q % size, size);

	  if (x == 1 && ret.count < retmin)
	    {
	      retmin = ret.count;
	      badflag = q;
	    }

	  s[q / size][q % size] = tmp;
	}
    }
  if (badflag != -1)
    {
      s[badflag / size][badflag % size] = 0;
      recursemin (s, cs - 1, size, ic);
    }
  else
    {
      {
	print (s, size);
	printf ("%d:%d\n", ic, cs);
	fflush (stdout); 
      }
    }
}

void
recursemax (int s[SIZE][SIZE], int cs, int size, int ic)
{
  /* same as recurse(), except remove entry where ret.count is biggest */
  int q, badflag = -1;
  int retmax = -1;
  if (cs < minsize)
    return;
  for (q = 0; q < size * size; q++)
    {
      if (s[q / size][q % size])
	{
	  int tmp = s[q / size][q % size], x;
	  struct bitmap ret;
	  s[q / size][q % size] = 0;
	  x = fill2 (s, cs - 1, 0, size);
	  ret = testone (s, q / size, q % size, size);

	  if (x == 1 && ret.count > retmax)
	    {
	      retmax = ret.count;
	      badflag = q;
	    }

	  s[q / size][q % size] = tmp;
	}
    }
  if (badflag != -1)
    {
      s[badflag / size][badflag % size] = 0;
      recursemax (s, cs - 1, size, ic);
    }
  else
    {
      {
	print (s, size);
	printf ("%d:%d\n", ic, cs);
	fflush (stdout); 
      }
    }
}


void
recursernd (int s[SIZE][SIZE], int cs, int size, int ic)
{
  /* same as recurse(), except remove entries randomly */
  int t[SIZE * SIZE], count = 0, q, badflag = 0;
  if (cs < minsize)
    return;
  memset (t, 0, sizeof (t));
  for (q = 0; q < size * size; q++)
    {
      if (s[q / size][q % size])
	{
	  int tmp = s[q / size][q % size], x;
	  s[q / size][q % size] = 0;
	  if (fill2 (s, cs - 1, 0, size) == 1)
	    {
	      badflag = 1;
	      t[count++] = q;
	    }
	  s[q / size][q % size] = tmp;
	}
    }
  if (badflag)
    {
      int rand1;
      rand1 = random () % count;
      s[t[rand1] / size][t[rand1] % size] = 0;

      recursernd (s, cs - 1, size, ic);
    }
  else
    {
      print (s, size);
      printf ("%d:%d\n", ic, cs);
    }
}


int
fill2 (int s[SIZE][SIZE], int level, int pos, int size)
{
  /* this function attempts to complete the 
     partial Latin square supplied in the array s.
     During the forcing process, if no element could
     possibly *be* in a position, then there
     is no point in continuing.  The function
     returns the number of possible Latin squares
     based on the given square.  */
  struct bitmap ret;
  int a = 0, b = 0, c, poss = 0;
  int min, minx, miny;
  /* go through each element in the Latin square */
  /* and test whether that element has a forced */
  /* completion. if not, move on to semi-strong */

  /* are we already finished? */
  if (level == size * size)
    return 1;

  /* try strong completion, then semistrong, then critical */
  /* now try all possibilities */
  min = size * size;
  a = b = 0;
  while (b * size + a < size * size)
    {
      if (!s[b][a])
	{
	  ret = testone (s, b, a, size);
	  if (ret.count == 0)
	    return 0;
	  if (min > ret.count)
	    {
	      minx = b;
	      miny = a;
	      min = ret.count;
	    }
	}
      a++;
      if (a == size)
	{
	  b++;
	  a = 0;
	}
    }
  b = minx;
  a = miny;
  ret = testone (s, b, a, size);
  for (c = 1; c <= size; c++)
    {
      if (ret.values[c])
	{
	  s[b][a] = c;
	  poss += fill2 (s, level + 1, b * size + a + 1, size);
	  s[b][a] = 0;
	}
      /* return immediately if > 1 */
      if (poss > 1)
	return poss;
    }
  /* if one square fails,
     forget about the rest */
  return poss;
}

void
display (int array[200][200])
{
  int i, j;
  for (i = 0; i < 200; i++)
    {
      int flag = 0;
      for (j = 0; j < 200; j++)
	{
	  if (array[i][j] > 0)
	    {
	      printf ("%d:%d:%d ", i, j, array[i][j]);
	      flag = 1;
	    }
	}
      if (flag)
	printf ("\n");
    }
  printf ("-----------------------------------------------------\n");
}


int
fill (int s[SIZE][SIZE], int level, int pos, int size)
{
  /* this function attempts to complete the
     partial Latin square supplied in the array s.
     During the forcing process, if no element could
     possibly *be* in a position, then there
     is no point in continuing.  The function
     returns the number of possible Latin squares
     based on the given square.  */
  struct bitmap ret;
  int a = 0, b = 0, c, changed, mod = 1, mod1 = 1, mod2 = 1, poss = 0, xb = 0;
  if (level == size * size)
    {
      {
	int s2[SIZE][SIZE];
	count++;
	memcpy (s2, s, sizeof (int) * SIZE * SIZE);
	if (icount (s, size) >= mininter)
	  {
	    if (recflag == 0)
	      recurse (s, size * size, size, icount (s, size));
	    if (recflag == 1)
	      recursemax (s, size * size, size, icount (s, size));
	    if (recflag == 2)
	      recursemin (s, size * size, size, icount (s, size));
	    if (recflag == 3)
	      recursernd (s, size * size, size, icount (s, size));
	  }
	memcpy (s, s2, sizeof (int) * SIZE * SIZE);
	/* if (array[icount(s,size)][0]==0)
	   {
	   array[icount(s,size)][0] = 1;
	   printf("%d %d\n",count,icount(s,size));
	   } */
	return 1;
      }
    }
  /* try strong completion, then semistrong, then critical */
  /* now try all possibilities */
  a = pos % size;
  b = pos / size;
  while (s[b][a])
    {
      a++;
      if (a == size)
	{
	  b++;
	  a = 0;
	}
    }


  ret = testone (s, b, a, size);
  for (c = 1; c <= size; c++)
    {
      if (ret.values[c])
	{
	  s[b][a] = c;
	  fill (s, level + 1, b * size + a + 1, size);
	  s[b][a] = 0;
	}
      /* return immediately if > 1 */
    }
  /* if one square fails,
     forget about the rest */
  return 0;
}
