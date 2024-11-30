# Code for finding smallest and largest critical sets in Latin squares (scs and lcs)

The concept of "critical sets in Latin squares" was invented by John Nelder in 1977, in a CSIRO publication.

The smallest and largest critical sets can be written as functions of the order of the Latin square "n".
So, although Nelder did not give a conjecture about "scs(n)" and "lcs(n)", we may guess that his conjecture was
scs(n) = floor(n^2/4) and lcs(n) = (n^2-n)/2.

Stinson and van Rees wrote a 1982 article on "Some Large Critical Sets" which showed critical sets of order 9 and size 39, and order 10 and size 55.
In my PhD (1998-2001) I found critical sets of order 9/size 44, order 10/size 57, order 11/size 70, order 12/size 90, order 14/size 118 and in 2002 I found a critical sets of order 9 and size 45.

The question is can we improve on that "lcs" bound? I would expect my conjecture from Bean and Mahmoodian is true -- lcs(n) ≤ n^2 - 3^(log_2 n)

So this would imply 
```
25 <= lcs(7) <= 27
lcs(8) = 37
45 <= lcs(9) <= 48
57 <= lcs(10) <= 61
70 <= lcs(11) <= 76
90 <= lcs(12) <= 92
95 <= lcs(13) <= 110 (30 November 2024, see below)
118 <= lcs(14) <= 130
146 <= lcs(15) <= 151
lcs(16) = 175
```
etc. (Does OEIS A080572 give a construction for a critical set of order 2^n+1?)

I would guess that lcs(7) = 25 - interesting question is how many of the 147 main classes of order 7 Latin squares have a critical set of size 25 -- in my thesis, I found 113 main classes with such critical sets.

So what can we do for lcs(9), lcs(10), lcs(11) etc?

## References

Richard Bean and E. S. Mahmoodian, [A new bound on the size of the largest critical set in a Latin square](https://arxiv.org/abs/math/0107159) arXiv:math/0107159 [math.CO], 2001.

C. Fu, H. Fu and W. Liao, A new construction for a critical set in special Latin squares, Proceedings of the Twenty-sixth Southeastern International Conference on Combinatorics, Graph Theory and Computing (Boca Raton, Florida, 1995), Congressus Numerantium, Vol. 110 (1995), pp. 161-166.

D. R. Stinson and G. H. J. van Rees, Some large critical sets, Proceedings of the Eleventh Manitoba Conference on Numerical Mathematics and Computing (Winnipeg, Manitoba, 1981), Congressus Numerantium, Vol. 34 (1982), pp. 441-456. 

J. Nelder, Critical sets in Latin squares, CSIRO Division of Math. and Stats. Newsletter, Vol. 38 (1977), p. 4.

## The On-Line Encyclopedia of Integer Sequences (OEIS) critical set sequences

[OEIS A002620 -- Quarter-squares: a(n) = floor(n/2)*ceiling(n/2). Equivalently, a(n) = floor(n^2/4)](http://oeis.org/A002620)

[OEIS A063437 -- Cardinality of largest critical set in any Latin square of order n](http://oeis.org/A063437)

[OEIS A080572 -- Number of ordered pairs (i,j), 0 <= i,j < n, for which (i & j) is nonzero, where & is the bitwise AND operator](http://oeis.org/A080572)

## Appendix

Order 13 size 95 critical set, completion has 125 intercalates

```
.............
.a.c.e.g..i.j
..ab..ef....l
.cba.gfe..jki
....abcd.iljk
.e......ab..g
..e..l.kc...f
.gf...lj..dce
.....abidhgfc
..l.bjdchefga
...gkcalefbhd
.l.fdkibgceah
.kilcdjafgheb
```
