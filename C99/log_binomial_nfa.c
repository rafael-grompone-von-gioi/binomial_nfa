/*----------------------------------------------------------------------------

  Copyright (c) 2007-2021 rafael grompone von gioi <grompone@gmail.com>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

  ----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/** ln(10) */
#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif /* !M_LN10 */

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/*----------------------------------------------------------------------------*/
/* Fatal error, print a message to standard-error output and exit.
 */
static void error(char * msg)
{
  fprintf(stderr,"error: %s\n",msg);
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/* Doubles relative error factor
 */
#define RELATIVE_ERROR_FACTOR 100.0

/*----------------------------------------------------------------------------*/
/* Compare doubles by relative error.

   The resulting rounding error after floating point computations
   depend on the specific operations done. The same number computed by
   different algorithms could present different rounding errors. For a
   useful comparison, an estimation of the relative rounding error
   should be considered and compared to a factor times EPS. The factor
   should be related to the cumulated rounding error in the chain of
   computation. Here, as a simplification, a fixed factor is used.
 */
static int double_equal(double a, double b)
{
  double abs_diff,aa,bb,abs_max;

  /* trivial case */
  if( a == b ) return TRUE;

  abs_diff = fabs(a-b);
  aa = fabs(a);
  bb = fabs(b);
  abs_max = aa > bb ? aa : bb;

  /* DBL_MIN is the smallest normalized number, thus, the smallest
     number whose relative error is bounded by DBL_EPSILON. For
     smaller numbers, the same quantization steps as for DBL_MIN
     are used. Then, for smaller numbers, a meaningful "relative"
     error should be computed by dividing the difference by DBL_MIN. */
  if( abs_max < DBL_MIN ) abs_max = DBL_MIN;

  /* equal if relative error <= factor x eps */
  return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * DBL_EPSILON);
}

/*----------------------------------------------------------------------------*/
/** Size of the table to store already computed inverse values.
 */
#define TABSIZE 100000

/*----------------------------------------------------------------------------*/
/** Computes log10(NFA).

    NFA stands for Number of False Alarms:
    @f[
        \mathrm{NFA} = NT \cdot B(n,k,p)
    @f]

    - NT       - number of tests
    - B(n,k,p) - tail of binomial distribution with parameters n,k and p:
    @f[
        B(n,k,p) = \sum_{j=k}^n
                   \left(\begin{array}{c}n\\j\end{array}\right)
                   p^{j} (1-p)^{n-j}
    @f]

    The value log10(NFA) is equivalent but more intuitive than NFA:
    -   1 corresponds to 10 mean false alarms
    -  -0 corresponds to 1 mean false alarm
    -  -1 corresponds to 0.1 mean false alarms
    -  -2 corresponds to 0.01 mean false alarms
    -  ...

    Used this way, the smaller the value, better the detection,
    and a logarithmic scale is used.

    @param n,k,p binomial parameters.
    @param logNT logarithm of Number of Tests

    The computation is based in the gamma function by the following
    relation:
    @f[
        \left(\begin{array}{c}n\\k\end{array}\right)
        = \frac{ \Gamma(n+1) }{ \Gamma(k+1) \cdot \Gamma(n-k+1) }.
    @f]
    We use efficient algorithms to compute the logarithm of
    the gamma function.

    To make the computation faster, not all the sum is computed, part
    of the terms are neglected based on a bound to the error obtained
    (an error of 10% in the result is accepted).
 */
double log_binomial_nfa(int n, int k, double p, double logNT)
{
  static double inv[TABSIZE];   /* table to keep computed inverse values */
  double tolerance = 0.1;       /* an error of 10% in the result is accepted */
  double log1term,term,bin_term,mult_term,bin_tail,err;
  double p_term = p / (1.0-p);
  int i;

  if( n<0 || k<0 || k>n || p<0.0 || p>1.0 )
    error("wrong n, k or p values in nfa()");

  if( n==0 || k==0 ) return logNT;
  if( n==k ) return logNT + (double)n * log10(p);

  log1term = lgamma((double)n+1.0) - lgamma((double)k+1.0)
           - lgamma((double)(n-k)+1.0)
           + (double)k * log(p) + (double)(n-k) * log(1.0-p);

  term = exp(log1term);
  if( term == 0.0 )                        /* the first term is almost zero */
    {
      if( (double)k > (double)n * p )      /* at begining or end of the tail? */
        return log1term / M_LN10 + logNT;  /* end: use just the first term */
      else
        return logNT;                      /* begin: the tail is roughly 1 */
    }

  bin_tail = term;
  for(i=k+1;i<=n;i++)
    {
      bin_term = (double)(n-i+1) * ( i<TABSIZE ?
                   (inv[i] ? inv[i] : (inv[i]=1.0/(double)i)) : 1.0/(double)i );
      mult_term = bin_term * p_term;
      term *= mult_term;
      bin_tail += term;
      if(bin_term<1.0)
        {
          /* when bin_term<1 then mult_term_j<mult_term_i for j>i.
             then, the error on the binomial tail when truncated at
             the i term can be bounded by a geometric serie of form
             term_i * sum mult_term_i^j.                            */
          err = term * ( ( 1.0 - pow( mult_term, (double)(n-i+1) ) ) /
                         (1.0-mult_term) - 1.0 );

          /* one wants an error at most of tolerance*final_result, or:
             tolerance * abs(-log10(bin_tail)-logNT).
             now, the error that can be accepted on bin_tail is
             given by tolerance*final_result divided by the derivative
             of -log10(x) when x=bin_tail. that is:
             tolerance * abs(-log10(bin_tail)-logNT) / (1/bin_tail)
             finally, we truncate the tail if the error is less than:
             tolerance * abs(-log10(bin_tail)-logNT) * bin_tail        */
          if( err < tolerance * fabs(-log10(bin_tail)-logNT) * bin_tail ) break;
        }
    }
  return log10(bin_tail) + logNT;
}

/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv)
{
  int n,k;
  double p,logNT,logNFA;

  if( argc != 5 ) error("use: log_binomial_nfa log10(NT) n k p");
  logNT = atof(argv[1]);
  n     = atoi(argv[2]);
  k     = atoi(argv[3]);
  p     = atof(argv[4]);

  logNFA = log_binomial_nfa(n,k,p,logNT);

  printf("logNT %g n %d k %d p %g logNFA %g NFA 10^%g\n",
         logNT, n, k, p, logNFA, logNFA );
}
/*----------------------------------------------------------------------------*/
