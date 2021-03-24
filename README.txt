Computation of a NFA related to a binomial distribution
=======================================================

Version 0.1 - March 24, 2021
by Rafael Grompone von Gioi <grompone@gmail.com>

This is a code to compute a binomial NFA as defined in

  "From Gestalt Theory to Image Analysis, a Probabilistic Approach" by
  A. Desolneux, L. Moisan, and J.M. Morel, Springer}, 2008.

There are two versions of the code in C language, one following ANSI C89
standard and the other the ANSI C99 standard.

To compile just execute 'make' on the corresponding directories, C89 or C99.

A typical execution is:

  $ ./log_binomial_nfa 10 100 80 0.1
  logNT 10 n 100 k 80 p 0.1 logNFA -50.1742 NFA 10^-50.1742
