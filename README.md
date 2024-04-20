# FGLM Algorithm

## Summary

This repository contains an implementation of the Faugère, Gianni, Lazard, and Mora (FGLM) algorithm for term order conversion of Gröbner bases. For theoretical details and applications, refer to our paper: [FGLM Algorithm Paper](FGLM%20Algorithm.pdf).

- Gröbner bases are a fundamental tool in computational algebra with applications ranging from solving systems of polynomial equations to cryptography and robotics path planning. However, the term ordering used to compute a Gröbner basis can greatly impact the computation time. Certain orderings like degree reverse lexicographic (degrevlex) are much more efficient than others like lexicographic (lex). Therefore, being able to convert a Gröbner basis from one ordering to another is very useful in practice.

## Code Implementation

Detailed worked examples are provided in the Appendix A section on pages 6-8 which correspond to the implementation provided in [FGLM.ipynb](FGLM.ipynb).  
