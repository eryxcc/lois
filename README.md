# LOIS

LOIS (Looping Over Infinite Sets) is a C++ library allowing iterating through certain
infinite sets, in finite time. 
The resulting language has an intuitive semantics, corresponding to execution of
infinitely many threads in parallel. This allows to merge the power of abstract
mathematical constructions into imperative programming. Infinite sets are internally 
represented using first order formulas over some underlying logical structure. 
To effectively handle such sets, we use and implement SMT solvers for various first
order theories. LOIS has applications in education, and in verification of infinite
state systems.

See [http://www.mimuw.edu.pl/~erykk/lois/](the official homepage of LOIS) for more
information.
