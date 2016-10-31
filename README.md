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

See [the official homepage of LOIS](http://www.mimuw.edu.pl/~erykk/lois/) for more
information.

# Subdirectories

The subdirectory `includes` contains the headers of the LOIS library.

The subdirectory `src` contains sources of the LOIS library itself.

The subdirectory `obj` contains generated object files.

The subdirectory `tests` contains programs using LOIS -- this includes the tutorial
from the technical documentation, automatic testing program (`autotest`),
a program comparing solvers (`soltest`), and some applications (e.g. `learning`).

The subdirectory `out` contains output from the executables.

The subdirectory `bin` contains generated binaries (both liblois.a and test executables).

