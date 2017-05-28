// Note that this example is parsed by
// parsetutorial.cpp (not included in the package) in order to include
// its fragments and output into the paper easily.

// This file has the following parts:
// - the "Tutorial" from the paper.
// - more examples from the paper. Output is also included in some cases.
// - example LOIS functions from the paper (no output included).

#define SECTION(x) cout << "SECTION " << x << endl;
#define FUNCTION(x) 


FUNCTION("initialization")

#include "../include/loisextra.h"

using namespace std;
using namespace lois;

#undef FILTERMAP
#undef FILTER
#undef MAP

#define FILTERMAP(x, set, condition, replace)                     \
  filtermap(set,                                                   \
    (std::function<rbool(const term&)>) ([=] (const __typeof(set.elts[0].a)& x) { return condition; },)     \
    (std::function<lpair<term,term>(const term&)>) ([=] (const __typeof(set.elts[0].a)& x) { return replace; })         \
    );

template<class T, class F, class M, class U> lsetof<U> 
  filtermap(const lsetof<T>& myset,
    std::function<rbool(const T&)> filter,
    std::function<U(const T&)> mapper
    ) {
    lsetof<U> result;
    for(const T& t: myset) If(filter(t)) result += mapper(t);
    return result;
    }

#define FILTER(x, set, condition) FILTERMAP(x, set, condition, x)
#define MAP(x, set, replace) FILTERMAP(x, set, ftrue, replace)

int main() {
  initLois();
  sym.useUnicode();
  Domain dA("A");  
  lsetof<term> A = dA.getSet();
  
  // vector<string> v;
  // typedef __typeof(v) vtype;
  // vtype::value_type vs = "hello";
  
  // typeof(A)::element elt = nullterm;
  
  std::function<rbool(const term& x)> filterer = 
    [&] (const term& x) { return true; };
  std::function<lpair<term,term> (const term& x)> mapper = 
    [&] (const term& x) { return make_lpair(x,x); };
  
  cout << 
    filtermap(A, 
      (std::function<rbool(const term&)>) ([&] (const term& x) { return true; }),
      (std::function<lpair<term,term> (const term&)>) ([&] (const term& x) { return make_lpair(x,x); })
      ) << endl;
  
  //MAP(x, A, make_lpair(x,x)) << endl;

  return 0;
  }
