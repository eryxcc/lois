#include "../include/loisextra.h"

namespace lois {

lsetof<lvector<term>> branchset(contextptr anccontext, contextptr nowcontext) {

  lsetof<lvector<term>> X;
  
  lvector<term> v;
  contextptr e = nowcontext;
  while(e != anccontext) {
    if(!e) throw env_exception();
    for(auto w: e->var) v.push_back(term(w));
    e = e->parent;
    }
  
  X.insert(v, anccontext);

  return X;
  }

/*
// evaluate a relation (given as a graph) on x
*/

// is the given function injective over X?
}
