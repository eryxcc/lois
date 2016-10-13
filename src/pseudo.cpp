#include "../include/loisextra.h"

namespace lois {

rset branchset(contextptr anccontext, contextptr nowcontext) {

  rset X = newSet();
  
  eltuple v;
  contextptr e = nowcontext;
  while(e != anccontext) {
    if(!e) throw env_exception();
    for(auto w: e->var) v.push_back(w);
    e = e->parent;
    }
  
  X->insert(v, anccontext);

  return X;
  }

lnum<int> cardinality(const lset& X) {
  lnum<int> answer;
  for(auto a: X) If(answer.isUndefined()) {
    auto butone = cardinality(X &~ newSet(a));
    for(int x: butone.singleton) answer.singleton += (x+1);
    }
  If(answer.isUndefined()) answer.singleton += 0;
  return answer;
  }

// evaluate a relation (given as a graph) on x
lelem eval(const lset& f, elem x) {
  lset ret;
  for(elpair p: lsetof<elpair> (f))
    If(x == p.first) ret += p.second;
  return extract(ret);
  }

// is the given function injective over X?
rbool isInjective(const lset& f, const lset& X) {
  lbool ret(ftrue);
  for (auto x:X) for(auto y:X) If(x!=y) {
    ret &= (eval(f,x) != eval(f,y));
    }
  return ret;
  }

}
