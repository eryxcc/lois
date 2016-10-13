#include "../include/loisinternal.h"

namespace lois {

// *** environment

const contextptr emptycontext;
contextptr currentcontext = emptycontext;


// The environment says which variables are currently iterated, 
// and under what constraints. The environment changes when 
// iterating over sets and enabling conditions.

// (Note: environment always has a single parent, so it is only allowed to 
// disable the last iterator/condition)

std::ostream& Context::display(std::ostream& os, contextptr upto) {
  if(this == &(*upto)) os << "|";
  else if(this == NULL) throw env_exception();
  else {
    for(auto w: var) os << "|" << w << sym.in << w->dom->name;
    if(!phi.isTrue()) os << "|" << phi;
    parent->display(os, upto);
    }
  return os;
  }

// the current environment

// std::ostream& operator << (std::ostream& os, contextptr e) { return env->display(os, e); }

// lift the rbool 'phi' from 'env' to 'anccontext' by quantifying over all the variables
rbool outenv(rbool phi, bool universal, contextptr anccontext, contextptr nowcontext) {
  auto qquery = std::make_shared<FormulaQ> (universal, phi);
  contextptr e = nowcontext;
  while(e != anccontext) {
    if(!e) throw env_exception();
    if(universal)
      qquery->right = (!e->phi) || qquery->right;
    else
      qquery->right = (e->phi) && qquery->right;
    for(auto w: e->var)
      if(qquery->right->uses(w)) qquery->var.push_back(w);
    e = e->parent;
    }
  
  qquery->setfv();
  
  return qquery->simplify(qquery);
  }

//==== condition ====

// a helper class to enable/disable constraints to the environment

void ArbCondition::disable() {
  if(ourcontext) {
    if(currentcontext == ourcontext) {
      currentcontext = ourcontext->parent, ourcontext = emptycontext;
      }
    else {
      throw condition_exception();
      }
    }
  }

bool ArbCondition::enable(rbool psi) { 
  disable();
  currentcontext = ourcontext = std::make_shared<Context> (currentcontext, psi);
  int res = solver->solveEnv();
  if(res == 1) { throw unsolvable_exception(); }
  return res == 2;
  }

// the usual case: a binary condition

struct Condition : ArbCondition {
  rbool phi;
  Condition(rbool _phi) : phi(_phi) {}
  operator bool() { return enable(phi); }
  bool operator !() { return enable(!phi); }
  };

// using formulae as iterators!

void CondIterator::go() {
  if(state == ciIfThen)
    if(!enable(phi)) state = ciEnd;
  if(state == ciIfThenElse)
    if(!enable(phi)) state = ciElse;
  if(state == ciElse)
    if(!enable(!phi)) state = ciEnd;
  if(state == ciWhile) {
    if(!enable(phi)) { disable(); state = ciWhileFailed; }
    }
  }
void CondIterator::operator ++ () { 
  if(state == ciIfThenElse) state = ciElse;
  else state = ciEnd;
  go();
  }

CondIterator begin(const formulamode& fm) {
  CondIterator it; it.state = fm.state; it.phi = fm.phi; 
  it.go(); 
  return it;
  }
CondIterator end(const formulamode& fm) { 
  CondIterator it; it.state = ciEnd; return it;
  }

BreakIterator begin(breakableIteratorType& b) { return BreakIterator(true); }
BreakIterator end(breakableIteratorType& b) { return BreakIterator(false); }

breakableIteratorType breakableIterator;

}