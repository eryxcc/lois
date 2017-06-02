#include "../include/loisinternal.h"

namespace lois {

//==== Element: element of the universe ====

// a Set (and two auxiliary classes: SimpleSet and EIterator)

extern std::shared_ptr<struct Simplifier> simplifier;

rbool simplifyWhat;

struct Simplifier : Formula {
  bool verify() { simplifyWhat->verify(); return false; }
  bool uses(vptr v) { return simplifyWhat->uses(v); }
  rbool negate() { return rbool(simplifier); }
  std::ostream& display (std::ostream &os) const { 
    return os << "[SIMPLIFY: " << simplifyWhat << "]";
    }
  rbool subst(const varsubstlist& l) const { simplifyWhat = simplifyWhat->subst(l); return rbool(simplifier); }
  rbool qsimplify2() { return rbool(simplifier); }
  int formulasize() { return 0; }
  void initDomains() { simplifyWhat->initDomains(); }
  void listFeatures(featurelist& f) { simplifyWhat->listFeatures(f); }
  };

std::shared_ptr<Simplifier> simplifier = std::make_shared<Simplifier> ();

int aggressive_simplification_tries = 0;

int simplifyVerbosity = 0;

lsetof<term> Domain::getSet() {
  auto x = newvar(this);
  varlist v;
  v.push_back(x);
  lsetof<term> res;
  term tx = term(x);
  res.elts.push_back(std::make_shared<SimpleSetOf<term>>( v, ftrue, tx));
  return res;
  }

void aggressive_simplify(varlist& selvar, rbool& selphi, contextptr& ain, std::function<bool(vptr)> a_isused, std::function<void(const varsubstlist&)> a_alpha) {
  if(selphi.isTrue()) return;
  if(!aggressive_simplification_tries) return;

  int size1 = selphi->formulasize();
  
  simplifyWhat = selphi;
  
  rbool phi ( std::make_shared<FormulaQ> (false, rbool(simplifier), selvar) );
  
  phi = outenv(phi, false, emptycontext, ain);

  if(notForInternalSolver(phi)) return;

  simplifyWhat->clearinfo();

  // cout << phi << endl;
  
  checkid++;
  phi->initDomains();

  triesleft = aggressive_simplification_tries; 
  try { phi->verify();
  } catch(unsolvable_exception) {
    triesleft = -1;
    }

  if(triesleft <= 0) {
    std::cout << "simplification failed for rbool of size " << size1 << std::endl;
    std::cout << phi << std::endl;
    std::cout << "NFI is : " << notForInternalSolver(phi) << std::endl;
    std::cout << "features are: ";

    featurelist fl;
    phi->listFeatures(fl);
    for(auto r: fl.relations) std::cout << (r->getName()); std::cout << std::endl;

    exit(1);
    return;
    }
  
  rbool orig = selphi;
  
  selphi = simplifyWhat->qsimplify();
  // cout << "simplification result: " << sel.phi << endl;

  int size3 = selphi->formulasize();
  
  if(size1 == size3) {
    // printf("Not Simplified: %d in %d tries\n", size1, aggressive_simplification_tries-triesleft);
    // cout << "orig: " << orig << endl;
    return;
    }

  if(simplifyVerbosity >= 1) 
     printf("Simplified: %d to %d in %d tries\n", size1, size3, aggressive_simplification_tries-triesleft);

  // std::cout << "FROM: " << phi << std::endl << "TO: " << selphi << std::endl;
  
  varlist var2;
  swap(var2, selvar);
  
  for(auto w: var2) 
    if(selphi->uses(w) || a_isused(w)) {
      term w1 = selphi->valueKnown(w, false);
      if(w1.p) {
        varsubstlist l = { make_pair(w, w1) };
        selphi = substitute(selphi,l);
        a_alpha(l);
        }
      else selvar.push_back(w);
      }
    
  aggressive_simplify(selvar, selphi, ain, a_isused, a_alpha);
  }

}
