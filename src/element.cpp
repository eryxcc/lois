#include "../include/loisinternal.h"

namespace lois {

//==== Element: element of the universe ====

const elem nullelem;

term& asa(const elem& x) { return as<term> (x); }

// a Set (and two auxiliary classes: SimpleSet and EIterator)

std::ostream& operator << (std::ostream & os, const SimpleSet& e) {
  os << e.a;
  bool comma = false;
  for(size_t i=0; i<e.var.size(); i++) {
    if(comma) os << sym.sfecomma; else os << sym.sfepipe;
    os << e.var[i]; os << sym.in; os << e.var[i]->dom->name;
    comma = true;
    }
  if(!e.phi.isTrue()) {
    if(comma) os << sym.sfecomma; else os << sym.sfepipe;
    os << e.phi;
    }
  return os;
  }

rbool ESet::operator == (Element& e) {
  auto e2 = dynamic_cast<ESet*>(&e);
  if(!e2) return ffalse;
  return subseteq(e2) && e2->subseteq(this);
  }

std::ostream& ESet::display (std::ostream &os) const { 
  if(elts.size() == 0) return os << sym.emptyset;
  os << sym.leftbrace;
  for(size_t i=0; i<elts.size(); i++) {
    if(i) os << sym.ssunion;
    os << elts[i];
    }
  return os << sym.rightbrace;
  }

elem ESet::alph(vptr v1, vptr v2) const { 
  auto s2 = std::make_shared<ESet> ();
  s2->elts.resize(elts.size());
  for(size_t i=0; i<elts.size(); i++) {
    const SimpleSet& e0(elts[i]);
    SimpleSet& e1(s2->elts[i]);
    e1.phi = e0.phi->alph(v1, v2);
    e1.a   = e0.a->alph(v1, v2);
    e1.var = e0.var;
    }
  return elem(s2);
  }

bool ESet::uses(vptr w) { 
  for(size_t i=0; i<elts.size(); i++) 
    if(elts[i].a->uses(w) || elts[i].phi->uses(w)) 
      return true;
  return false;
  }

const elem& EIterator::operator * () {
  if(!at.p) throw iterator_exception();
  return at;
  }
void EIterator::operator ++ () { 
  disconnectIterator();
  // if(index >= 0) 
  index++;
  connectIterator();
  }
void EIterator::connectIterator() {
  if(index >= s->elts.size() || index < 0) return;
  const SimpleSet& sel(s->elts[index]);
  at = sel.a;
  // add variables and alph!
  currentcontext = ourcontext = std::make_shared<Context> (currentcontext, sel.phi);
  for(auto w: sel.var) {
    auto v = newvar(w->dom);
    at = at->alph(w, v);
    currentcontext->phi = currentcontext->phi->alph(w, v);
    currentcontext->var.push_back(v);
    }
  // sometimes the environment is inconsistent with phi
  if(solver->solveEnv() != 2) ++(*this);
  }
void EIterator::disconnectIterator() {
  if(!at.p) return;
  if(currentcontext != ourcontext) throw iterator_exception();
  // printf("%p: disconnect\n", this);
  at = nullelem; currentcontext = ourcontext->parent;
  }

rset asSet(elem a) {
  auto sa = std::dynamic_pointer_cast<ESet> (a.p);
  if(!sa) throw as_exception();
  return sa;
  }

bool isSet(elem a) {
  return std::dynamic_pointer_cast<ESet> (a.p) ? true : false;
  }

struct EIterator begin(elem a) {
  return asSet(a)->begin(); 
  }
struct EIterator end(elem a) { 
  return asSet(a)->end();
  }

struct EIterator ESet::begin() const {
  return EIterator(this, 0);
  }
  
struct EIterator ESet::end() const {
  return EIterator(this, -1);
  }

bool EIterator::operator != (const EIterator& other) {
  if(index < 0 || index >= s->elts.size())
    return !(other.index < 0 || other.index >= other.s->elts.size());
  return index != other.index;
  }

rbool ESet::subseteq (ESet* e) {
  lbool a(ftrue); for(auto x: *this) a &= e->hasmember(x); return a;
  }

rbool ESet::hasmember(elem e) {
  lbool a(ffalse); for(auto x: *this) a |= (e == x); return a;
  }

extern std::shared_ptr<struct Simplifier> simplifier;

rbool simplifyWhat;

struct Simplifier : Formula {
  bool verify() { simplifyWhat->verify(); return false; }
  bool uses(vptr v) { return simplifyWhat->uses(v); }
  rbool negate() { return rbool(simplifier); }
  std::ostream& display (std::ostream &os) const { 
    return os << "[SIMPLIFY: " << simplifyWhat << "]";
    }
  rbool alph(vptr v1, vptr v2) const { simplifyWhat = simplifyWhat->alph(v1, v2); return rbool(simplifier); }
  rbool qsimplify2() { return rbool(simplifier); }
  int formulasize() { return 0; }
  void initDomains() { simplifyWhat->initDomains(); }
  void listFeatures(featurelist& f) { simplifyWhat->listFeatures(f); }
  };

std::shared_ptr<Simplifier> simplifier = std::make_shared<Simplifier> ();

int aggressive_simplification_tries = 0;

int simplifyVerbosity = 0;

void aggressive_simplify(varlist& selvar, rbool& selphi, elem& a, contextptr& ain) {
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
    if(selphi->uses(w) || a->uses(w)) {
      vptr w1 = selphi->valueKnown(w, false);
      if(w1) selphi = selphi->alph(w, w1), a = a->alph(w, w1);
      else selvar.push_back(w);
      }
    
  aggressive_simplify(selvar, selphi, a, ain);
  }


void ESet::insert(elem a, contextptr ain, contextptr nowenv) {
  
  // cout << "insert: " << a << " env=" << currentcontext << " setenv="<<ain << endl;
  
  varlist selvar;
  rbool selphi = ftrue;
  
  varlist to_quantify;
  
  contextptr ainlocal = nowenv;  

  while(ainlocal != ain) {
    if(!ainlocal) throw env_exception();
    rbool aphi = ainlocal->phi;
    selphi = selphi && aphi;
    for(auto w: ainlocal->var) 
      if(selphi->uses(w) || a->uses(w)) {
        vptr w1 = selphi->valueKnown(w, false);
        if(w1) selphi = selphi->alph(w, w1), a = a->alph(w, w1);
        else if(a->uses(w)) 
          selvar.push_back(w);
        else to_quantify.push_back(w);
        }
    ainlocal = ainlocal->parent;
    }
  if(to_quantify.size()) {
    auto phi = std::make_shared<FormulaQ> (false, selphi, to_quantify);
    selphi = phi->simplify(phi);
    }
#ifdef AGGSYM
  aggressive_simplify(selvar, selphi, a, ain);
#endif

  elts.emplace_back(selvar, selphi, a);
  }

rset newSet(elem x) { 
  rset S = newSet(); S->insert(x, currentcontext); return S; 
  }

rset newSet(elem x, elem y) { 
  rset S = newSet(); S->insert(x, currentcontext); S->insert(y, currentcontext); return S; 
  }

rset newSet(std::initializer_list<elem> l) {
  rset S = newSet(); for(elem it: l) S->insert(it, currentcontext); return S;
  }

rsetof<term> Domain::getSet() {
  auto x = newvar(this);
  std::shared_ptr<ESet> D = std::make_shared<ESet> ();
  D->elts.resize(1);
  SimpleSet& e(D->elts[0]);
  e.phi = ftrue;
  e.a = x;
  e.var.push_back(x);
  rset reps = rset(D);
  return rsetof<term> (reps);
  }

rset emptyset;

}
