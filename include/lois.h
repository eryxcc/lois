// LOIS: Looping Over Infinite Sets

#ifndef _lois_h_
#define _lois_h_

#define AGGSYM

#include <vector>
#include <string>
#include <ostream>
#include <iostream>
#include <memory>

namespace lois {

// --- LOIS uses the following internally for pretty indenting.

extern int indent;

struct indenter { int i; int doit; };

extern indenter iindent, ido, ieol, iunindent;

std::ostream& operator << (std::ostream& os, indenter a);

struct autoindenter { autoindenter() { indent += 2; } ~autoindenter() { indent -= 2; } };

// --- Error handling: exceptions which are thrown by LOIS for
// incorrect programs.

struct loisexception : std::exception {
  virtual const char* what() const noexcept { return "LOIS exception"; };
  };

struct noterms_exception : loisexception {};
struct domain_exception : loisexception {};
struct as_exception : loisexception {};
struct iterator_exception : loisexception {};
struct env_exception : loisexception {};
struct condition_exception : loisexception {};
struct unsatisfiable_exception : loisexception {};
struct unsolvable_exception : loisexception {};
struct subdomain_exception : loisexception {};
struct inconsistency_exception : loisexception {};
struct assignment_exception : loisexception {};

struct other_exception : loisexception {
  std::string s;
  other_exception(std::string _s) : s(_s) {}
  virtual const char* what() const noexcept { return s.c_str(); };
  ~other_exception() noexcept(true) {}
  };

// -- query languages allowed --
// outlan defines the language to use when displaying formulas.
// Additionally, L0, SMT, SMT_INC, CVC3, and SPASS constants are also used 
// in function Relation::worksWith to tell whether the given relation works
// with the given solver.

enum OutputLanguage { L0, CVC3, SMT, SPASS, SPASS_SYMBOLS, SPASS_AXIOMS, SMT_INC };

extern OutputLanguage outlan;

// --- the structure 'sym' lists the symbols which are used for rendering of
// various mathematical notions. One can switch between using the unicode, LaTeX,
// or ASCII symbols (and potentially implement new symbol sets).

struct symbol2 { virtual std::string asString() = 0; };

struct symbolvariable : symbol2 { 
  std::shared_ptr<std::string> s; 
  void set(std::string _s) { *s = _s; }
  symbolvariable() { s = std::make_shared<std::string> ("[undefined symbol]"); }
  symbolvariable(const std::string& _s) { s = std::make_shared<std::string> (_s); }
  symbolvariable(const char *_s) { s = std::make_shared<std::string> (_s); }
  virtual std::string asString() { return *s; } 
  };

struct symbol {
  std::shared_ptr<symbol2> sym;
  std::string asString() { return sym->asString(); }
  symbol(symbolvariable v) { sym = std::make_shared<symbolvariable> (v); }
  symbol(const std::string& s) { sym = std::make_shared<symbolvariable> (s); }
  symbol(const char *s) { sym = std::make_shared<symbolvariable> (s); }
  };

extern bool debuglois;

inline std::ostream& operator << (std::ostream& os, symbol a) { return os << a.asString(); }

struct symbols {
  symbolvariable
    // basic logic
    exists, forall, _and, _or, eq, neq, in, 
    // basic sets
    emptyset, ssunion, leftbrace, rightbrace, sfepipe, sfecomma, pseudo,
    // relational and functional symbols
    leq, geq, greater, less, max, min, plus, times, minus, divide,
    edge, noedge, arrow, notarrow;
  
  void useUnicode();
  void useLaTeX();
  void useASCII();
  };

extern symbols sym;

// --- LOIS uses the shared pointers to handle its objects easily.

struct Formula;
struct term;
struct Subdomain;
struct rbool;

template<class T> struct lsetof;

typedef std::shared_ptr<struct Variable> vptr; 
typedef std::shared_ptr<struct SubRelation> subrelptr;
typedef std::shared_ptr<struct TermVariable> tvptr;
typedef std::vector<vptr> varlist;
typedef std::shared_ptr<struct Context> contextptr;

typedef std::pair<vptr, term> varsubst;
typedef std::vector<varsubst> varsubstlist;

// --- rbool, or the first order formulae

extern rbool ftrue, ffalse;

struct rbool {
  std::shared_ptr<Formula> p;
  rbool(std::shared_ptr<Formula> f) : p(f) {}
  rbool(bool b) { p = (b?ftrue:ffalse).p; }
  rbool() {}
  Formula* operator -> () const { return &(*p); }
  bool isTrue() const { return p == ftrue.p; }
  bool isFalse() const { return p == ffalse.p; }
  };

// binary operations: negation, and, or
rbool operator ! (rbool a);
rbool operator && (rbool a, rbool b);
rbool operator || (rbool a, rbool b);

// display the formula
std::ostream& operator << (std::ostream& os, rbool a);

// equality and inequality of elements of the universe, given by the terms
rbool operator == (const term& a, const term& b);
rbool operator != (const term& a, const term& b);

// make a quantified formula (both existentially and universally)
// va - variables, r - inner formula, t - true if universal
rbool makeq(bool t, rbool r, const varlist& va);

// -- domains & relations --

struct Domain {
  // name of the domain, to be displayed when printing contexts
  std::string name;
  // private
  int curcheck;
  // get the domain as a lsetof<term>
  struct lsetof<term> getSet();

  Domain(std::string s) : name(s), curcheck(-1) {}
  };

// Relation objects are created to add extra structure (depending on the 
// subclass) to the domain a Relation can have multiple relational and
// functional symbols, identified by id's

struct Relation {
  virtual subrelptr genSub() = 0;
  // display the symbol of the relation (this->id)
  virtual std::ostream& display(std::ostream& os, int id, const term& t1, const term& t2) = 0;
  // display the symbol of application of the function (this->id)
  virtual std::ostream& displaytermbin(std::ostream& os, const struct TermBinary* t) 
    { throw noterms_exception(); return os; }
  virtual std::ostream& displayconst(std::ostream& os, const struct TermConst* t) 
    { throw noterms_exception(); return os; }
  // give the id of the negation of this->id
  virtual int negate(int id) = 0;
  // make a rbool representing the binary relation (rel->id) on a and b
  virtual rbool binform(int id, const term& a, const term& b) = 0;
  // is this relation available in the given solver?
  virtual bool worksWith(OutputLanguage lan) { return false; }
  // get the name of this relation
  virtual std::string getName() { return "[unknown relation]"; }
  virtual int checktermconst(int id) { throw noterms_exception(); }
  // get the isomorphism formula
  rbool isomorphic(std::vector<std::pair<vptr, vptr>> mapping) { throw loisexception(); }
  };

// make a rbool representing the binary relation (rel->id) on terms a and b
rbool makebinary(struct Relation *rel, int id, const term& a, const term& b);

// -- terms --
const vptr nullvptr;

// prefix used for the auto-generated variable names
extern std::string autoprefix;

// term
struct term {
  tvptr p;
  term(tvptr v);
  term() { p = NULL; }
  void display(std::ostream& os);
  vptr asVar() const;
  };

term substitute(const term&, const varsubstlist&);
rbool substitute(const rbool& phi, const varsubstlist&);

template<class T> T substitute(const T& x, vptr var, term& val) {
  varsubstlist l(1, make_pair(var, val));
  return substitute(x, l);
  }

extern const term nullterm;

// syntactic equality of terms (not yet implemented fully!)
inline bool eqterm(const term& a, const term& b) {
  return a.asVar() && a.asVar() == b.asVar();
  }

// terms are internally shared pointers to a TermVariable (which could be
// either a Term (TermBinary, TermConst) or a Variable)

struct TermVariable {
  int curcheck;
  virtual int getValue() = 0;
  virtual Domain* getDom() = 0;
  virtual std::ostream& display (std::ostream& os) const = 0;
  virtual term subst(const term& ths, const varsubstlist& l) const = 0;
  virtual bool uses(vptr v) = 0;
  virtual vptr subFindv(tvptr) = 0;
  virtual void listFeatures(struct featurelist& f) = 0;
  };

struct Variable : TermVariable {
  Domain *dom;
  int curcheck;
  std::shared_ptr<Subdomain> sub;
  vptr sublink;
  mutable int id;
  int value;
  std::string name;
  Variable(Domain *d) : dom(d) { id = 0; curcheck=-1; }
  Variable(Domain *d, const std::string&  _n) : dom(d) { name = _n; curcheck=-1; }
  virtual std::ostream& display (std::ostream& os) const;
  int getValue() { return value; }
  Domain *getDom() { return dom; }
  virtual term subst(const term& ths, const varsubstlist& l) const { 
    for(auto& p: l)
      if(ths.p == p.first) return p.second;
    return ths;
    }
  virtual bool uses(vptr v) { return &*v == this; }
  virtual vptr subFindv(tvptr);
  virtual void listFeatures(struct featurelist& f);
  };

struct TermBinary : TermVariable {
  Relation *r;
  int op;
  term left, right;
  vptr sublink;
  TermBinary(Relation *_r, int _op, const term& _left, const term& _right) : 
    r(_r), left(_left), right(_right), op(_op) {
    curcheck = -1;
    if(left.p->getDom() != right.p->getDom()) throw domain_exception();
    }
  Domain *getDom() { return left.p->getDom(); }
  virtual std::ostream& display (std::ostream& os) const { 
    return r->displaytermbin(os, this);
    }
  int getValue(); 
  virtual term subst(const term& ths, const varsubstlist& l) const {
    return term(std::make_shared<TermBinary>(r, op, substitute(left,l), substitute(right,l)));
    }
  virtual bool uses(vptr v) { return left.p->uses(v) || right.p->uses(v); }
  virtual vptr subFindv(tvptr);
  virtual void listFeatures(struct featurelist& f);
  };

struct TermConst : TermVariable {
  Relation *r;
  vptr v; // we need a variable here for technical reasons
  int op;
  TermConst(Relation *_r, int _op, Domain *d) : 
    r(_r), op(_op) {
    curcheck = -1; v = std::make_shared<Variable> (d);
    }
  Domain *getDom() { return v->dom; }
  virtual std::ostream& display (std::ostream& os) const { 
    return r->displayconst(os, this);
    }
  int getValue(); 
  virtual term subst(const term& ths, const varsubstlist& l) const {
    return ths;
    // return term(std::make_shared<TermConst>(r, op, v->dom));
    }
  virtual bool uses(vptr w) { return v==w; }
  virtual vptr subFindv(tvptr);
  virtual void listFeatures(struct featurelist& f);
  };

std::ostream& operator << (std::ostream& os, vptr a);

// create a new variable in the given domain
inline vptr newvar(Domain *d) { 
  return std::make_shared<Variable> (d);
  }

// create a new variable in the given domain, with the given name
inline vptr newvar(Domain *d, const std::string& s) {
  return std::make_shared<Variable> (d, s); 
  }

// isused: is the variable v used in the second parameter?
bool isused(vptr v, vptr w);
bool isused(vptr v, rbool phi);

inline std::ostream& operator << (std::ostream& os, const term& a) { 
  if(!a.p) os << "[null term]";
  else a.p->display(os); 
  return os;
  }

// -- context, or the LOIS stack --
// the stack is implemented as a linked list, allowing the inner contexts
// to simply point to a part of the stack

// nulleptr is the initial (empty) stack
extern const contextptr emptycontext;

// currentcontext is the pointer to the most recent varset on the stack 
extern contextptr currentcontext;

struct Context {
  // the rest of the stack
  contextptr parent;
  // constraint (only one rbool is necessary)
  rbool phi;
  // new variables here
  varlist var;
  // constructor and display
  Context(contextptr par, rbool fi, varlist v = varlist()) : parent(par), phi(fi), var(v) {}
  };

std::ostream& displaycontext(std::ostream& os, contextptr what, contextptr upto);

// display the varset difference, FROM (current) currentcontext TO (parameter) e
// (for example, you should do "<< emptycontext" to view the whole current stack)
inline std::ostream& operator << (std::ostream& os, contextptr e) { 
  return displaycontext(os, currentcontext, e); 
  }

// display the varset difference, FROM e2.first (e.g. currentcontext) TO e2.second (e.g. emptycontext)
inline std::ostream& operator << (std::ostream& os, std::pair<contextptr, contextptr> e2) { 
  return displaycontext(os, e2.first, e2.second); 
  }

// lift the formula 'phi' from 'currentcontext' to 'anccontext' by applying all the
// new parts of the context which have been pushed between them (i.e., quantifying
// over the variables which are consistent with the formulae, and adding the 
// constraints as appropritae -- quantification is universal or existential, based
// on the 'universal' choice)

rbool outenv(rbool phi = ftrue, bool universal = false, contextptr anccontext = emptycontext, contextptr nowcontext = currentcontext);

// does the given object use the new variables from 'nowcontext' to 'anccontext'
template<class T> bool usesenv(const T& obj, contextptr anccontext = emptycontext, contextptr nowcontext = currentcontext) {
  contextptr e = nowcontext;
  while(e != anccontext) {
    if(!e) throw env_exception();
    for(auto w: e->var) if(isused(w, obj)) return true;
    e = e->parent;
    }
  return false;
  }

// -- solvers --

// solver: verify whether the given formula is satisfiable

// 0 = not satisfiable
// 1 = unknown
// 2 = satisfiable

struct Solver {
  // solve a formula
  virtual int solve(rbool phi) = 0;

  // solve the current environment.
  // the default implementation is 
  virtual int solveEnv() { return solve(outenv()); }
  // but better implementations can be provided to 
  // take advantage of incremental solving.
  };

extern std::string smtLogic;
std::string smtsort(); // Real or Int

typedef std::shared_ptr<Solver> solveptr;
extern solveptr solver;

// dummy: always return the same value
solveptr solverDummy(int i);

// crash: cannot solve, so crash LOIS
solveptr solverCrash();

// stack of solvers
solveptr solverStack(solveptr s, solveptr fallback);

inline solveptr operator || (solveptr a, solveptr b) { return solverStack(a, b); }

// solve only the trivial cases:
solveptr solverBasic();

// LOIS's default solver. t is the time allowed for solving
// (measured in the number of tries), v is verbosity.
solveptr solverExhaustive(int t, bool v);

// verbose:: output the formula first
solveptr solverVerbose(std::string m, solveptr s);

// named: for each formula, output the solver's name, result, and time
solveptr solverNamed(std::string n, solveptr s);

// diagnostic: compare solvers and check for inconsistencies
solveptr solverCompare(std::initializer_list<solveptr> p);

// external solvers:
solveptr solverSMT();
solveptr solverSMT(std::string t);
solveptr solverCVC();
solveptr solverCVC(std::string t);
solveptr solverSPASS();
solveptr solverSPASS(std::string t);

solveptr solverIncremental(std::string t);

// use the default solver sequence (i - the number of tries from which output 
// diagnostic information, j - the number of tries from which crash)

void useDefaultSolver(int i = 500, int j = 2000000000);

// a helper class to enable/disable constraints in the environment
// enable = push a rbool constraint on the stack, disable = pop it

struct ArbCondition {
  contextptr ourcontext;
  ArbCondition()  { ourcontext = emptycontext; }
  void disable();
  bool enable(rbool psi);
  template<class T> bool enable(T vars, rbool psi) {
    disable();
    currentcontext = ourcontext = std::make_shared<Context> (currentcontext, psi);
    for(auto v: vars) currentcontext->var.push_back(v);
    int res = solver->solveEnv();
    if(res == 1) throw unsolvable_exception();
    return res == 2;
    }
  virtual ~ArbCondition() { disable(); }
  };

// -- macros implementing the If and While construct
// this is extremely technical

enum ciState { ciIfThen, ciIfThenElse, ciWhile, ciElse, ciEnd, ciWhileFailed };

extern struct breakableIteratorType {} breakableIterator;

struct BreakIterator {
  bool b;
  BreakIterator(bool _b) : b(_b) {}
  bool& operator* () { return b; }
  bool operator != (const BreakIterator& other) { return b != other.b; }
  void operator ++ () {}
  };

struct BreakIterator begin(breakableIteratorType& b);
struct BreakIterator end(breakableIteratorType& b);

struct formulamode {
  rbool phi;
  ciState state;
  };

inline formulamode fmIf(rbool x) { return formulamode{x, ciIfThen}; }
inline formulamode fmIfe(rbool x) { return formulamode{x, ciIfThenElse}; }
inline formulamode fmWhile(rbool x) { return formulamode{x, ciWhile}; }

struct CondIterator begin(const formulamode& fm);
struct CondIterator end(const formulamode& fm);

struct CondIterator : ArbCondition {
  rbool phi;
  ciState state;
  const bool operator * () { return state != ciWhile && state != ciElse; }
  void go();
  void operator ++ ();
  CondIterator() {}
  CondIterator(const CondIterator& src) = delete;
  CondIterator(CondIterator&& src) {
    phi = src.phi; state = src.state; ourcontext = src.ourcontext; src.ourcontext = emptycontext;
    }  
  bool operator != (const CondIterator& other) { return state != other.state; }
  };

// we do if(!t); else to prevent adding "else" (which would not work),
// and to prevent the "unused variable" warning
#define If(x) for(bool t ## __LINE__ : fmIf(x)) if(!t ## __LINE__); else
#define Ife(x) for(bool t ## __LINE__ : fmIfe(x)) if(t ## __LINE__)
#define While(x) \
  for(bool& b ## __LINE__ : breakableIterator) \
    for(bool t ## __LINE__ : fmWhile(x)) \
      if(t##__LINE__) b##__LINE__=false; else

// lbool: lvalue formula

struct lbool {
  // inner context
  contextptr parent;
  // content
  rbool phi;

  lbool(bool b = false) : parent(currentcontext), phi(b?ftrue:ffalse) {}
  lbool(rbool ph) : parent(currentcontext), phi(ph) {}
  operator rbool () { return phi; }
  lbool& operator &= (rbool psi) {
    phi = phi && outenv(psi, true, parent); return (*this);
    }
  lbool& operator |= (rbool psi) {
    phi = phi || outenv(psi, false, parent); return (*this);
    }
  lbool& operator = (const rbool& rval) {
    if(usesenv(rval, parent)) throw assignment_exception();
    rbool alpha = outenv(ftrue, false, parent);
    phi = (alpha && rval) || ((!alpha) && phi);
    return (*this);
    }
  lbool& operator &= (bool psi) { return (*this) &= rbool(psi); }
  lbool& operator |= (bool psi) { return (*this) |= rbool(psi); }
  lbool& operator = (bool rval) { return (*this) = rbool(rval); }
  };

// -- specific relations --

// order relation:

#define ID_GREATER 1
#define ID_LEQ     2
#define ID_EQUIV   3
#define ID_MIN     4
#define ID_MAX     5

// SPASS does not accept Unicode symbols in relation names, escape them somehow
std::string spassescape(symbol t);

struct RelOrder : Relation {
  symbol opgreater, opleq, opmax, opmin;
  RelOrder(symbol _greater, symbol _leq, symbol _max, symbol _min) : 
    opgreater(_greater), opleq(_leq), opmax(_max), opmin(_min) { }
  subrelptr genSub();
  rbool less(const term& a, const term& b) { return binform(ID_GREATER, b, a); }
  rbool leq(const term& a, const term& b) { return binform(ID_LEQ, a, b); }
  std::ostream& display(std::ostream& os, int id, const term& t1, const term& t2) {
    if(outlan == SPASS) {
      if(id == ID_GREATER) 
        return os << spassescape(opgreater) << "(" << t1 << ", " << t2 << ")";
      if(id == ID_LEQ) 
        return os << "not(" << spassescape(opgreater) << "(" << t1 << ", " << t2 << "))";
      }
    if(outlan == SMT || outlan == SMT_INC) {
      if(id == ID_GREATER) 
        return os << "(> " << t1 << " " << t2 << ")";
      if(id == ID_LEQ) 
        return os << "(>= " << t2 << " " << t1 << ")";
      }
    return os << t1 << (id == ID_GREATER ? opgreater : opleq) << t2;
    }
  int negate(int id) { return 3-id; }
  rbool binform(int id, const term& a, const term& b);

  term max (const term& a, const term& b) {
    return term(std::make_shared<TermBinary> (this,ID_MAX,a,b));
    }

  term min (const term& a, const term& b) {
    return term(std::make_shared<TermBinary> (this,ID_MIN,a,b));
    }

  std::ostream& displaytermbin(std::ostream& os, const TermBinary* t) {
    int id = t->op;
    if(outlan == SPASS) {
      if(id == ID_MIN) return os << spassescape(opmin) << "(" << t->left << ", " << t->left << ")";
      if(id == ID_MAX) return os << spassescape(opmax) << "(" << t->right << ", " << t->right << ")";
      }
    return os << "(" << t->left << (id==ID_MIN?opmin:opmax) << t->right << ")";
    }

  bool worksWith(OutputLanguage lan) { return lan == L0 || lan == SPASS || lan == SMT_INC; }
  virtual std::string getName() { return "[order " + opgreater.asString() + "]"; }
  };

// ordered fields (Int and Real), which work with SMT solvers, but not with the
// internal solver

#define ID_PLUS    6
#define ID_TIMES   7
#define ID_MINUS   8
#define ID_DIVIDE  9

struct RelOrderedField : RelOrder {
  symbol opplus, optimes, opminus, opdivide;

  RelOrderedField(
    symbol _greater, symbol _leq, symbol _max, symbol _min,
    symbol _plus, symbol _times, symbol _minus, symbol _divide) : 
    RelOrder(_greater, _leq, _max, _min),
    opplus(_plus), optimes(_times),
      opminus(_minus), opdivide(_divide) { }

  term plus (const term& a, const term& b) {
    return term(std::make_shared<TermBinary> (this,ID_PLUS,a,b));
    }

  term times (const term& a, const term& b) {
    return term(std::make_shared<TermBinary> (this,ID_TIMES,a,b));
    }

  term minus (const term& a, const term& b) {
    return term(std::make_shared<TermBinary> (this,ID_MINUS,a,b));
    }

  term divide (const term& a, const term& b) {
    return term(std::make_shared<TermBinary> (this,ID_DIVIDE,a,b));
    }

  term constant(Domain *domain, int i);

  subrelptr genSub();

  std::ostream& displaytermbin(std::ostream& os, const TermBinary* t); 
  std::ostream& displayconst(std::ostream& os, const TermConst* t); 

  bool worksWith(OutputLanguage lan) { return lan == SMT_INC; }
  virtual std::string getName() { return "[real " + opplus.asString() + "]"; }

  virtual int checktermconst(int id) { throw unsolvable_exception(); }
  };

struct RelInt : RelOrderedField {

  RelInt(
    symbol _greater, symbol _leq, symbol _max, symbol _min,
    symbol _plus, symbol _times, symbol _minus, symbol _divide) : 
    RelOrderedField(_greater, _leq, _max, _min, _plus, _times, _minus, _divide) {}
  };

struct RelReal : RelOrderedField {

  RelReal(
    symbol _greater, symbol _leq, symbol _max, symbol _min,
    symbol _plus, symbol _times, symbol _minus, symbol _divide) : 
    RelOrderedField(_greater, _leq, _max, _min, _plus, _times, _minus, _divide) {}
  };

// the order used for the usual operators (<,>,<=,>=)
extern RelOrder *mainOrder;

inline rbool operator < (const term& a, const term& b) { return mainOrder->less(a, b); }
inline rbool operator > (const term& a, const term& b) { return mainOrder->less(b, a); }
inline rbool operator <= (const term& a, const term& b) { return mainOrder->leq(a, b); }
inline rbool operator >= (const term& a, const term& b) { return mainOrder->leq(b, a); }

// use this to use another order locally
// for example, if you create RelReal or RelInt, use pushorder to use it while
// it is defined

struct pushorder {
  RelOrder *orig;
  pushorder(RelOrder *r) { orig = mainOrder; mainOrder = r; }
  ~pushorder() { mainOrder = orig; }
  };

// binary relation:

#define ID_BINARY 1
#define ID_NOBINARY 2

enum loopmode { lmNoLoops, lmAllLoops, lmPossibleLoops };
enum symmode  { smSymmetric, smAsymmetric, smAntisymmetric };

struct RelBinary : Relation {
  symbol opinrel, opnotinrel;
  loopmode lm;
  symmode sm;
  RelBinary(symbol _inrel, symbol _notinrel, loopmode l, symmode s) : opinrel(_inrel), opnotinrel(_notinrel), lm(l), sm(s) { }

  subrelptr genSub();
  std::ostream& display(std::ostream& os, int id, const term& t1, const term& t2) {
    if(outlan == SPASS) {
      if(id == ID_BINARY) return os << spassescape(opinrel) << "(" << t1 << ", " << t2 << ")";
      if(id == ID_LEQ) return os << "not(" << spassescape(opinrel) << "(" << t1 << ", " << t2 << "))";
      }
    return os << t1 << (id == ID_BINARY ? opinrel : opnotinrel) << t2;
    }

  int negate(int id) { return 3-id; }

  rbool binform(int id, const term& a, const term& b);
  
  rbool operator () (const term& a, const term& b) { return binform(ID_BINARY, a, b); }

  virtual bool worksWith(OutputLanguage lan) { return lan == L0; }

  virtual std::string getName() { return "[binary " + opinrel.asString() + "]"; }
  };

#define ID_UN_EQ -1
#define ID_UN_NEQ -2

// unary relation:

struct RelUnary : Relation {
  symbol oprel;
  int nogroups;
  RelUnary(symbol rel, int n) : oprel(rel), nogroups(n) { }

  subrelptr genSub();
  std::ostream& display(std::ostream& os, int id, const term& t1, const term& t2);
  int negate(int id);
  rbool binform(int id, const term& a, const term& b);
  
  // is 'a' in the group number 'v' (0-based)?
  rbool operator () (const term& a, int v) { return binform(1<<v, a, a); }

  // are 'a' and 'b' in the same group?  
  rbool together(const term& a, const term& b) { return binform(ID_UN_EQ, a, b); }

  virtual bool worksWith(OutputLanguage lan) { return lan == L0; }
  virtual std::string getName() { return "[unary " + oprel.asString() + "]"; }
  };

// infinite homogeneous tree:

#define ID_ANCESTOR_OR_EQUAL 1
#define ID_NOT_ANCESTOR_OR_EQUAL 2

struct RelTree : Relation {
  symbol opanceq, opnotanceq, oplca;
  RelTree(symbol rel, symbol notrel, symbol join)
    : opanceq(rel), opnotanceq(notrel), oplca(join) { }

  subrelptr genSub();

  std::ostream& display(std::ostream& os, int id, const term& t1, const term& t2) {
    if(outlan == SPASS) {
      if(id == ID_ANCESTOR_OR_EQUAL) return os << spassescape(opanceq) << "(" << t1 << ", " << t2 << ")";
      if(id == ID_NOT_ANCESTOR_OR_EQUAL) return os << "not(" << spassescape(opanceq) << "(" << t1 << ", " << t2 << "))";
      }
    return os << t1 << (id == ID_ANCESTOR_OR_EQUAL ? opanceq : opanceq) << t2;
    }

  std::ostream& displaytermbin(std::ostream& os, const TermBinary *t) {
    return os << "(" << t->left << oplca << t->right << ")";
    }

  int negate(int id) { 
    return 3-id;
    }

  virtual bool worksWith(OutputLanguage lan) { return lan == L0; }

  rbool binform(int id, const term& a, const term& b);
  
  rbool anceq (const term& a, const term& b) { 
    return binform(ID_ANCESTOR_OR_EQUAL, a, b); 
    }

  term lca (const term& a, const term& b) {
    return term(std::make_shared<TermBinary> (this,0,a,b));
    }

  virtual std::string getName() { return "[tree " + opanceq.asString() + "]"; }
  };

// -- sets --


// --- ElementOf ---

// represent the type T as a subclass of Element

// a simple set (set builder expression)

template<class T> struct SimpleSetOf {
  varlist var; // list of variables in the context (rhs of the set builder expression)
  rbool phi; // the constraint (one is enough)
  T a;
  SimpleSetOf() {}
  SimpleSetOf(const SimpleSetOf& src) {
    phi = src.phi; a = src.a; var = src.var;
    }
  SimpleSetOf(varlist _var, rbool _phi, T _a) : var(_var), phi(_phi), a(_a) {}
  SimpleSetOf(SimpleSetOf&& src) {
    phi = src.phi; a = src.a; var = src.var;
    src.var.resize(0);
    }
  };

// iterator for ESet (very technical)

template<class T> struct EIteratorOf {
  const lsetof<T> &s;
  int index;
  std::shared_ptr<T> at;
  contextptr ourcontext;
  T& operator * ();
  void operator ++ ();
  void connectIterator();
  void disconnectIterator();
  EIteratorOf(const lsetof<T> &_s, int id) : s(_s), index(id) { 
    connectIterator();
    }
  EIteratorOf(const EIteratorOf<T>& src) = delete;
  EIteratorOf(EIteratorOf<T>&& src) : s(src.s) {
    index = src.index; at = src.at; ourcontext = src.ourcontext;
    }
  ~EIteratorOf() { 
    disconnectIterator(); 
    }
  bool operator != (const EIteratorOf<T>& other);
  };

// create a set with the given contents

template<class T> lsetof<T> newSet(T x) {
  lsetof<T> res; res.insert(x, currentcontext); return res;   
  }

template<class T> lsetof<T> newSet(T x, T y) {
  lsetof<T> res; 
  res.insert(x, currentcontext); 
  res.insert(y, currentcontext); 
  return res;   
  }

template<class T> lsetof<T> newSet(std::initializer_list<T> l) {
  lsetof<T> res;
  for(auto z: l) res.insert(z, currentcontext);
  return res;
  }

inline bool isused(vptr v, const term& a) { return a.p->uses(v); }

// negation structure: 
// a technical construct required to enable writing the set difference as (A&&~B)
//-------------------------------------------------------------------------------

template<class T> struct negated {
  T original;
  negated(T t) { original = t; }
  };

template<class T> T operator ~ (negated<T> A) { return A.original; }

// set/setof lvalues
//-------------------

template<class T> struct lsetof {
  typedef T element;
  contextptr ain; // environment at the time of definition

  std::vector<std::shared_ptr<SimpleSetOf<T>>> elts;
  std::ostream& display (std::ostream &os) const;
  void insert(T a, contextptr ain, contextptr nowcontext = currentcontext);
  bool uses(vptr w);
  rbool subseteq (const lsetof<T>& e) const;
  rbool hasmember (const T& e) const;
  rbool operator == (const lsetof<T>& e) const {
    return subseteq(e) && e.subseteq(*this);
    }
  rbool operator != (const lsetof<T>& e) const { return !((*this) == e); }
  rbool operator <= (const lsetof<T>& e) const {
    return subseteq(e);
    }
  rbool operator < (const lsetof<T>& e) const {
    return subseteq(e) && !e.subseteq(*this);
    }
  
  lsetof() { ain = currentcontext; }

  lsetof<T> removeall();
  lsetof<T> removeallnonset();

  lsetof<T>& operator += (const T& y) { insert(y, ain); return *this; }
  lsetof<T>& operator |= (const lsetof<T>& y) { for(const T& e: y) insert(e, ain); return *this; }
  lsetof<T>& operator -= (const T& y);
  lsetof<T>& operator &= (const lsetof<T>& y);
  lsetof<T>& operator &= (const negated<lsetof<T>>& y);

  EIteratorOf<T> begin() const { return EIteratorOf<T> (*this, 0); }
  EIteratorOf<T> end() const { return EIteratorOf<T> (*this, 2000000000); }

  lsetof<T>& operator = (const lsetof<T>& rval) {
    if(debuglois) std::cout << ido << "assign " << rval << " ain = " << ain << " context = " << emptycontext << std::endl;
    autoindenter i1;
    if(usesenv(rval, ain)) throw assignment_exception();
    removeall(); (*this) |= rval;
    return *this;
    }

  lsetof<T>(const lsetof<T> &x) : ain(currentcontext) {
    if(debuglois) std::cout << ido << "copying " << x << " ain = " << emptycontext << std::endl;
    autoindenter i2;
    (*this) |= x;
    }

  lsetof<T>(lsetof<T> &&x) : ain(currentcontext) {
    if(x.ain == currentcontext) {
      swap(elts, x.elts);
      }
    else {
      (*this) |= x;
      }
    }
  };

template<class T> rbool isempty(lsetof<T>& what) {
  lbool res = ftrue;
  for(auto e: what) res = ffalse;
  return res;
  }

template<class T> lsetof<T> substitute(const lsetof<T>& where, const varsubstlist& l) {
  lsetof<T> res;
  for(auto elt: where.elts)
    res.elts.emplace_back(std::make_shared<SimpleSetOf<T>> (
      elt->var,
      substitute(elt->phi, l),
      substitute(elt->a, l))
      );
  return res;
  }

template<class T> bool isused(vptr w, const lsetof<T>& where) { 
  for(auto e: where.elts)
    if(isused(w, e->a) || isused(w, e->phi)) 
      return true;
  return false;
  }

template<class T>
  negated<lsetof<T>> operator ~ (const lsetof<T>& A) { return negated<lsetof<T>> (A); }

template<class T> std::ostream& operator << (std::ostream & os, const SimpleSetOf<T>& e) {
  os << e.a; // os << e.a;
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

template<class T> std::ostream& lsetof<T>::display (std::ostream &os) const { 
  if(elts.size() == 0) return os << sym.emptyset;
  os << sym.leftbrace;
  for(size_t i=0; i<elts.size(); i++) {
    if(i) os << sym.ssunion;
    os << *elts[i];
    }
  return os << sym.rightbrace;
  }

template<class T> 
std::ostream& operator << (std::ostream& os, const lsetof<T>& A) { return A.display(os); }

template<class T> T& EIteratorOf<T>::operator * () {
  return *at;
  }
template<class T> void EIteratorOf<T>::operator ++ () { 
  disconnectIterator();
  // if(index >= 0) 
  index++;
  connectIterator();
  }
template<class T> void EIteratorOf<T>::connectIterator() {
  if(index >= s.elts.size() || index < 0) return;
  const SimpleSetOf<T>& sel(*s.elts[index]);
  // add variables and alph!

  if(debuglois) {
    std::cout << ido << "connect " << sel.a << " phi=" << sel.phi;
    for(auto v: sel.var) std::cout << " v:"<<v;
    std::cout << " in context="<<emptycontext << std::endl;
    }
  autoindenter i1;

  varsubstlist varpairs;
  varlist vlist;

  for(auto w: sel.var) {
    vptr v = newvar(w->dom);
    varpairs.emplace_back(w, term(v));
    vlist.push_back(v);
    }

  currentcontext = ourcontext = std::make_shared<Context> (currentcontext, 
    substitute(sel.phi, varpairs), vlist);
  
  // sometimes the environment is inconsistent with phi
  if(debuglois) std::cout << ido << "got for conscheck: " << *at << std::endl << ido << "  IN " << emptycontext << std::endl;
  int res = solver->solveEnv();
  
  if(debuglois) std::cout << ido << "checking consistency: " << res << std::endl;
  
  if(res == 2) at = std::make_shared<T> (substitute(sel.a, varpairs));
  else ++(*this);

  if(debuglois) {
    if(at == NULL) std::cout << ido << "no result\n" << std::endl;
    else std::cout << ido << "result: " << *at << " in context="<<emptycontext << std::endl;
    }
  }
template<class T> void EIteratorOf<T>::disconnectIterator() {
  // printf("%p: disconnect\n", this);
  if(ourcontext) {
    if(currentcontext != ourcontext) throw iterator_exception();
    currentcontext = ourcontext->parent;
    ourcontext = NULL;
    at = NULL;
    }
  }

template<class T> bool EIteratorOf<T>::operator != (const EIteratorOf<T>& other) {
  if(index < 0 || index >= s.elts.size())
    return !(other.index < 0 || other.index >= other.s.elts.size());
  return index != other.index;
  }

template<class T> rbool lsetof<T>::subseteq (const lsetof<T>& e) const {
  lbool a(ftrue); for(auto& x: *this) a &= e.hasmember(x); return a;
  }

template<class T> rbool lsetof<T>::hasmember(const T& e) const {
  lbool a(ffalse); for(auto& x: *this) a |= (e == x); return a;
  }

void aggressive_simplify(varlist& selvar, rbool& selphi, contextptr& ain, std::function<bool(vptr)> isused, std::function<void(const varsubstlist&)> alpha);

template<class T> void aggressive_simplify(varlist& selvar, rbool& selphi, T& a, contextptr& ain) {
  aggressive_simplify(selvar, selphi, ain,
    [&] (vptr v) { return isused(v,a); },
    [&] (const varsubstlist& l) { a = substitute(a,l); }
    );    
  }

term valueKnown(rbool& phi, vptr v, bool negated);
rbool makequantifier(bool f, rbool& phi, varlist to_quantify);

template<class T> void lsetof<T>::insert(T a, contextptr ain, contextptr nowenv) {

  if(debuglois) std::cout << ido << "insert " << a << ieol << "  ENV " << emptycontext << ieol << "  AIN "<< ain << std::endl;
  autoindenter i1;
  
  varlist selvar;
  rbool selphi = ftrue;
  
  varlist to_quantify;
  
  contextptr ainlocal = nowenv;  

  while(ainlocal != ain) {
    if(!ainlocal) throw env_exception();
    rbool aphi = ainlocal->phi;
    selphi = selphi && aphi;
    for(auto w: ainlocal->var) 
      if(isused(w, selphi) || isused(w, a)) {
        term w1 = valueKnown(selphi, w, false);
        if(debuglois) std::cout << ido << "checking variable: " << w << std::endl;
        
        if(w1.p) {
          if(debuglois) std::cout << ido << "value known to be " << w1 << std::endl;
          selphi = substitute(selphi, w, w1), a = substitute(a, w, w1);
          }
        else if(isused(w, a)) {
          if(debuglois) std::cout << ido << "value of " << w << " not known, added to selvar" << std::endl;
          selvar.push_back(w);
          }
        else {
          if(debuglois) std::cout << ido << "value of " << w << " not known, added to to_quantify" << std::endl;
          to_quantify.push_back(w);
          }
        }
    ainlocal = ainlocal->parent;
    }
  if(to_quantify.size()) {
    if(debuglois) std::cout << ido << "phi (1) = " << selphi << std::endl;
    if(debuglois) for(vptr v: to_quantify) std::cout << ido << "quantifier to add: " << v << std::endl;
    selphi = makequantifier(false, selphi, to_quantify);
    if(debuglois) std::cout << ido << "phi (2) = " << selphi << std::endl;
    }
#ifdef AGGSYM
  aggressive_simplify(selvar, selphi, a, ain);
#endif

  if(debuglois) {
    for(vptr v: selvar) 
    std::cout << ido << "| selvar = " << v << std::endl;
    std::cout << ido << "| selphi = " << selphi << std::endl;
    std::cout << ido << "| a      = " << a << std::endl;
    }

  elts.emplace_back(std::make_shared<SimpleSetOf<T>>(selvar, selphi, a));
  }


// initialize everything

void initLois();

}

#endif

