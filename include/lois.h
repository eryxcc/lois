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

struct indenter { int i; bool doit; };

extern indenter iindent, ieol, iunindent;

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
struct elem;
struct lset;
struct Subdomain;
struct rbool;
struct Element;
struct ESet;

template<class T> struct lsetof;

typedef std::shared_ptr<struct Variable> vptr; 
typedef std::shared_ptr<struct SubRelation> subrelptr;
typedef std::shared_ptr<struct TermVariable> tvptr;
typedef std::vector<vptr> varlist;
typedef std::shared_ptr<struct Context> contextptr;

// convert an elem to a term
term& asa(const elem& x);

// convert an elem to type T, throw as_exception if it does not contain type T
template<class T> T& as(elem a);

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
  void display(std::ostream& os);
  vptr asVar() const;
  };

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
  virtual term alphterm(const term& ths, vptr v1, vptr v2) const = 0;
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
  virtual term alphterm(const term& ths, vptr v1, vptr v2) const { 
    if(this == &(*v1)) return term(v2); else return ths;
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
  virtual term alphterm(const term& ths, vptr v1, vptr v2) const {
    return term(std::make_shared<TermBinary>(r, op, left.p->alphterm(left,v1,v2), right.p->alphterm(right,v1,v2)));
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
  virtual term alphterm(const term& ths, vptr v1, vptr v2) const {
    return term(std::make_shared<TermConst>(r, op, v->dom));
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

// alpha: substitute the variable (v1 is replaced by v2)
vptr alpha(vptr vthis, vptr v1, vptr v2);
term alpha(term tthis, vptr v1, vptr v2);

// isused: is the variable v used in the second parameter?
bool isused(vptr vthis, vptr v);
bool isused(vptr v, rbool phi);

inline std::ostream& operator << (std::ostream& os, const term& a) { 
  a.p->display(os); return os;
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
  Context(contextptr par, rbool fi) : parent(par), phi(fi) {}
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

// the set of extensions of anccontext-valuation to nowcontext
lset branchset(contextptr anccontext = emptycontext, contextptr nowcontext = currentcontext);

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
  rbool operator () (const elem& a, const elem& b) {
    return binform(ID_BINARY, asa(a), asa(b));
    }

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
  rbool operator () (const elem& a, int v) { return (*this) (asa(a), v); }

  // are 'a' and 'b' in the same group?  
  rbool together(const term& a, const term& b) { return binform(ID_UN_EQ, a, b); }
  rbool together(elem& a, elem& b) { return binform(ID_UN_EQ, asa(a), asa(b)); }

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

  rbool anceq (const elem& a, const elem& b) { 
    return binform(ID_ANCESTOR_OR_EQUAL, asa(a), asa(b)); 
    }
  
  term lca (const term& a, const term& b) {
    return term(std::make_shared<TermBinary> (this,0,a,b));
    }

  term lca (const elem& a, const elem& b) {
    return term(std::make_shared<TermBinary> (this,0,asa(a),asa(b)));
    }

  virtual std::string getName() { return "[tree " + opanceq.asString() + "]"; }
  };

// -- sets --

typedef std::pair<elem, elem> elpair;
typedef std::vector<elem> eltuple;

struct elem {
  std::shared_ptr<struct Element> p; 
  elem(std::shared_ptr<struct Element> e) : p(e) {}
  elem() {}
  Element* operator -> () const { return &(*p); }
  elem(const lset& s);
  elem(std::shared_ptr<struct ESet> s);
  elem(int);
  elem(vptr);
  elem(term);
  elem(elpair);
  elem(eltuple);
  };

extern const elem nullelem;

// elem is internally a shared_ptr to Element

struct Element {
  // destructor
  virtual ~Element() {}
  // display on a stream
  virtual std::ostream& display (std::ostream &os) const = 0;
  // alpha-convert this from 'v1' to 'v2' (i.e., substitute variable)
  virtual elem alph(vptr v1, vptr v2) const = 0;
  // is the free variable 'v' used by this element?
  virtual bool uses(vptr v) = 0;
  virtual rbool operator == (Element& e) = 0;
  virtual rbool operator != (Element& e) = 0;
  virtual rbool operator <  (Element& e) = 0;
  virtual rbool operator <= (Element& e) = 0;
  };

inline rbool operator == (elem x, elem y) { return (*x.p) == (*y.p); }
inline rbool operator != (elem x, elem y) { return (*(x.p)) != (*(y.p)); }
inline rbool operator <  (elem x, elem y) { return (*x.p) < (*y.p); }
inline rbool operator <= (elem x, elem y) { return (*(x.p)) <= (*(y.p)); }
inline rbool operator >  (elem x, elem y) { return (*y.p) < (*x.p); }
inline rbool operator >= (elem x, elem y) { return (*(y.p)) <= (*(x.p)); }

inline bool isused(vptr v, elem e) { return e->uses(v); }

inline std::ostream& operator << (std::ostream& os, elem a) { return a->display(os); }

// --- ElementOf ---

// represent the type T as a subclass of Element

template<class T> struct ElementOf: Element {
  T data;
  explicit ElementOf(T&& tmp, int) : data(tmp) {}
  explicit ElementOf(const T& tmp, int) : data(tmp) {}
  ElementOf(T tmp) : data(tmp) {}
  std::ostream& display (std::ostream &os) const { return os << data; }
  virtual elem alph(vptr v1, vptr v2) const { 
    return elem(std::make_shared<ElementOf> (alpha(data, v1, v2)));
    }
  virtual bool uses(vptr v) { return isused(v, data); }
  rbool operator == (Element& e) {
    auto e2 = dynamic_cast<ElementOf<T>*>(&e);
    if(!e2) return ffalse;
    return data == e2->data;
    }
  rbool operator != (Element& e) {
    auto e2 = dynamic_cast<ElementOf<T>*>(&e);
    if(!e2) return ftrue;
    return data != e2->data;
    }
  rbool operator < (Element& e) {
    auto e2 = dynamic_cast<ElementOf<T>*>(&e);
    if(!e2) return ffalse;
    return data < e2->data;
    }
  rbool operator <= (Element& e) {
    auto e2 = dynamic_cast<ElementOf<T>*>(&e);
    if(!e2) return ffalse;
    return data <= e2->data;
    }
  };

// encapsulate T as ElementOf<T>

template<class T> elem elof(const T& x) { return elem(std::make_shared<ElementOf<T>> (x)); }
template<class T> elem elof(T&& x) { return elem(std::make_shared<ElementOf<T>> (x)); }

// ElementOf<T> back to T
template<class T> T& as(elem a) {
  auto sa = std::dynamic_pointer_cast<ElementOf<T>> (a.p);
  if(!sa) throw as_exception();
  return sa->data;
  }

// is this element of type T?
template<class T> bool is(elem a) {
  auto sa = std::dynamic_pointer_cast<ElementOf<T>> (a.p);
  return sa ? true : false;
  }

// a simple set (set builder expression)

struct SimpleSet {
  varlist var; // list of variables in the context (rhs of the set builder expression)
  rbool phi; // the constraint (one is enough)
  elem a;
  SimpleSet() {}
  SimpleSet(const SimpleSet& copy_from_me) = delete;
  SimpleSet(varlist _var, rbool _phi, elem _a) : var(_var), phi(_phi), a(_a) {}
  SimpleSet(SimpleSet&& src) {
    phi = src.phi; a = src.a; var = src.var;
    src.var.resize(0);
    }
  };

// a set

struct ESet : Element {
  ESet() { }
  std::vector<struct SimpleSet> elts;
  virtual struct EIterator begin() const;
  virtual struct EIterator end() const;
  std::ostream& display (std::ostream &os) const;
  virtual elem alph(vptr v1, vptr v2) const;
  virtual void insert(elem a, contextptr ain, contextptr nowcontext = currentcontext);
  virtual bool uses(vptr w);
  rbool subseteq (ESet* e);
  rbool hasmember (elem e);
  rbool operator == (Element& e);
  rbool operator != (Element& e) { return !((*this) == e); }
  rbool operator <= (Element& e) {
    auto e2 = dynamic_cast<ESet*>(&e);
    if(!e2) return ffalse;
    return subseteq(e2);
    }
  rbool operator < (Element& e) {
    auto e2 = dynamic_cast<ESet*>(&e);
    if(!e2) return ffalse;
    return subseteq(e2) && !e2->subseteq(this);
    }
  };

inline elem::elem(std::shared_ptr<struct ESet> s) : p(s) {}

// iterator for ESet (very technical)

struct EIterator {
  const ESet *s;
  int index;
  elem at;
  contextptr ourcontext;
  const elem& operator * ();
  void operator ++ ();
  void connectIterator();
  void disconnectIterator();
  EIterator(const ESet *_s, int id) : s(_s), index(id) { 
    at = nullelem; if(s) connectIterator();
    }
  EIterator(const EIterator& src) = delete;
  EIterator(EIterator&& src) {
    s = src.s; index = src.index; at = src.at; ourcontext = src.ourcontext;
    src.at = nullelem;
    }
  ~EIterator() { 
    disconnectIterator(); 
    }
  bool operator != (const EIterator& other);
  };

EIterator begin(elem a);
EIterator end(elem a);

// create a set with the given contents
inline std::shared_ptr<ESet> newInternalSet() { return std::make_shared<ESet> (); }

lset newSet(elem x);
lset newSet(elem x, elem y);
lset newSet(std::initializer_list<elem> l);

inline bool isused(vptr v, const term& a) { return a.p->uses(v); }

inline elem alpha(elem x, vptr v1, vptr v2) { return x.p->alph(v1, v2); }

extern lset emptyset;

bool isSet(elem a);

lset asSet(elem a);

// setof

template<class T> struct EIteratorOf {
  EIterator it;
  EIteratorOf(EIterator&& i) : it(NULL, 0) {
    it.s = i.s; it.index = i.index; it.at = i.at; 
    it.ourcontext = i.ourcontext; i.at = nullelem;
    }
  EIteratorOf(const EIterator& i) = delete;
  EIteratorOf& operator ++ () { ++it; return *this; }
  T& operator* () { return as<T> (*it); }
  EIteratorOf(EIteratorOf&& src) : it(NULL, 0) { 
    it.s = src.it.s; it.index = src.it.index; it.at = src.it.at; 
    it.ourcontext = src.it.ourcontext; src.it.at = nullelem;
    }
  bool operator != (const EIteratorOf<T> & x) { return it != x.it; }
  };

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

struct lset {
  std::shared_ptr<struct ESet> p;
  contextptr ain; // environment at the time of definition
  
  ESet* operator -> () const { return &(*p); }
  
  lset() { p = newInternalSet(); ain = currentcontext; }
  lset(std::shared_ptr<struct ESet> e) : ain(currentcontext), p(e) {}
  explicit lset(elem e) : ain(currentcontext), p(std::dynamic_pointer_cast<ESet>(e.p)) { if(!p) throw as_exception(); }

  lset removeall();
  lset removeallnonset();

  lset& operator += (elem y) { p->insert(y, ain); return *this; }
  lset& operator |= (const lset& y) { for(elem e: y) p->insert(e, ain); return *this; }
  lset& operator -= (elem y);
  lset& operator &= (const lset& y);
  lset& operator &= (negated<lset> y);

  struct EIterator begin() const { return p->begin(); }
  virtual struct EIterator end() const { return p->end(); }

  template<class T>
  lset(const lsetof<T>& o) { ain = currentcontext; p = newInternalSet(); (*this) |= o.orig; }

  lset& operator = (const lset& rval) {
    if(usesenv(rval, ain)) throw assignment_exception();
    removeall(); (*this) |= rval;
    return *this;
    }

  lset(const lset &x) : ain(currentcontext) {
    p = newInternalSet();
    (*this) |= x;
    }

  lset(lset &&x) : ain(currentcontext) {
    if(x.ain == currentcontext) {
      printf("quick move used\n");
      p = x.p;
      }
    else {
      p = newInternalSet();
      (*this) |= x;
      }
    }
  };

inline elem::elem(const lset& s) : p(s.p) {}
inline rbool operator == (const lset& x, const lset& y) { return (*x.p) == (*y.p); }
inline rbool operator != (const lset& x, const lset& y) { return (*x.p) != (*y.p); }
inline bool isused(vptr v, const lset& e) { return e->uses(v); }
inline lset newSet() { return lset(); }

inline negated<lset> operator ~ (const lset& A) { return negated<lset>(A); }
template<class T>
  negated<lsetof<T>> operator ~ (const lsetof<T>& A) { return negated<lsetof<T>> (A); }

inline std::ostream& operator << (std::ostream& os, const lset& a) { return a->display(os); }

template<class T> struct lsetof {
  lset orig;
  
  lsetof() { }
  lsetof(const lsetof<T>& r) : orig(r.orig) {}
  lsetof(const lset& r) : orig(r) {}

  lsetof<T> removeall() { return lsetof<T> (orig.removeall()); }

  lsetof<T>& operator += (T y) { orig += elof<T>(y); return *this; }
  lsetof<T>& operator |= (const lsetof<T>& x) { orig |= x.orig; return (*this); }

  EIteratorOf<T> begin() const { return EIteratorOf<T> (orig.begin()); }
  EIteratorOf<T> end() const { return EIteratorOf<T> (orig.end()); }
  rbool isEmpty() { return orig == emptyset; }
  
  operator lset () const { return orig; }
  };

template<class T> std::ostream& operator << (std::ostream& os, const lsetof<T>& A) {
  A.orig.p->display(os); return os;
  }

template<class T> lsetof<T> alpha(lsetof<T> A, vptr v1, vptr v2) {
  return lsetof<T> (A.orig.alph(v1, v2));
  }

template<class T> bool isused(vptr v, lsetof<T> A) { return A.orig->uses(v); }

template<class T> lsetof<T> newSetOf() { return lsetof<T> (newSet()); }
template<class T> lsetof<T> newSetOf(T x) { return lsetof<T> (newSet(elof<T> (x))); }
template<class T> lsetof<T> newSetOf(T x, T y) { return lsetof<T> (newSet(elof<T> (x), elof<T> (y))); }

// initialize everything

void initLois();

}

#endif

