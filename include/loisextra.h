// this file includes features which are non the strict core of LOIS
// anyway, most LOIS programs would use them

#ifndef _loisextra_h_
#define _loisextra_h_

#include "lois.h"

namespace lois {

// the basic elements
//--------------------

// integer constant

// inline elem newConst(int name) { return newEl<int> (name); }
inline bool isused(vptr v, int name) { return false; }
inline int alpha(int name, vptr v1, vptr v2) { return name; }

// element pair

inline bool isused(vptr v, elpair p) { return isused(v, p.first) || isused(v, p.second); }
inline elpair alpha(elpair p, vptr v1, vptr v2) { 
  return std::make_pair(alpha(p.first, v1, v2), alpha(p.second, v1, v2));
  }

inline rbool operator == (elpair a, elpair b) { return a.first == b.first && a.second == b.second; }
inline rbool operator != (elpair a, elpair b) { return a.first != b.first || a.second != b.second; }
inline rbool operator <  (elpair a, elpair b) { return a.first < b.first || (a.first == b.first && a.second < b.second); }
inline rbool operator <= (elpair a, elpair b) { return a.first < b.first || (a.first == b.first && a.second <= b.second); }

inline std::ostream& operator << (std::ostream& os, elpair p) { return os << "(" << p.first << "," << p.second << ")"; }

// element tuple

inline bool isused(vptr v, eltuple t) { 
  for(elem e: t) if(e->uses(v)) return true;
  return false;
  }

inline eltuple alpha(eltuple t, vptr v1, vptr v2) {
  eltuple v; for(elem e:t) v.push_back(alpha(e, v1, v2)); return v;
  }

rbool operator == (eltuple a, eltuple b);
inline rbool operator != (eltuple a, eltuple b) { return !(a==b); }

rbool operator <  (eltuple a, eltuple b);
rbool operator <= (eltuple a, eltuple b);

std::ostream& operator << (std::ostream& os, eltuple t);

// negation structure: 
// a technical construct required to enable writing the set difference as (A&&~B)
//-------------------------------------------------------------------------------

template<class T> struct negated {
  T original;
  negated(T t) { original = t; }
  };

inline negated<rset> operator ~ (rset A) { return negated<rset>(A); }
template<class T>
  negated<rsetof<T>> operator ~ (rsetof<T> A) { return negated<rsetof<T>> (A); }
template<class T> T operator ~ (negated<T> A) { return A.original; }

// set/setof lvalues
//-------------------

struct lset {
  rset orig;
  contextptr ain; // environment at the time of definition
  
  lset() { orig = newSet(); ain = currentcontext; }

  rset removeall();
  rset removeallnonset();

  lset& operator += (elem y) { orig->insert(y, ain); return *this; }
  lset& operator |= (rset y) { for(elem e: y) orig->insert(e, ain); return *this; }
  lset& operator -= (elem y);
  lset& operator &= (rset y);
  lset& operator &= (negated<rset> y);

/*lset& operator += (int i) { return (*this) += elof<int> (i); }
  lset& operator += (elpair x) { return (*this) += elof<elpair> (x); }
  lset& operator += (eltuple x) { return (*this) += elof<eltuple> (x); }
  lset& operator += (rset x) { return (*this) += elem(x); }
*/
  operator rset() const { lset copy; copy |= orig; return copy.orig; }
  operator elem() { lset copy; copy |= orig; return copy.orig; }

  struct EIterator begin() const { return orig.p->begin(); }
  virtual struct EIterator end() const { return orig.p->end(); }

  lset(rset o) { ain = currentcontext; orig = newSet(); (*this) |= o; }
  template<class T>
  lset(rsetof<T> o) { ain = currentcontext; orig = newSet(); (*this) |= o.orig; }

  lset& operator = (const rset& rval) {
    if(usesenv(rval, ain)) throw assignment_exception();
    removeall(); (*this) |= rval;
    return *this;
    }

  lset& operator = (const lset& rval) { return (*this) = rval.orig; }
  lset(const lset &x) { orig = newSet(); ain = currentcontext; (*this) |= x; }
  };

inline rbool operator == (const lset& a, const lset& b) { return a.orig == b.orig; }
inline rbool operator != (const lset& a, const lset& b) { return a.orig != b.orig; }
inline rbool operator == (const lset& a, const rset& b) { return a.orig == b; }
inline rbool operator != (const lset& a, const rset& b) { return a.orig != b; }
inline rbool operator == (const rset& a, const lset& b) { return a == b.orig; }
inline rbool operator != (const rset& a, const lset& b) { return a != b.orig; }

inline std::ostream& operator << (std::ostream& os, const lset& s) { return os <<
    s.orig; }

template<class T> struct lsetof {
  lset orig;
  
  lsetof() { }
  lsetof(rsetof<T> r) : orig(r.orig) {}

  rsetof<T> removeall() { return rsetof<T> (orig.removeall()); }

  lsetof<T>& operator += (T y) { orig += elof<T>(y); return *this; }
  lsetof<T>& operator |= (rsetof<T> x) { orig |= x.orig; return (*this); }

  EIteratorOf<T> begin() { return EIteratorOf<T> (orig.begin()); }
  EIteratorOf<T> end() { return EIteratorOf<T> (orig.end()); }
  rbool isEmpty() { return orig == emptyset; }
  
  operator rsetof<T> () const { return rsetof<T> (rset(orig)); }
  rsetof<T> r() const { return rsetof<T> (rset(orig)); }
  operator lset () { return orig; }
  operator rset () const { return rset(orig); }
  };

template<class T> std::ostream& operator << (std::ostream& os, lsetof<T> A) {
  A.orig.orig->display(os); return os;
  }

template<class T> rbool operator == (const lsetof<T>& a, const lsetof<T>& b) { return a.orig == b.orig; }
template<class T> rbool operator != (const lsetof<T>& a, const lsetof<T>& b) { return a.orig != b.orig; }
template<class T> rbool operator == (const lsetof<T>& a, const rsetof<T>& b) { return a.orig == b; }
template<class T> rbool operator != (const lsetof<T>& a, const rsetof<T>& b) { return a.orig != b; }
template<class T> rbool operator == (const rsetof<T>& a, const lsetof<T>& b) { return a == b.orig; }
template<class T> rbool operator != (const rsetof<T>& a, const lsetof<T>& b) { return a != b.orig; }

// basic set theoretic operations
//--------------------------------

inline rbool memberof(elem x, rset y) {
  lbool af(ffalse); for(elem a: y) af |= (x == a); return af;
  }

inline rbool memberof(elem x, const lset& y) {
  lbool af(ffalse); for(elem a: y) af |= (x == a); return af;
  }

template<class T> inline rbool memberof(T x, rsetof<T> y) {
  lbool af(ffalse); for(T a: y) af |= (x == a); return af;
  }

template<class T> inline rbool memberof(T x, const lsetof<T>& y) {
  lbool af(ffalse); for(T a: y.r()) af |= (x == a); return af;
  }

inline rbool subseteq(rset x, rset y) {
  lbool af(ftrue); for(elem a: x) af &= memberof(a, y); return af;
  }

inline rset operator & (rset x, rset y) { 
  lset xy;
  for(elem e: x) If(memberof(e, y)) xy += e;
  return xy;
  }

inline rset operator & (rset x, negated<rset> y) { 
  lset xy;
  for(elem e: x) If(!memberof(e, y.original)) xy += e;
  return xy;
  }

inline rset operator | (rset x, rset y) { 
  lset xy;
  for(elem e: x) xy += e;
  for(elem e: y) xy += e;
  return xy;
  }

template<class T> inline rsetof<T> operator | (rsetof<T> x, rsetof<T> y) { 
  lsetof<T> xy;
  for(elem e: x) xy += e;
  for(elem e: y) xy += e;
  return xy;
  }

// remove the repetitions: each element will be listed only once in the optimized set
rset optimize(rset x);

// optimize the elements of given type
rset optimizeType(rset x, rset type);

// the Cartesian product
rset operator * (rset x, rset y);

// the Cartesian product of more than 2
rset cartesian(std::initializer_list<elem> l, eltuple v = eltuple());

// macros
//--------

#define FORALL(x, set, condition) \
  ([=]() -> rbool { lbool ret(ftrue); for(auto x: (set)) ret &= (condition); return ret; }())
#define EXISTS(x, set, condition) \
  ([=]() -> rbool { lbool ret(ffalse); for(auto x: (set)) ret |= (condition); return ret; }())
#define FILTERMAP(x, set, condition, replace) \
  ([=]() -> rset { lset S; for(auto x: (set)) If(condition) S += (replace); return S; }())
#define FILTER(x, set, condition) FILTERMAP(x, set, condition, x)
#define MAP(x, set, replace) FILTERMAP(x, set, ftrue, replace)

// axioms & declareatom
//----------------------

// declare an axiom
struct axiom : ArbCondition {
  axiom(rbool phi);
  };

// declare an atom
struct declareatom : ArbCondition {
  elem at;
  declareatom(Domain *d, const std::string& s);
  operator elem() { return at; }
  term asTerm() { return as<term> (at); }
  };

// singleton construction: [r/l]elem
//-----------------------------------

extern struct extracttype {} EXTRACT;

struct relem {
  rset singleton;
  relem() { singleton = newSet(); }
  relem(elem b) { singleton = newSet(b); }
  relem(extracttype, rset X) : singleton(X) {}
  rbool isUndefined() { return singleton == emptyset; }
  };

// extract the element from a singleton set, and vice versa
inline relem extract(rset e) { return relem(EXTRACT, e); }
inline rset newSet(relem el) { return el.singleton; }

inline rbool operator == (relem a, relem b) 
  { return a.singleton == b.singleton; }

inline rbool operator != (relem a, relem b) 
  { return a.singleton != b.singleton; }

inline rbool operator > (relem a, relem b) {
  lbool p = true;
  for(elem x: a.singleton) for(elem y: b.singleton) p &= x>y;
  return p;
  }

inline std::ostream& operator << (std::ostream& os, relem e) {
  return os << sym.pseudo << e.singleton;
  }

struct lelem {
  lset singleton;
  lelem() {}
  lelem(elem b) { singleton |= newSet(b); }
  lelem(extracttype, rset X) : singleton(X) {}
  rbool isUndefined() { return singleton == emptyset; }  

  lelem(relem b) { singleton |= b.singleton; }
  operator relem() const { return extract(rset(singleton)); }
  
  lelem& operator = (const elem& rval) {
    if(usesenv(rval, singleton.ain)) throw assignment_exception();
    singleton.removeall();
    singleton += rval;
    return *this;
    }
  lelem& operator = (const relem& rval) {
    if(usesenv(rval.singleton, singleton.ain)) throw assignment_exception();
    singleton.removeall();
    singleton |= rval.singleton;
    return *this;
    }  
  lelem& operator += (const elem& rval) {
    auto& s(singleton);
    if(usesenv(rval, s.ain)) throw assignment_exception();
    for(size_t i=0; i<s.orig->elts.size(); i++) if(isSet(s.orig->elts[i].a))
      asSet(s.orig->elts[i].a)->insert(rval, s.ain);
    for(auto a: s.removeallnonset()) {
      // asSet(a)->ain = currentcontext;
      throw loisexception();
      // a += rval;
      singleton += a;
      }
    return *this;
    }
  };

lset& operator += (lset& A, const relem& x);

// singleton construction: [r/l]num
//----------------------------------

template<class T> struct rnum;
template<class T> struct lnum;

template<class T> struct rnum {
  rsetof<T> singleton;
  rnum() : singleton(newSetOf<T>()) { }
  rnum(T t) : singleton(newSetOf<T>(t)) { }
  rnum(extracttype, rsetof<T> x) : singleton(x) { }  
  rbool isUndefined() { return singleton.isEmpty(); }  
  };

template<class T> inline rnum<T> extractnum(rsetof<T> e) { return rnum<T>(EXTRACT, e); }

template<class T> inline rbool operator == (rnum<T> left, T right) {
  return left.singleton == newSetOf<T> (right);
  }

template<class T> inline rbool operator == (rnum<T> left, rnum<T> right) {
  return left.singleton == right.singleton;
  }

template<class T> inline std::ostream& operator << (std::ostream& os, rnum<T> e) {
  return os << sym.pseudo << e.singleton;
  }

rnum<int> cardinality(rset X);

template<class T> struct lnum {
  lsetof<T> singleton;
  lnum() {}
  lnum(T t) : singleton(newSetOf<T>(t)) { }
  lnum(extracttype, rsetof<T> x) : singleton(x) { }
  rbool isUndefined() { return singleton.isEmpty(); }  
  operator rnum<T>() const { return extractnum(singleton.r()); }

  lnum<T>& operator = (const T& rval) {
    singleton.removeall();
    singleton += rval;
    return *this;
    }

  lnum<T>& operator = (const rnum<T>& rval) {
    singleton.removeall();
    singleton |= rval.singleton;
    return *this;
    }

  lnum<T>& operator += (rnum<T> B) {
    for(T a: singleton.removeall()) 
    for(T b: B.singleton) 
    for(int c: cardinality(branchset(singleton.orig.ain)).singleton)
      singleton += (a + c * b);
    return *this;
    }

  lnum<T>& operator += (T b) {
    for(T a: singleton.removeall())
    for(int c: cardinality(branchset(singleton.orig.ain)).singleton)
      singleton += (a + c * b);
    return *this;
    }

  lnum<T>& operator ++ (int) {
    (*this) += 1;
    return *this;
    }
  };

template<class T> inline rnum<T> operator + (rnum<T> A, int b) {
  lnum<T> result;
  for(T a: A.singleton) result.singleton += (a+b);
  return result;
  }

template<class T> inline rnum<T> operator + (rnum<T> A, rnum<T> B) {
  lnum<T> result;
  for(T a: A.singleton) for(T b: B.singleton) 
    result.singleton += (a+b);
  return result;
  }

template<class T> inline std::ostream& operator << (std::ostream& os, lnum<T> e) {
  return os << sym.pseudo << e.singleton;
  }

template<class T> inline rbool operator == (lnum<T> left, T right) {
  return left.singleton == newSetOf<T> (right);
  }

template<class T> inline rbool operator == (lnum<T> left, rnum<T> right) {
  return left.singleton == right.singleton;
  }

template<class T> inline rbool operator == (lnum<T> left, lnum<T> right) {
  return left.singleton == right.singleton;
  }

// evaluate a relation (given as a graph) on x
relem eval(rset f, elem x);

// is the given function injective over X?
rbool isInjective(rset f, rset X);

inline rbool imply(rbool a, rbool b) { return b || !a; }
inline rbool equivalent(rbool a, rbool b) { return imply(a,b) && imply(b,a); }

// get the orbit of 'what', with respect to the relation 'rels'; 
// contexts can be used to specify which elements are fixed
rset getorbit(elem what, std::initializer_list<Relation*> rels, contextptr nowcontext = currentcontext, contextptr anccontext = emptycontext);

// get the set of singletons of elements of X, represented as a single set-builder
// expression
rset getsingletonset(rset X);

struct EIteratorParallel {
  EIterator it;
  EIteratorParallel(EIterator&& i) : it(NULL, 0) {
    it.s = i.s; it.index = i.index; it.at = i.at; 
    it.ourcontext = i.ourcontext; i.at = nullelem;
    }
  EIteratorParallel(const EIterator& i) = delete;
  EIteratorParallel& operator ++ () { ++it; return *this; }
  relem operator* () { return extract(asSet(*it)); }
  EIteratorParallel(EIteratorParallel&& src) : it(NULL, 0) { 
    it.s = src.it.s; it.index = src.it.index; it.at = src.it.at; 
    it.ourcontext = src.it.ourcontext; src.it.at = nullelem;
    }
  bool operator != (const EIteratorParallel& x) { return it != x.it; }
  };

struct fullypseudoparallel {
  rset over;
  fullypseudoparallel(rset o) : over(getsingletonset(o)) {}
  EIteratorParallel begin() { return EIteratorParallel(over->begin()); }
  EIteratorParallel end() { return EIteratorParallel(over->end()); }
  };

// use the namespace orderedfield_ops to take the +, -, *, / from the structure
// given as mainField (RelReal or RelInt); also integer constants taken as LOIS elements
// will be interpreted as elements of mainDomain / mainField instead of the usual C++
// ints. When using this, you would probably also want to change mainOrder.

namespace orderedfield_ops {
  extern RelOrderedField *mainField;
  extern Domain *mainDomain;

  inline term operator+ (elem x, elem y) { return mainField->plus(as<term>(x),as<term>(y)); }
  inline term operator- (elem x, elem y) { return mainField->minus(as<term>(x),as<term>(y)); }
  inline term operator* (elem x, elem y) { return mainField->times(as<term>(x),as<term>(y)); }
  inline term operator/ (elem x, elem y) { return mainField->divide(as<term>(x),as<term>(y)); }
  }

}

#endif

