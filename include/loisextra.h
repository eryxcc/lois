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

// string constant

inline bool isused(vptr v, const std::string& name) { return false; }
inline std::string alpha(const std::string& name, vptr v1, vptr v2) { return name; }

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

template<class T> rbool operator == (const lsetof<T>& a, const lsetof<T>& b) { return a.orig == b.orig; }
template<class T> rbool operator != (const lsetof<T>& a, const lsetof<T>& b) { return a.orig != b.orig; }

// basic set theoretic operations
//--------------------------------

inline rbool memberof(elem x, const lset& y) {
  lbool af(ffalse); for(elem a: y) af |= (x == a); return af;
  }

template<class T> inline rbool memberof(T x, const lsetof<T>& y) {
  lbool af(ffalse); for(T a: y) af |= (x == a); return af;
  }

inline rbool subseteq(const lset& x, const lset& y) {
  lbool af(ftrue); for(elem a: x) af &= memberof(a, y); return af;
  }

inline lset operator & (const lset& x, const lset& y) { 
  lset xy;
  for(elem e: x) If(memberof(e, y)) xy += e;
  return xy;
  }

inline lset operator & (const lset& x, negated<lset> y) { 
  lset xy;
  for(elem e: x) If(!memberof(e, y.original)) xy += e;
  return xy;
  }

inline lset operator | (const lset& x, const lset& y) { 
  lset xy;
  for(elem e: x) xy += e;
  for(elem e: y) xy += e;
  return xy;
  }

template<class T> inline lsetof<T> operator | (const lsetof<T>& x, const lsetof<T>& y) { 
  lsetof<T> xy;
  for(elem e: x) xy += e;
  for(elem e: y) xy += e;
  return xy;
  }

// remove the repetitions: each element will be listed only once in the optimized set
lset optimize(const lset& x);

// optimize the elements of given type
lset optimizeType(const lset& x, const lset& type);

// the Cartesian product
lset operator * (const lset& x, const lset& y);

// the Cartesian product of more than 2
lset cartesian(std::initializer_list<elem> l, eltuple v = eltuple());

// the Cartesian power
lset cartesianPower(const lset& x, int power, eltuple v = eltuple());

// macros
//--------

#define FORALL(x, set, condition) \
  ([=]() -> rbool { lbool ret(ftrue); for(auto x: (set)) ret &= (condition); return ret; }())
#define EXISTS(x, set, condition) \
  ([=]() -> rbool { lbool ret(ffalse); for(auto x: (set)) ret |= (condition); return ret; }())
#define FILTERMAP(x, set, condition, replace) \
  ([=]() -> lset { lset S; for(auto x: (set)) If(condition) S += (replace); return S; }())
#define FILTER(x, set, condition) FILTERMAP(x, set, condition, x)
#define MAP(x, set, replace) FILTERMAP(x, set, ftrue, replace)
    
#define SETOF_1(expr, decl1, condition) \
  ([=]() -> lset { lset S; for(auto decl1) If(condition) S += (expr); return S; }())
#define SETOF_2(expr, decl1, decl2, condition) \
  ([=]() -> lset { lset S; for(auto decl1) for(auto decl2) If(condition) S += (expr); return S; }())
#define SETOF_3(expr, decl1, decl2, decl3, condition) \
  ([=]() -> lset { lset S; for(auto decl1) for(auto decl2) for(auto decl3) If(condition) S += (expr); return S; }())
#define SETOF_4(expr, decl1, decl2, decl3, decl4, condition) \
  ([=]() -> lset { lset S; for(auto decl1) for(auto decl2) for(auto decl3) for(auto decl4) If(condition) S += (expr); return S; }())
#define SETOF_5(expr, decl1, decl2, decl3, decl4, decl5, condition) \
  ([=]() -> lset { lset S; for(auto decl1) for(auto decl2) for(auto decl3) for(auto decl4) for(auto decl5) If(condition) S += (expr); return S; }())

#define GET_MACRO(U,  _1,_2,_3,_4,_5, NAME,...) NAME
#define SETOF(expr, ...) GET_MACRO(__VA_ARGS__, SETOF_5,SETOF_4, SETOF_3, SETOF_2, SETOF_1)(expr, __VA_ARGS__)
    

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

struct lelem {
  lset singleton;
  lelem() {}
  lelem(elem b) { singleton |= newSet(b); }
  lelem(extracttype, const lset& X) : singleton(X) {}
  rbool isUndefined() { return singleton == emptyset; }  

  lelem(const lelem& b) { singleton |= b.singleton; }
  
  lelem& operator = (const elem& rval) {
    if(usesenv(rval, singleton.ain)) throw assignment_exception();
    singleton.removeall();
    singleton += rval;
    return *this;
    }
  lelem& operator = (const lelem& rval) {
    if(usesenv(rval.singleton, singleton.ain)) throw assignment_exception();
    singleton.removeall();
    singleton |= rval.singleton;
    return *this;
    }  
  lelem& operator += (const elem& rval) {
    auto& s(singleton);
    if(usesenv(rval, s.ain)) throw assignment_exception();
    for(size_t i=0; i<s.p->elts.size(); i++) if(isSet(s.p->elts[i].a))
      asSet(s.p->elts[i].a)->insert(rval, s.ain);
    for(auto a: s.removeallnonset()) {
      // asSet(a)->ain = currentcontext;
      throw loisexception();
      // a += rval;
      singleton += a;
      }
    return *this;
    }
  };

inline lbool memberof(const lelem& x, const lset& A) {
  return subseteq(x.singleton, A);
  }

lset& operator += (lset& A, const lelem& x);

// extract the element from a singleton set, and vice versa
inline lelem extract(const lset& e) { return lelem(EXTRACT, e); }
inline const lset& newSet(const lelem& el) { return el.singleton; }

inline rbool operator == (const lelem& a, const lelem& b) 
  { return a.singleton == b.singleton; }

inline rbool operator != (const lelem& a, const lelem& b) 
  { return a.singleton != b.singleton; }

inline rbool operator > (const lelem& a, const lelem& b) {
  lbool p = true;
  for(elem x: a.singleton) for(elem y: b.singleton) p &= x>y;
  return p;
  }

inline std::ostream& operator << (std::ostream& os, const lelem& e) {
  return os << sym.pseudo << e.singleton;
  }

// singleton construction: [r/l]num
//----------------------------------

template<class T> struct lnum;

template<class T> inline lnum<T> extractnum(const lsetof<T>& e) { return lnum<T>(EXTRACT, e); }

template<class T> inline rbool operator == (const lnum<T>& left, T right) {
  return left.singleton == newSetOf<T> (right);
  }

template<class T> inline rbool operator == (const lnum<T>& left, const lnum<T>& right) {
  return left.singleton == right.singleton;
  }

template<class T> inline std::ostream& operator << (std::ostream& os, const lnum<T>& e) {
  return os << sym.pseudo << e.singleton;
  }

lnum<int> cardinality(const lset& X);

template<class T> struct lnum {
  lsetof<T> singleton;
  lnum() {}
  lnum(T t) : singleton(newSetOf<T>(t)) { }
  lnum(extracttype, const lsetof<T>& x) : singleton(x) { }
  rbool isUndefined() { return singleton.isEmpty(); }  
  operator lnum<T>() const { return extractnum(singleton); }

  lnum<T>& operator = (const T& rval) {
    singleton.removeall();
    singleton += rval;
    return *this;
    }

  lnum<T>& operator = (const lnum<T>& rval) {
    singleton.removeall();
    singleton |= rval.singleton;
    return *this;
    }

  lnum<T>& operator += (const lnum<T>& B) {
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

template<class T> inline lnum<T> operator + (const lnum<T>& A, int b) {
  lnum<T> result;
  for(T a: A.singleton) result.singleton += (a+b);
  return result;
  }

template<class T> inline lnum<T> operator + (const lnum<T>& A, const lnum<T>& B) {
  lnum<T> result;
  for(T a: A.singleton) for(T b: B.singleton) 
    result.singleton += (a+b);
  return result;
  }

// evaluate a relation (given as a graph) on x
lelem eval(const lset& f, elem x);

// is the given function injective over X?
rbool isInjective(const lset& f, const lset& X);

inline rbool imply(rbool a, rbool b) { return b || !a; }
inline rbool equivalent(rbool a, rbool b) { return imply(a,b) && imply(b,a); }

// get the orbit of 'what', with respect to the relation 'rels'; 
// contexts can be used to specify which elements are fixed
lset getorbit(elem what, std::initializer_list<Relation*> rels, contextptr nowcontext = currentcontext, contextptr anccontext = emptycontext);

// get the set of singletons of elements of X, represented as a single set-builder
// expression
lset getsingletonset(const lset& X);

struct EIteratorParallel {
  EIterator it;
  EIteratorParallel(EIterator&& i) : it(NULL, 0) {
    it.s = i.s; it.index = i.index; it.at = i.at; 
    it.ourcontext = i.ourcontext; i.at = nullelem;
    }
  EIteratorParallel(const EIterator& i) = delete;
  EIteratorParallel& operator ++ () { ++it; return *this; }
  lelem operator* () { return extract(asSet(*it)); }
  EIteratorParallel(EIteratorParallel&& src) : it(NULL, 0) { 
    it.s = src.it.s; it.index = src.it.index; it.at = src.it.at; 
    it.ourcontext = src.it.ourcontext; src.it.at = nullelem;
    }
  bool operator != (const EIteratorParallel& x) { return it != x.it; }
  };

struct fullypseudoparallel {
  lset over;
  fullypseudoparallel(lset o) : over(getsingletonset(o)) {}
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

// statically typed pair
// std::pair assumes that the result of comparison is of type 'bool',
// and that's why we cannot use it

template<class T, class U> struct lpair {
  T first;
  U second;
  lpair(T f, U s) : first(f), second(s) {}
  operator elem() const { return elem(elof<lpair<T,U>> (*this)); }
  };

template<class T, class U> lpair<T,U> make_lpair(T t, U u) { return lpair<T,U>(t,u); }

template<class T, class U> 
  bool isused(vptr v, lpair<T,U> p) { return isused(v, p.first) || isused(v, p.second); }

template<class T, class U>  inline lpair<T,U> alpha(lpair<T,U> p, vptr v1, vptr v2) { 
  return make_lpair(alpha(p.first, v1, v2), alpha(p.second, v1, v2));
  }

template<class T, class U> 
  rbool operator == (lpair<T,U> a, lpair<T,U> b) { return a.first == b.first && a.second == b.second; }
template<class T, class U> 
  rbool operator != (lpair<T,U> a, lpair<T,U> b) { return a.first != b.first || a.second != b.second; }
template<class T, class U> 
  rbool operator <  (lpair<T,U> a, lpair<T,U> b) { return a.first < b.first || (a.first == b.first && a.second < b.second); }
template<class T, class U> 
  rbool operator <= (lpair<T,U> a, lpair<T,U> b) { return a.first < b.first || (a.first == b.first && a.second <= b.second); }

template<class T, class U> 
  std::ostream& operator << (std::ostream& os, lpair<T,U> p) { return os << "(" << p.first << "," << p.second << ")"; }

}

#endif

