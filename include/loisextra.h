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
inline int substitute(int name, const varsubstlist& l) { return name; }

// string constant

inline bool isused(vptr v, const std::string& name) { return false; }
inline std::string substitute(const std::string& name, const varsubstlist& l) { return name; }

// statically typed pair
// std::pair assumes that the result of comparison is of type 'bool',
// and that's why we cannot use it

template<class T, class U> struct lpair {
  T first;
  U second;
  lpair() {}
  lpair(T f, U s) : first(f), second(s) {}
  };

template<class T, class U> lpair<T,U> make_lpair(T t, U u) { return lpair<T,U>(t,u); }

template<class T, class U> 
  bool isused(vptr v, lpair<T,U> p) { return isused(v, p.first) || isused(v, p.second); }

template<class T, class U>  inline lpair<T,U> substitute(const lpair<T,U>& p, const varsubstlist& l) { 
  return make_lpair(substitute(p.first, l), substitute(p.second, l));
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

// statically typed tuple
// std::vector assumes that the result of comparison is of type 'bool',
// and that's why we cannot use it

template<class T> struct lvector {
  std::vector<T> inner;
  lvector(const std::vector<T>& a) : inner(a) {}
  lvector(std::vector<T>&& a) : inner(a) {}
  lvector() {}
  lvector(const std::initializer_list<T>& l) : inner(l) {}
  void push_back(const T& t) { inner.push_back(t); }
  typename std::vector<T>::const_iterator begin() const { return inner.begin(); }
  typename std::vector<T>::const_iterator end() const { return inner.end(); }
  typename std::vector<T>::iterator begin() { return inner.begin(); }
  typename std::vector<T>::iterator end() { return inner.end(); }
  size_t size() const { return inner.size(); }
  T& operator [] (int i) { return inner[i]; }
  const T& operator [] (int i) const { return inner[i]; }
  };

template<class T>  bool isused(vptr v, lvector<T> p) { 
  for(auto& el: p.inner) if(isused(v, el)) return true;
  return false;
  }

template<class T>  lvector<T> substitute(const lvector<T>& p, const varsubstlist& l) { 
  lvector<T> res;
  for(auto& a: p) res.push_back(substitute(a,l));
  return res;
  }

template<class T> 
  rbool operator == (const lvector<T>& a, const lvector<T>& b) { 
    lbool res = ftrue;
    if(a.size() != b.size()) { res = ffalse; return res;}
    for(int i=0; i<(int) a.size(); i++)
      res &= (a[i] == b[i]);
    return res; 
    }
template<class T> rbool operator != (const lvector<T>& a, const lvector<T>& b) { 
  return !(a == b); }
  
template<class T> rbool operator < (const lvector<T>& a, const lvector<T>& b) { 
  lbool res = ffalse;
  lbool eq = ftrue;
  for(int i=0;; i++) {
    if(i >= (int) a.size()) return res;
    if(i >= (int) b.size()) return res || eq;
    res |= (eq && (a[i] < b[i]));
    eq &= a[i] == b[i];
    }
  }

template<class T> rbool operator <= (const lvector<T>& a, const lvector<T>& b) { 
  lbool res = ffalse;
  lbool eq = ftrue;
  for(int i=0;; i++) {
    if(i >= (int) b.size()) return res || eq;
    if(i >= (int) a.size()) return res;
    res |= (eq && (a[i] < b[i]));
    eq &= a[i] == b[i];
    }
  }

template<class T> 
  std::ostream& operator << (std::ostream& os, const lvector<T>& p) { 
    bool first = true;
    os << "[";    
    for(auto& t: p) { if(!first) os << ","; first = false; os << t; }
    os << "]";
    return os;
    }

// basic set theoretic operations
//--------------------------------

template<class T> rbool memberof(const T& x, const lsetof<T> & y) {
  return y.hasmember(x);
  }

template<class T> rbool subseteq(const lsetof<T>& x, const lsetof<T>& y) {
  return x.subseteq(y);
  }

template<class T> lsetof<T> operator & (const lsetof<T>& x, const lsetof<T>& y) { 
  lsetof<T> xy;
  for(const T& e: x) If(memberof(e, y)) xy += e;
  return xy;
  }

template<class T> lsetof<T> operator & (const lsetof<T>& x, negated<lsetof<T>> y) { 
  lsetof<T> xy;
  for(const T& e: x) If(!memberof(e, y.original)) xy += e;
  return xy;
  }

template<class T> lsetof<T> operator | (const lsetof<T>& x, const lsetof<T>& y) { 
  lsetof<T> xy;
  for(const T& e: x) xy += e;
  for(const T& e: y) xy += e;
  return xy;
  }

struct pushcontext {
  contextptr orig;
  pushcontext(contextptr p) { orig = currentcontext; currentcontext = p; }
  ~pushcontext() { currentcontext = orig; }
  };

template<class T> lsetof<T> lsetof<T>::removeall() {
  
  lsetof<T> xyes;
  rbool phi = outenv(true, false, ain);

  swap(xyes.elts, elts);

  pushcontext pc(ain);
  for(const T& a: xyes) Ife(!phi) {
    insert(a, ain);
    }

  return xyes;
  }

template<class T> lsetof<T>& lsetof<T>::operator -= (const T& y) {
  lsetof<T> x1 = removeall();
  insert(y, ain);
  lsetof<T> y1 = removeall();
  for(const T& e: x1) If(!memberof(e, y1)) insert(e, ain);
  return (*this);
  }

template<class T> lsetof<T>& lsetof<T>::operator &= (const negated<lsetof<T>>& y) {
  lsetof<T> x1 = removeall();
  (*this) |= y.original;
  lsetof<T> y1 = removeall();
  for(const T& e: x1) If(!memberof(e, y1)) insert(e, ain);
  return (*this);
  }

template<class T> lsetof<T>& lsetof<T>::operator &= (const lsetof<T>& y) {
  // does not work that easily:
  // for(elem e: removeall(x)) If(e != y) x += e;
  lsetof<T> x1 = removeall();
  (*this) += y;
  lsetof<T> y1 = removeall();
  for(const T& e: x1) {
    lbool f(ftrue);
    for(const T& ye: y1) {
      f &= memberof(e, ye);
      }
    If(f) (*this) += e;
    }
  return (*this);
  }

// remove the repetitions: each element will be listed only once in the optimized set
template<class T> lsetof<T> optimize(const lsetof<T>& x) {
  lsetof<T> y;
  for(const T& e: x) If(!memberof(e, y)) y += e;
  return y;
  }

// optimize the elements of given type
template<class T> lsetof<T> optimizeType(const lsetof<T>& x, const lsetof<T>& type) {
  return (type & x) | (x &~ type);
  }

// the Cartesian product
template<class T, class U>
lsetof<lpair<T,U>> operator * (const lsetof<T> & x, const lsetof<U> & y) {
  lsetof<lpair<T,U>> xy;
  for(const T& e: x) for(const U& f: y) xy += make_lpair(e,f);
  return xy;
  }

// the Cartesian product of more than 2
// lset cartesian(std::initializer_list<elem> l, eltuple v = eltuple());

// the Cartesian power
// template<class T> lsetof<std::vector<T>> cartesianPower(const lsetof<T>& x, int power, std::vector<T> v = std::vector<T>());

// macros
//--------

#define FORALL(x, set, condition) \
  ([&]() -> rbool { lbool ret(ftrue); for(auto x: (set)) ret &= (condition); return ret; }())
#define EXISTS(x, set, condition) \
  ([&]() -> rbool { lbool ret(ffalse); for(auto x: (set)) ret |= (condition); return ret; }())
#define FILTERMAP(x, set, condition, replace, T) \
  ([&]() -> lsetof<T> { lsetof<T> S; for(auto x: (set)) If(condition) S += (replace); return S; }())
#define FILTER(x, set, condition, T) FILTERMAP(x, set, condition, x, T)
#define MAP(x, set, replace, T) FILTERMAP(x, set, ftrue, replace, T)

// axioms & declareatom
//----------------------

// declare an axiom
struct axiom : ArbCondition {
  axiom(rbool phi);
  };

// declare an atom
struct declareatom : ArbCondition {
  term at;
  declareatom(Domain *d, const std::string& s);
  operator term() { return at; }
  };

// singleton construction: [r/l]elem
//-----------------------------------

extern struct extracttype {} EXTRACT;

template<class T> struct lelemof {
  lsetof<T> singleton;
  lelemof() {}
  lelemof(const T& b) { singleton |= newSet(b); }
  lelemof(extracttype, const lsetof<T>& X) : singleton(X) {}
  rbool isUndefined() { return isempty(singleton); }  

  lelemof(const lelemof<T>& b) { singleton |= b.singleton; }
  
  lelemof<T>& operator = (const T& rval) {
    if(usesenv(rval, singleton.ain)) throw assignment_exception();
    singleton.removeall();
    singleton += rval;
    return *this;
    }
  lelemof<T>& operator = (const lelemof<T>& rval) {
    if(usesenv(rval.singleton, singleton.ain)) throw assignment_exception();
    singleton.removeall();
    singleton |= rval.singleton;
    return *this;
    }
  
  EIteratorOf<T> begin() const { return EIteratorOf<T> (singleton, 0); }
  EIteratorOf<T> end() const { return EIteratorOf<T> (singleton, singleton.elts.size()); }

  /* lelemof<T>& operator += (const T& rval) {
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
    } */
  };

template<class T> inline lbool memberof(const lelemof<T>& x, const lsetof<T>& A) {
  return subseteq(x.singleton, A);
  }

template<class T> lsetof<T>& operator += (lsetof<T>& A, const lelemof<T>& x);

// extract the element from a singleton set, and vice versa
template<class T> lelemof<T> extract(const lsetof<T>& e) { return lelemof<T>(EXTRACT, e); }
template<class T> const lsetof<T>& newSet(const lelemof<T>& el) { return el.singleton; }

template<class T> inline rbool operator == (const lelemof<T>& a, const lelemof<T>& b) 
  { return a.singleton == b.singleton; }

template<class T> inline rbool operator != (const lelemof<T>& a, const lelemof<T>& b) 
  { return a.singleton != b.singleton; }

template<class T> inline rbool operator == (const lelemof<T>& a, const T& b) 
  { return a.singleton == newSet(b); }

template<class T> inline rbool operator != (const lelemof<T>& a, const T& b) 
  { return a.singleton != newSet(b); }

template<class T> rbool operator > (const lelemof<T>& a, const lelemof<T>& b) {
  lbool p = true;
  for(const T& x: a.singleton) for(const T& y: b.singleton) p &= x>y;
  return p;
  }

template<class T> inline std::ostream& operator << (std::ostream& os, const lelemof<T>& e) {
  return os << sym.pseudo << e.singleton;
  }

// the set of extensions of anccontext-valuation to nowcontext
lsetof<lvector<term>> branchset(contextptr anccontext = emptycontext, contextptr nowcontext = currentcontext);

template<class T> inline lelemof<int> cardinality(const lsetof<T>& X);

inline lelemof<int>& operator += (lelemof<int>& v, const lelemof<int>& x) {
  for(int a: v.singleton.removeall()) 
  for(int b: x.singleton) 
  for(int c: cardinality(branchset(v.singleton.ain)).singleton)
    v.singleton += (a + c * b);
  return v;
  }

inline lelemof<int>& operator += (lelemof<int>& v, int b) {
  for(int a: v.singleton.removeall())
  for(int c: cardinality(branchset(v.singleton.ain)).singleton)
    v.singleton += (a + c * b);
  return v;
  }

inline lelemof<int>& operator ++ (lelemof<int>& v, int) {
  const int i = 1;
  v += i;
  return v;
  }

inline lelemof<int> operator + (const lelemof<int>& A, int b) {
  lelemof<int> result;
  for(int a: A.singleton) result.singleton += (a+b);
  return result;
  }

inline lelemof<int> operator + (const lelemof<int>& A, const lelemof<int>& B) {
  lelemof<int> result;
  for(int a: A.singleton) for(int b: B.singleton) 
    result.singleton += (a+b);
  return result;
  }

template<class T> inline lelemof<int> cardinality(const lsetof<T>& X) {
  lelemof<int> answer;
  for(auto a: X) If(answer.isUndefined()) {
    auto butone = cardinality(X &~ newSet(a));
    for(int x: butone.singleton) answer.singleton += (x+1);
    }
  If(answer.isUndefined()) answer.singleton += 0;
  return answer;
  }

// evaluate a relation (given as a graph) on x
template<class T, class U> lelemof<U> eval(const lsetof<lpair<T,U>>& f, const T& x) {
  lsetof<U> ret;
  for(auto& p: f) If(x == p.first) ret += p.second;
  return extract(ret);
  }

// is the given function injective over X?
template<class T, class U>
rbool isInjective(const lsetof<lpair<T,U>>& f, const lsetof<T>& X) {
  lbool ret(ftrue);
  for (auto& x:X) for(auto& y:X) If(x!=y) {
    ret &= (eval(f,x) != eval(f,y));
    }
  return ret;
  }

inline rbool imply(rbool a, rbool b) { return b || !a; }
inline rbool equivalent(rbool a, rbool b) { return imply(a,b) && imply(b,a); }

// get the orbit of 'what', with respect to the relation 'rels'; 
// contexts can be used to specify which elements are fixed
template <class T> lsetof<T> getorbit(T what, std::initializer_list<Relation*> rels, contextptr nowcontext = currentcontext, contextptr anccontext = emptycontext);

// get the set of singletons of elements of X, represented as a single set-builder
// expression
template <class T> lsetof<lsetof<T>> getsingletonset(const lsetof<T>& X);

template<class T> struct EIteratorParallel {
  EIteratorOf<T> it;
  EIteratorParallel(EIteratorOf<T>&& i) : it(NULL, 0) {
    it.s = i.s; it.index = i.index; it.at = i.at; 
    it.ourcontext = i.ourcontext; i.at = NULL;
    }
  EIteratorParallel(const EIteratorParallel<T>& i) = delete;
  EIteratorParallel& operator ++ () { ++it; return *this; }
  lelemof<T> operator* () { return extract(*it); }
  EIteratorParallel(EIteratorParallel<T>&& src) : it(NULL, 0) { 
    it.s = src.it.s; it.index = src.it.index; it.at = src.it.at; 
    it.ourcontext = src.it.ourcontext; src.it.at = NULL;
    }
  bool operator != (const EIteratorParallel<T>& x) { return it != x.it; }
  };

template<class T> struct fullypseudoparallel {
  lsetof<T> over;
  fullypseudoparallel(lsetof<T> o) : over(getsingletonset(o)) {}
  EIteratorParallel<T> begin() { return EIteratorParallel<T>(over->begin()); }
  EIteratorParallel<T> end() { return EIteratorParallel<T>(over->end()); }
  };

// use the namespace orderedfield_ops to take the +, -, *, / from the structure
// given as mainField (RelReal or RelInt); also integer constants taken as LOIS elements
// will be interpreted as elements of mainDomain / mainField instead of the usual C++
// ints. When using this, you would probably also want to change mainOrder.

namespace orderedfield_ops {
  extern RelOrderedField *mainField;
  extern Domain *mainDomain;

  inline term operator+ (term x, term y) { return mainField->plus(x,y); }
  inline term operator- (term x, term y) { return mainField->minus(x,y); }
  inline term operator* (term x, term y) { return mainField->times(x,y); }
  inline term operator/ (term x, term y) { return mainField->divide(x,y); }
  }

}

#endif

