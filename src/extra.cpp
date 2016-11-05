#include "../include/loisextra.h"

// needed for phi->alph in getsingletonset
#include "../include/loisinternal.h"


namespace lois {

inline rbool operator == (eltuple a, eltuple b) { 
  if(a.size() != b.size()) return ffalse;
  rbool f = ftrue;
  for(size_t i=0; i<a.size(); i++) if(!f.isFalse()) f = f && a[i] == b[i];
  return f;
  }

rbool fless(eltuple a, eltuple b, size_t at) {
  if(at >= b.size()) return ffalse;
  if(at >= a.size()) return ftrue;
  return a[at] < b[at] || (a[at] == b[at] && fless(a, b, at+1));
  }

rbool flesseq(eltuple a, eltuple b, size_t at) {
  if(at >= a.size() && at >= b.size()) return ftrue;
  if(at >= b.size()) return ffalse;
  if(at >= a.size()) return ftrue;
  return a[at] < b[at] || (a[at] == b[at] && flesseq(a, b, at+1));
  }

rbool operator <  (eltuple a, eltuple b) { return fless(a, b, 0); }
rbool operator <= (eltuple a, eltuple b) { return flesseq(a, b, 0); }

std::ostream& operator << (std::ostream& os, eltuple t) { 
  os << "[";
  for(size_t i=0; i<t.size(); i++) {
    if(i) os << ", ";
    os << t[i];
    }
  return os << "]";
  }

elem::elem(eltuple x) : p(elem(elof<eltuple> (x)).p) {}
elem::elem(elpair x) : p(elem(elof<elpair> (x)).p) {}
elem::elem(std::string x) { p = elem(elof<std::string> (x)).p; }
elem::elem(int x) { 
  using namespace orderedfield_ops;
  if(mainField && mainDomain)
    p = elem(elof<term> (mainField->constant(mainDomain, x))).p;
  else
    p = elem(elof<int> (x)).p;
  }
elem::elem(term x) : p(elem(elof<term> (x)).p) {}
elem::elem(vptr x) : p(elem(elof<term> (term(x))).p) {}

namespace orderedfield_ops {
  RelOrderedField *mainField;
  Domain *mainDomain;
  }

struct pushcontext {
  contextptr orig;
  pushcontext(contextptr p) { orig = currentcontext; currentcontext = p; }
  ~pushcontext() { currentcontext = orig; }
  };


lset lset::removeall() {
  
  lset xyes  = std::make_shared<ESet> ();
  rbool phi = outenv(true, false, ain);

  swap(xyes->elts, p->elts);

  pushcontext pc(ain);
  for(elem a: xyes) Ife(!phi) {
    p->insert(a, ain);
    }

  return xyes;
  }

elem removeallnonset(elem x) {
  /*
  lset xs(x);
  auto xset  = std::make_shared<ESet> ();

  auto xyes  = std::make_shared<ESet> ();
  rbool phi = outenv(true, false, xs->ain);
  
  for(size_t i=0; i<xs->elts.size(); i++) 
    if(isSet(xs->elts[i].a))
      xset->elts.emplace_back(xs->elts[i].var, xs->elts[i].phi, xs->elts[i].a);
    else
      xyes->elts.emplace_back(xs->elts[i].var, xs->elts[i].phi, xs->elts[i].a);

  swap(xset->elts, xs->elts);

  env = xs->ain;

  for(elem a: xyes) Ife(!phi) {
    xs->insert(a);
    }

  env = xyes->ain;

  return elem(xyes); */
  throw loisexception();
  }

lset& lset::operator -= (elem y) {
  lset x1 = removeall();
  p->insert(y, ain);
  lset y1 = removeall();
  for(elem e: x1) If(!memberof(e, y1)) p->insert(e, ain);
  return (*this);
  }

lset& lset::operator &= (negated<lset> y) {
  lset x1 = removeall();
  (*this) |= y.original;
  lset y1 = removeall();
  for(elem e: x1) If(!memberof(e, y1)) p->insert(e, ain);
  return (*this);
  }

lset& lset::operator &= (const lset& y) {
  // does not work that easily:
  // for(elem e: removeall(x)) If(e != y) x += e;
  lset x1 = removeall();
  (*this) += y;
  lset y1 = removeall();
  for(elem e: x1) {
    lbool f(ftrue);
    for(elem ye: y1) {
      f &= memberof(e, asSet(ye));
      }
    If(f) (*this) += e;
    }
  return (*this);
  }

// remove the repetitions: each element will be listed only once in the optimized set
lset optimize(const lset& x) {
  lset y;
  for(elem e: x) If(!memberof(e, y)) y += e;
  return y;
  }

// optimize the elements of given type
lset optimizeType(const lset& x, const lset& type) {
  return (type & x) | (x &~ type);
  }

// the Cartesian product
lset operator * (const lset& x, const lset& y) {
  lset xy;
  for(elem e: x) for(elem f: y) xy += elpair(e,f);
  return xy;
  }

// the Cartesian product of more than 2
template<class T> void cartesianRange(T beg, T end, eltuple& v, lset& cart) {
  if(beg == end) cart += v;
  else for(elem a: *beg) {
    v.push_back(a);
    beg++;
    cartesianRange(beg, end, v, cart);
    v.pop_back();
    beg--;
    }
  }

lset cartesian(std::initializer_list<elem> l, eltuple v) {
  lset cart;
  cartesianRange(l.begin(), l.end(), v, cart);
  return cart;
  }

axiom::axiom(rbool phi) {
  if(!enable(phi)) throw unsatisfiable_exception();
  }

declareatom::declareatom(Domain *d, const std::string& s) {
  vptr v = newvar(d, s);
  at = v;
  if(!enable(varlist({v}), ftrue)) 
    throw unsatisfiable_exception();
  }

lset& operator += (lset& A, const lelem& x) {
  for(elem e: newSet(x)) A += e;
  return A;
  }

// getorbit
lset getorbit(elem what, std::initializer_list<Relation*> rels, contextptr nowcontext, contextptr anccontext) {

  lset X = newSet();
  
  std::vector<std::pair<vptr, vptr>> v;
  rbool phi = true;

  contextptr e = nowcontext;

  while(e != anccontext) {
    if(!e) throw env_exception();
    for(auto w: e->var) {
      v.push_back(make_pair(w, newvar(w->dom)));
      }
    phi = phi && e->phi;
    e = e->parent;
    }
  
  for(auto p: v) {
    // phi = phi->alph(p.first, p.second);
    what = alpha(what, p.first, p.second);
    }
  
  for(int i=0; i<v.size(); i++) for(int j=0; j<i; j++) {
    term ai = term(v[i].first);
    term aj = term(v[j].first);
    term bi = term(v[i].second);
    term bj = term(v[j].second);
    phi = phi && equivalent(ai==aj, bi==bj);
    }
    
  for(auto r: rels) phi = phi && r->isomorphic(v);
    
  e = std::make_shared<Context> (currentcontext, phi);
  for(auto p: v) e->var.push_back(p.second);

  X->insert(what, currentcontext, e);
  
  return X;
  }

// singleton

Domain *singletonDomain;
RelUnary *singletonRelation;

std::string singletonRelationName = "SR";

lset getsingletonset(const lset& X) {

  if(!singletonDomain) {
    singletonDomain = new Domain("S");
    singletonRelation = new RelUnary(singletonRelationName, 2);
    }
    
  lset Result = newSet();
  elem elY = newSet();
  std::shared_ptr<ESet> Y = std::dynamic_pointer_cast<ESet> (elY.p);
  std::vector<vptr> v;
  vptr v0 = newvar(singletonDomain);
  
  contextptr eorig = currentcontext;
  
  int id = 0;
  for(elem x: X) {
    std::cout << "x = " << x << std::endl;
    
    std::vector<vptr> vx;
    contextptr e = currentcontext;
    
    if(singletonRelation->nogroups <= id)
      singletonRelation->nogroups++;
  
    rbool phi = (*singletonRelation) (v0, id);
    id++;
    while(e != eorig) {
      if(!e) throw env_exception();
      for(auto w: e->var) {
        std::cout << "add var " << w << std::endl;
        vx.push_back(w);
        }
      std::cout << "add phi " << e->phi << std::endl;
      phi = phi && e->phi;
      e = e->parent;
      }

    std::vector<bool> used(v.size(), false);
    
    for(vptr wx: vx) {
      vptr w;
      for(int i=0; i<int(v.size()); i++) 
        if(v[i]->dom == wx->dom && !used[i]) {
          w = v[i]; used[i] = true; goto found;
          }
      w = newvar(wx->dom);
      std::cout << "add v-var " << w << std::endl;
      v.push_back(w); used.push_back(true);
      found:
      std::cout << "translate " << wx << " to " << w << std::endl;
      x = alpha(x, wx, w);
      phi = phi->alph(wx, w);
      }
    
    Y->elts.push_back(SimpleSet(varlist(), phi, x));
    }

  // std::cout << "Y is " << Y << std::endl;
  SimpleSet s(varlist(), true, elY);
  std::cout << "s loaded" << std::endl;
  s.var.push_back(v0);
  for(vptr w: v) s.var.push_back(w);
  Result->elts.push_back(std::move(s));
    
  return Result;
  }

}
