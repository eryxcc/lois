#include "../include/loisextra.h"

// needed for phi->alph in getsingletonset
#include "../include/loisinternal.h"
namespace lois {

namespace orderedfield_ops {
  RelOrderedField *mainField;
  Domain *mainDomain;
  }


/*
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
*/

/*
void cartesianPowerAux(const lset& x, eltuple& v, int left, lset& cart) {
  if(!left) cart += v;
  else for(elem a: x) {
    v.push_back(a);
    cartesianPowerAux(x, v, left-1, cart);
    v.pop_back();
    }
  }

lset cartesianPower(const lset& x, int power, eltuple v) {
  lset cart;
  cartesianPowerAux(x, v, power, cart);
  return cart;
  }
*/

axiom::axiom(rbool phi) {
  if(!enable(phi)) throw unsatisfiable_exception();
  }

declareatom::declareatom(Domain *d, const std::string& s) {
  vptr v = newvar(d, s);
  at = term(v);
  if(!enable(varlist({v}), ftrue)) 
    throw unsatisfiable_exception();
  }

/*
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
*/

}
