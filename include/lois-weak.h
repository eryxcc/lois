// this file defines the struct 'elem' for weakly typed sets

#ifndef _lois_weak_h_
#define _lois_weak_h_

#include <memory>
#include "loisextra.h"

namespace lois {

struct elem {
  std::shared_ptr<struct Element> p; 
  elem(std::shared_ptr<struct Element> e) : p(e) {}
  elem() {}
  Element* operator -> () const { return &(*p); }
//  elem(int);
//  elem(std::string);
//  elem(vptr);
//  elem(term);
  };

extern const elem nullelem;

// elem is internally a shared_ptr to Element

struct Element {
  // destructor
  virtual ~Element() {}
  // display on a stream
  virtual std::ostream& display (std::ostream &os) const = 0;
  // alpha-convert this from 'v1' to 'v2' (i.e., substitute variable)
  virtual elem subst(const varsubstlist& l) const = 0;
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

inline std::ostream& operator << (std::ostream& os, elem a) { os << "[."; a->display(os); return os << ".]"; }

template<class T> struct ElementOf: Element {
  T data;
//  explicit ElementOf(T&& tmp, int) : data(tmp) {}
//  explicit ElementOf(const T& tmp, int) : data(tmp) {}
  explicit ElementOf(T&& tmp) : data(tmp) {}
  explicit ElementOf(const T& tmp) : data(tmp) {}
  std::ostream& display (std::ostream &os) const { return os << data; }
  virtual elem subst(const varsubstlist& l) const { 
    return elem(std::make_shared<ElementOf> (substitute(data, l)));
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

template<class T> elem elof(T x) { return elem(std::make_shared<ElementOf<T>> (x)); }
// template<class T> elem elof(T&& x) { return elem(std::make_shared<ElementOf<T>> (x)); }

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

inline elem substitute(elem x, const varsubstlist& l) { return x.p->subst(l); }

typedef lsetof<elem> lset;
typedef lelemof<elem> lelem;

}

#endif
