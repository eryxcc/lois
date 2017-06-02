// definitions used internally by LOIS -- most users need not read that

#ifndef _loisinternal_h_
#define _loisinternal_h_

#include "lois.h"

#include <set>

namespace lois {

struct SubRelation {
  Relation *rel;
  SubRelation(Relation *r) : rel(r) {}
  virtual bool exists(struct FormulaQ* q, int at, std::shared_ptr<struct Subdomain>& dom, int val, int rid) = 0;
  virtual bool check(int id, int val1, int val2) = 0;
  virtual int checkterm(int id, int val1, int val2) { 
    throw noterms_exception(); 
    }
  };

struct Subdomain {
  int nextval;
  std::vector<subrelptr> rels;
  int numEq, numQuantifier;
  Subdomain() { numEq = numQuantifier = nextval = 0; }
  void useRelation(Relation *r);
  };

struct featurelist {
  std::set<Domain*> domains;
  std::set<Relation*> relations;
  };

struct Formula { 
  // list free variables
  // varlist fv;
#ifdef AGGSYM
  int qty[2];
#endif
  // display the rbool
  virtual std::ostream& display (std::ostream &os) const = 0;
  virtual std::ostream& display (std::ostream &os, int prio) { return display(os); }
  // substitute 'v' for 't' in this formula
  virtual rbool subst(const varsubstlist& l) const = 0;
  // negate the rbool
  virtual rbool negate() = 0;
  // verify whether the rbool is true.
  // shared variables have values given by their 'value' property, 0<=value<nextval
  virtual bool verify() = 0;
  // if the variable 'v' used in this rbool?
  virtual bool uses(vptr v) = 0;
  virtual term valueKnown(vptr v, bool negated) { return nullterm; }
  virtual void initDomains() = 0;
  virtual void listFeatures(featurelist&) = 0;
#ifdef AGGSYM
  virtual void clearinfo() { qty[0] = qty[1] = 0; }
  virtual int formulasize() = 0;
  rbool qsimplify() {
    if(!qty[0]) return ftrue;
    if(!qty[1]) return ffalse;
    rbool ret = qsimplify2();
    return ret;
    }
  virtual rbool qsimplify2() = 0;
#endif
  };

struct FormulaQ : Formula {
  bool type;
  rbool right;
  varlist var;
  void setfv() { } // vl_split(fv, right->fv, var);
  FormulaQ(bool t, rbool r) : type(t), right(r) { setfv(); }
  FormulaQ(bool t, rbool r, const varlist& va) : type(t), right(r), var(va) { setfv(); }
  FormulaQ(bool t, rbool r, std::initializer_list<vptr> va) : type(t), right(r) {
    for(auto v: va) var.push_back(v); setfv();
    }
  std::ostream& display (std::ostream &os) const;
  virtual std::ostream& display (std::ostream &os, int prio) { 
    if(prio) { os << "("; display(os); return os << ")"; }
    return display(os);
    }
  rbool subst(const varsubstlist& l) const { 
    auto q = std::make_shared<FormulaQ> (type, right->subst(l), var);
    return q->simplify(q);
    }
  rbool negate() { 
    return rbool(std::make_shared<FormulaQ> (!type, !right, var));
    }
  rbool simplify(std::shared_ptr<FormulaQ> def);
  bool verifyAt(int at);
  bool generateAll(int at, std::shared_ptr<Subdomain>& dom, int val, int rid);
  bool verify() {
    bool b = verifyAt(var.size()-1); qty[b]++; return b; 
    }
  bool uses(vptr v) { return right->uses(v); }
  term valueKnown(vptr v, bool negated);
  rbool qsimplify2() {
    auto q = std::make_shared<FormulaQ> (type, right->qsimplify(), var);
    return q -> simplify(q);
    }
  void clearinfo() { qty[0] = qty[1] = 0; right->clearinfo(); }
  int formulasize() { return 1 + right->formulasize(); };
  void initDomains();
  void listFeatures(featurelist& f) { right->listFeatures(f); }
  };

extern int triesleft, checkid;

// is 'phi' too hard for the internal solver?
bool notForInternalSolver(rbool phi);

}

#endif
