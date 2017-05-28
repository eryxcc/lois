#include "../include/loisinternal.h"

namespace lois {

void subJoin(std::shared_ptr<Subdomain>& s1, std::shared_ptr<Subdomain>& s2);
std::shared_ptr<Subdomain> subJoin(const term& t1, const term &t2);
std::shared_ptr<Subdomain> subGet(vptr v);
std::shared_ptr<Subdomain> subGet(const term& t);

int triesleft, vids, checkid;

#define CountOnce(x) (qty[0]++?0:x)

std::ostream& operator << (std::ostream& os, rbool a) { return a->display(os); }

// a shorthand for negation
rbool operator ! (rbool a) { return a->negate(); }

// a fixed rbool (either true or false)
struct FormulaFixed : Formula {
  bool type;
  FormulaFixed(bool t) : type(t) {}
  std::ostream& display (std::ostream &os) const { 
    // fixed formulae usually are simplified and thus not sent to a solver
    if(outlan == CVC3) return os << (type ? "(1=1)" : "(1=0)");
    return os << (type ? "true" : "false"); 
    }
  virtual rbool subst(const varsubstlist& l) const { return type ? ftrue : ffalse; }
  rbool negate() { return type ? ffalse : ftrue; }
  bool verify() { 
#ifdef AGGSYM
    qty[type]++; 
#endif
    return type; }
  virtual bool uses(vptr v) { return false; }
  void initDomains() {}
  void listFeatures(featurelist& f) {}
#ifdef AGGSYM
  rbool qsimplify2() { return type ? ftrue : ffalse; }
  int formulasize() { return 1; };
#endif
  };

bool isused(vptr v, rbool phi) { return phi.p->uses(v); }

rbool ftrue(std::make_shared<FormulaFixed> (true));
rbool ffalse(std::make_shared<FormulaFixed> (false));

// a Boolean conjunctive, type is true for "AND", false for "OR"

struct FormulaBi : Formula {
  bool type;
  rbool left, right;
  vptr lastusecheck; int usestat;  
  FormulaBi(bool t, rbool l, rbool r) : type(t), left(l), right(r) {
    // vl_merge(fv, left->fv, right->fv);
    }
  std::ostream& display (std::ostream &os) const { 
    if(outlan == CVC3)
      return os << "(" << left << (type ? " AND " : " OR ") << right << ")"; 
    if(outlan == SMT || outlan == SMT_INC)
      return os << "(" << (type ? "and " : "or ") << iindent << left << ieol << right << ieol << ")" << iunindent; 
    if(outlan == SPASS)
      return os << (type ? "and " : "or ") << "(" << iindent << left << "," << ieol << right << ieol << ")" << iunindent; 

    left->display(os, type ? 1:2);
    os << (type ? sym._and : sym._or);
    return right->display(os, type ? 1:2);
    }
  virtual std::ostream& display (std::ostream &os, int prio) { 
    if(prio != (type?1:2)) { os << "("; display(os); return os << ")"; }
    return display(os); 
    }

  virtual rbool subst(const varsubstlist& l) const { 
    if(l.size() == 1 && l[0].first == lastusecheck && usestat < 3) {
      if(usestat == 1 && type)
        return substitute(left, l) && right;
      if(usestat == 1 && !type)
        return substitute(left, l) || right;
      if(usestat == 2 && type)
        return left && substitute(right, l);
      if(usestat == 2 && !type)
        return left || substitute(right, l);
      }
    if(type) return substitute(left, l) && substitute(right, l);
    return substitute(left, l) || substitute(right, l);
    // make_shared<FormulaBi> (type, left->alph(v1,v2), right->alph(v1,v2));
    }
  rbool negate() {
    return rbool(std::make_shared<FormulaBi> (!type, !left, !right));
    }
  bool verify() { 
    bool b;
    if(type) b= left->verify() && right->verify(); 
    else b= left->verify() || right->verify();
#ifdef AGGSYM
    qty[b]++; 
#endif
    return b;
    }
  bool uses(vptr v) { 
    if(v == lastusecheck) return usestat;
    lastusecheck = v; usestat = 0;
    if(left->uses(v)) usestat |= 1;
    if(right->uses(v)) usestat |= 2;
    return usestat; }
  term valueKnown(vptr v, bool negated) { 
    // non-negated: AND => at least one
    if(negated == type) {
      // OR: both need the same value
      term v1 = left->valueKnown(v, negated);
      term v2 = right->valueKnown(v, negated);
      if(v1.p == v2.p) return v1;
      return nullterm;
      }
    else {
      term v1 = left->valueKnown(v, negated);
      if(v1.p) return v1;
      return right->valueKnown(v, negated);
      }
    }
#ifdef AGGSYM
  rbool qsimplify2() { 
    if(type) return left->qsimplify() && right->qsimplify();
    else return left->qsimplify() || right->qsimplify();
    }
  void clearinfo() { qty[0] = qty[1] = 0; left->clearinfo(); right->clearinfo(); }
  int formulasize() { return left->formulasize() + 1 + right->formulasize(); };
#endif
  void initDomains() { left->initDomains(); right->initDomains(); }
  void listFeatures(featurelist& f) { left->listFeatures(f); right->listFeatures(f); }
  };

// an equality (type=true) or inequality (type=false)

struct FormulaEq : Formula {
  bool type;
  term t1, t2;
  FormulaEq(bool t, term v1, term v2) : type(t), t1(v1), t2(v2) {
    // todo
    // if(t1 == t2) aerror("ta sama zmienna", t1, t2);
    /* if(v1->sortid < v2->sortid) {
      fv.push_back(v1);
      fv.push_back(v2);
      }
    else {
      fv.push_back(v2);
      fv.push_back(v1);
      } */
    }
  std::ostream& display (std::ostream &os) const { 
    if((outlan == SMT || outlan == SMT_INC) && type)
      return os << "(= " << t1 << " " << t2 << ")"; 
    if((outlan == SMT || outlan == SMT_INC) && !type)
      return os << "(not (= " << t1 << " " << t2 << "))";

    if(outlan == SPASS && type)
      return os << "equal(" << t1 << ", " << t2 << ")"; 
    if(outlan == SPASS && !type)
      return os << "not(equal(" << t1 << ", " << t2 << "))";

    return os << t1 << (type ? sym.eq : sym.neq) << t2;
    }
  virtual rbool subst(const varsubstlist& l) const { 
    if(type)
      return substitute(t1, l) == substitute(t2, l);
    return substitute(t1, l) != substitute(t2, l);
    }
  rbool negate() { return rbool(std::make_shared<FormulaEq> (!type, t1, t2)); }
  bool verify() { 
    bool b;
    if(type) b= t1.p->getValue() == t2.p->getValue(); 
    else b = t1.p->getValue() != t2.p->getValue();
#ifdef AGGSYM
    qty[b]++; 
#endif
    return b;
    }
  bool uses(vptr v) { return t1.p->uses(v) || t2.p->uses(v); }
  term valueKnown(vptr v, bool negated) { 
    // non-negated: AND => at least one
    if(negated == !type) {
      if(t1.asVar() == v && t2.asVar()) return term(t2.asVar());
      if(t2.asVar() == v && t1.asVar()) return term(t1.asVar());
      }
    return nullterm;
    }
#ifdef AGGSYM
  rbool qsimplify2() {
    return rbool(std::make_shared<FormulaEq> (type, t1, t2)); 
    }
  int formulasize() { return 1; };
#endif
  void initDomains() { subJoin(t1, t2)->numEq++; }
  void listFeatures(featurelist& f) { t1.p->listFeatures(f); t2.p->listFeatures(f); }
  };

struct FormulaBinary : Formula {
  int id;
  Relation *rel;
  term t1, t2;
  FormulaBinary(Relation *r, int i, term v1, term v2) : rel(r), id(i), t1(v1), t2(v2) {
    /* if(v1->sortid < v2->sortid) {
      fv.push_back(v1);
      fv.push_back(v2);
      }
    else {
      fv.push_back(v2);
      fv.push_back(v1);
      } */
    }
  std::ostream& display (std::ostream &os) const { 
    rel->display(os, id, t1, t2);
    return os;
    }
  rbool subst(const varsubstlist& l) const { 
    return rel->binform(id, substitute(t1,l), substitute(t2,l));
    }
  rbool negate() { return rel->binform(rel->negate(id), t1, t2); }
  bool verify() { 
    std::shared_ptr<Subdomain> d = subGet(t1);
    if(d != subGet(t2)) throw subdomain_exception();
    int i = 0;
    while(d->rels[i]->rel != rel) i++;
    bool b = d->rels[i]->check(id, t1.p->getValue(), t2.p->getValue());
#ifdef AGGSYM
    qty[b]++; 
#endif
    return b;
    }
  bool uses(vptr v) { return t1.p->uses(v) || t2.p->uses(v); }
  term valueKnown(vptr v, bool negated) { 
    return nullterm;
    }
#ifdef AGGSYM
  rbool qsimplify2() {
    return rel->binform(id, t1, t2);
    }
  int formulasize() { return 1; };
#endif
  void initDomains() { 
    subJoin(t1, t2)->useRelation(rel); 
    }
  void listFeatures(featurelist& f) { 
    t1.p->listFeatures(f); t2.p->listFeatures(f); 
    f.relations.insert(rel);
    }
  };

std::ostream& FormulaQ::display (std::ostream &os) const { 
  if(!var.size()) return os << right;
  if(outlan == CVC3) {
    os << (type ? "(FORALL(" : "(EXISTS(");
    bool comma = false;
    for(auto v: var) {
      if(comma) os << ", ";
      comma = true;
      os << v << ":" << v->dom->name;
      }
    return os << "): " << right << ")";
    }
  if(outlan == SMT) {
    os << (type ? "(forall " : "(exists ");
    for(auto v: var) {
      os << "(" << v << " " << v->dom->name << ") ";
      }
    return os << iindent << right << ieol << ")" << iunindent;
    }
  if(outlan == SMT_INC) {
    os << (type ? "(forall (" : "(exists (");
    for(auto v: var) {
      os << "(" << v << " " << smtsort() << ") ";
      }
    return os << ")" << iindent << right << ieol << ")" << iunindent;
    }
  if(outlan == SPASS) {
    os << (type ? "forall([ " : "exists([ ");
    bool had = false;
    for(auto v: var) {
      if(had) os << ", ";
      os << v->dom->name << "(" << v << ")"; had = true;
      }
    return os << "], " << iindent << right << ieol << ")" << iunindent;
    }
  os << (type ? sym.forall : sym.exists);
  for(auto v: var) {
    os << v << sym.in << v->dom->name;
    /* os << "[" << subFindv(v);
    auto dom = subGet(v);
    if(dom->numOrder) os << "!";
    if(dom->numEdge) os << "?";
    // os << dom->numEq << ":" << dom->numNeq << ":" << dom->numQuantifier << ":" << dom->numOrder;
    os << "]"; */
    os << " ";
    }
  return os << right;
  }

void decomposeBin(rbool b, bool bt, std::vector<rbool>& v) {
  auto dec = std::dynamic_pointer_cast<FormulaBi> (b.p);
  if(dec && dec->type == bt) {
    decomposeBin(dec->left, bt, v);
    decomposeBin(dec->right, bt, v);
    }
  else
    v.push_back(b);
  }

rbool FormulaQ::simplify(std::shared_ptr<FormulaQ> def) {
  
  using namespace std;
  // we know the value of something, maybe?
  for(size_t i=0; i<var.size(); i++) {
    vptr v = var[i];
    if(!right->uses(v)) {
      var[i] = var[var.size()-1];
      var.pop_back(); i--;
      }
    else {
      term v1 = right->valueKnown(v, type);
      if(v1.p) {
        right = substitute(right, v, v1);
        var[i] = var[var.size()-1];
        var.pop_back(); i--;
        }
      }
    }
  
  // no variables
  if(!var.size()) return right;    
  // quantified true/false
  auto fix = std::dynamic_pointer_cast<FormulaFixed> (right.p);
  if(fix) return right;
  
  /*
  // quantified binary, one side does not use the vars
  auto bin = std::dynamic_pointer_cast<FormulaBi> (right.p);
  if(bin) {
    bool leftuses = false, rightuses = false;
    for(auto v: var) leftuses |= bin->left->uses(v), rightuses |= bin->right->uses(v);

    if(!leftuses && bin->type == true) return bin->left && makeq(type, bin->right, var);
    if(!rightuses && bin->type == true) return bin->right && makeq(type, bin->left, var);
    if(!leftuses && bin->type == false) return bin->left || makeq(type, bin->right, var);
    if(!rightuses && bin->type == false) return bin->right || makeq(type, bin->left, var);
    } */
  
  // quantified (in)equality
  // as long as we assume that !right->unused(all vars),
  // and that they are variables,
  // we know the value
  auto eq = std::dynamic_pointer_cast<FormulaEq> (right.p);
  if(eq && eq->t1.asVar() && eq->t2.asVar()) return type ? ffalse : ftrue;
  
  // take the independent components out of the quantifier
  auto bin = std::dynamic_pointer_cast<FormulaBi> (right.p);
  if(bin) {
    bool btype = bin->type;
    std::vector<rbool> phi;
    std::vector<int> vused;
    decomposeBin(right, btype, phi);
    bool suc = false;
    for(size_t k=0; k<phi.size(); k++) {
      size_t li = 0;
      for(size_t i=0; i<var.size(); i++)
        if(phi[k]->uses(var[i])) li = i+1;
      if(li < var.size()) suc = true;
      vused.push_back(li);
      }
    
    if(suc) {
      for(int k=1; k<phi.size(); ) {
        if(k && vused[k] < vused[k-1]) 
          std::swap(vused[k], vused[k-1]), std::swap(phi[k], phi[k-1]), k--;
        else
          k++;
        }
      
      int at = phi.size()-1;
      rbool nphi = phi[at];
      int lev = vused[at];
  
/*    std::cout << "SIMPLIFYING: ";
      display(std::cout);
      std::cout << std::endl;
      
      for(int k=0; k<phi.size(); k++)
        std::cout << k << "." << vused[k] << ") " << phi[k] << std::endl; */
      
      while(true) {
        at--;
        int nlev = at >= 0 ? vused[at] : 0;
        if(lev > nlev) {
          varlist v;
          while(lev > nlev) { lev--; v.push_back(var[lev]); }
          nphi = makeq(type, nphi, v);
          }
        if(at < 0) {
/*        std::cout << "SIMPLIFIED: ";
          display(std::cout);
          std::cout << " TO: " << nphi << std::endl; */
          return nphi;
          }
        else if(btype) nphi = phi[at] && nphi;
        else nphi = phi[at] || nphi;
        }
      }
    }
  
  return rbool(def);
  }

bool FormulaQ::verifyAt(int at) {
  if(at == -1) return right->verify();
  auto dom = subGet(var[at]);
  int& nextval(dom->nextval);
  triesleft -= nextval;
  if(triesleft < 0) return false;
  for(int z=0; z<nextval; z++) { 
    var[at]->value = z; 
    if(verifyAt(at-1) != type) return !type;
    }
  var[at]->value = nextval;
  nextval++;
  bool b = generateAll(at, dom, nextval-1, 0);
  nextval--;
  return b;
  }

bool FormulaQ::generateAll(int at, std::shared_ptr<Subdomain>& dom, int val, int rid) {
  if(rid >= int(dom->rels.size())) {
    if(val < dom->nextval-1) 
      return generateAll(at, dom, val+1, 0);
    else
      return verifyAt(at-1);
    }
  else return dom->rels[rid]->exists(this, at, dom, val, rid);
  }

term FormulaQ::valueKnown(vptr v, bool negated) {
  term v1 = right->valueKnown(v, negated);
  if(!v1.p) return v1;
  for(auto qx: var) if(qx == v1.asVar()) {
    // should not happen if simplified!
    std::cout << "should not happen: v= " <<v << " rbool= ";
    display(std::cout); std::cout << std::endl;
    term vk = right->valueKnown(qx, type);
    if(!vk.p) std::cout << "value unknown!" << std::endl;
    else std::cout << "valueKnown is " << vk <<std::endl;
    return nullterm;
    }
  return v1;
  }

void FormulaQ::initDomains() { 
  for(auto q: var) subGet(q)->numQuantifier++;
  right->initDomains();
  }


rbool makeq(bool t, rbool r, const varlist& va) {
  auto q = std::make_shared<FormulaQ>  (t,r,va);
  return q->simplify(q);
  }

// AND and OR, more effective than using FormulaBi due to simplifications
rbool operator && (rbool a, rbool b) {
  if(a.isTrue()) return b;
  if(b.isTrue()) return a;
  if(a.isFalse() || b.isFalse()) return ffalse;
  if(a.p == b.p) return a;

  auto ae = std::dynamic_pointer_cast<FormulaEq>(a.p);
  auto be = std::dynamic_pointer_cast<FormulaEq>(b.p);
  if(ae) if(be) {
    if(eqterm(ae->t1, be->t1) && eqterm(ae->t2, be->t2))
      return ae->type == be->type ? a : ffalse;
    if(eqterm(ae->t1, be->t2) && eqterm(ae->t2, be->t1))
      return ae->type == be->type ? a : ffalse;
    }

  // if(a->formulasize() > b->formulasize()) swap(a, b);
  
  return rbool(std::make_shared<FormulaBi> (true, a, b));
  }

rbool operator || (rbool a, rbool b) {
  if(a.isFalse()) return b;
  if(b.isFalse()) return a;
  if(a.isTrue() || b.isTrue()) return ftrue;
  if(a.p == b.p) return a;

  auto ae = std::dynamic_pointer_cast<FormulaEq>(a.p);
  auto be = std::dynamic_pointer_cast<FormulaEq>(b.p);
  if(ae) if(be) {
    if(eqterm(ae->t1, be->t1) && eqterm(ae->t2, be->t2)) {
      return ae->type == be->type ? a : ftrue;
      }
    if(eqterm(ae->t1, be->t2) && eqterm(ae->t2, be->t1)) {
      return ae->type == be->type ? a : ftrue;
      }
    }

  return rbool(std::make_shared<FormulaBi> (false, a, b));
  }

rbool operator == (const term& a, const term& b) { 
  if(a.p->getDom() != b.p->getDom()) return ffalse;
  if(eqterm(a, b)) return ftrue;
  else return rbool(std::make_shared<FormulaEq> (true, a, b)); 
  }

rbool operator != (const term& a, const term& b) { 
  if(a.p->getDom() != b.p->getDom()) return ftrue;
  if(eqterm(a, b)) return ffalse;
  return rbool(std::make_shared<FormulaEq> (false, a, b)); 
  }

rbool makebinary(Relation *rel, int id, const term& a, const term& b) {
  return rbool(std::make_shared<FormulaBinary> (rel,id,a,b));
  }

void Subdomain::useRelation(Relation *r) {
  for(size_t i=0; i<rels.size(); i++) 
    if(rels[i]->rel == r) return;
  rels.push_back(r->genSub());
  }

void autoname(std::ostream& os, int id) {
  if(id >= 26) autoname(os, id/26-1); 
  os << char('a' + (id%26));
  }

int sortids = 0;
std::string autoprefix = "";

std::ostream& Variable::display (std::ostream& os) const {
  if(!id) id = ++vids;
  if(outlan == SMT || outlan == SMT_INC) os << "?";
  if(outlan == SPASS) os << "lois";
  if(name != "") os << name;
  else { os << autoprefix; autoname(os, id-1); }
  return os;
  }

std::shared_ptr<Subdomain> nullsub;

// Find & Union
void subInit(vptr v) {
  if(v->curcheck < checkid) {
    v->curcheck = checkid;
    v->sublink = v;
    v->sub = nullsub;
    }
  }

int recdepth;

vptr Variable::subFindv(tvptr v) {
  vptr vv = std::dynamic_pointer_cast<Variable> (v);
  subInit(vv);
  vv->sublink = (vv->sublink == vv) ? vv : subFindv(vv->sublink);
  return vv->sublink;
  }

void subJoin(std::shared_ptr<Subdomain>& s1, std::shared_ptr<Subdomain>& s2);

vptr TermBinary::subFindv(tvptr v) {
  if(v->curcheck < checkid) {
    vptr v1 = left.p->subFindv(left.p);
    vptr v2 = right.p->subFindv(right.p);
    if(v1 == nullvptr) return v2;
    if(v2 == nullvptr) return v1;
    if(v1->sublink != v2->sublink) {
      v2->sublink = v1->sublink;
      subJoin(v1->sub, v2->sub);
      }
    if(!v1->sub) v1->sub = std::make_shared<Subdomain> ();
    v1->sub->useRelation(r);
    sublink = v1->sublink;
    }
  return sublink;
  }

vptr TermConst::subFindv(tvptr inv) {
  vptr vv = std::dynamic_pointer_cast<Variable> (inv);
  if(!vv) vv = v;
  subInit(vv);
  vv->sublink = (vv->sublink == vv) ? vv : subFindv(vv->sublink);
  return vv->sublink;
  }

void TermBinary::listFeatures(featurelist& f) {
  f.relations.insert(r);
  left.p->listFeatures(f);
  right.p->listFeatures(f);
  }

void TermConst::listFeatures(featurelist& f) {
  f.relations.insert(r);
  }

void Variable::listFeatures(featurelist& f) {
  f.domains.insert(dom);
  }

std::shared_ptr<Subdomain> subGet(vptr v) {
  v = v->subFindv(v);
  if(!v->sub) v->sub = std::make_shared<Subdomain> ();
  return v->sub;
  }

std::shared_ptr<Subdomain> subGet(const term& t) {
  vptr v = t.p->subFindv(t.p);
  if(v == nullvptr) return nullsub;
  if(!v->sub) v->sub = std::make_shared<Subdomain> ();
  return v->sub;
  }

void subJoin(std::shared_ptr<Subdomain>& s1, std::shared_ptr<Subdomain>& s2) {
  if(s2 == nullsub) return;
  if(s1 == nullsub) { swap(s1, s2); return; }
  s1->numEq += s2->numEq;
  s1->numQuantifier += s2->numQuantifier;
  for(size_t i=0; i<s2->rels.size(); i++) 
    s1->useRelation(s2->rels[i]->rel);
  s2 = nullsub;
  }

std::shared_ptr<Subdomain> subJoin(const term& t1, const term &t2) {
  vptr v1 = t1.p->subFindv(t1.p);
  vptr v2 = t2.p->subFindv(t2.p);
  if(v1->sublink != v2->sublink) {
    v2->sublink = v1->sublink;
    subJoin(v1->sub, v2->sub);
    }
  if(!v1->sub) v1->sub = std::make_shared<Subdomain> ();
  return v1->sub;
  }

int TermBinary::getValue() { 
  std::shared_ptr<Subdomain> d = subGet(left);
  if(d != subGet(right)) throw subdomain_exception();
  int i = 0;
  while(d->rels[i]->rel != r) i++;
  int r = d->rels[i]->checkterm(op, left.p->getValue(), right.p->getValue());
  return r;
  }

int TermConst::getValue() {
  return r->checktermconst(op);
  }

vptr term::asVar() const { return std::dynamic_pointer_cast<Variable> (p); }

term::term(tvptr v) : p(v) {}

std::ostream& operator << (std::ostream& os, vptr a) { return a->display(os); }

term substitute(vptr vthis, const varsubstlist& l) {
  for(auto& p: l) if(vthis == p.first) return p.second;
  return term(vthis);
  }

// alpha-convert the term 'tthis' from 'v1' to 'v2'
term substitute(const term& tthis, const varsubstlist& l) {
  return tthis.p->subst(tthis, l);
  }

bool isused(vptr vthis, vptr v) {
  return vthis == v;
  }

/*
void vl_merge(varlist& vres, const varlist& v1, const varlist& v2) {
  vres.clear();
  size_t i = 0;
  size_t j = 0;
  while(i < v1.size() && j<v2.size())
    if(v1[i]->sortid < v2[j]->sortid) vres.push_back(v1[i++]); else 
    if(v1[i]->sortid > v2[j]->sortid) vres.push_back(v2[j++]); else 
    vres.push_back(v1[i++]), j++;
  while(i < v1.size()) vres.push_back(v1[i++]);
  while(j < v2.size()) vres.push_back(v2[j++]);
  }

// remove all variables in v2 from v1
void vl_split(varlist& vres, const varlist& v1, const varlist& v2) {
  vres.clear();
  size_t i = 0;
  size_t j = 0;
  while(i < v1.size() && j<v2.size())
    if(v1[i]->sortid < v2[j]->sortid) vres.push_back(v1[i++]); else 
    if(v1[i] == v2[j]) i++, j++; else 
    j++;
  while(i < v1.size()) vres.push_back(v1[i++]);
  }
*/

rbool makequantifier(bool f, rbool& phi, varlist to_quantify) {
  auto sphi = std::make_shared<FormulaQ> (f, phi, to_quantify);
  return sphi->simplify(sphi);
  }

rbool substitute(const rbool& phi, const varsubstlist& l) {
  return phi->subst(l);
  }

term valueKnown(rbool& phi, vptr v, bool negated) {
  return phi->valueKnown(v, negated);
  }

const term nullterm;

}


