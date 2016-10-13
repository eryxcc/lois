#include "../include/loisinternal.h"

namespace lois {

// ORDER RELATION

#define ID_GREATER 1
#define ID_LEQ     2
#define ID_EQUIV   3

#define ID_BINARY 1
#define ID_NOBINARY 2

struct SubOrder : SubRelation {
  std::vector<int> position;

  SubOrder(Relation *r) : SubRelation(r) {}
  
  bool exists(FormulaQ* q, int at, std::shared_ptr<Subdomain>& dom, int val, int rid) {
    if(int(position.size()) < val+1) position.resize(val+1);
    for(int i=0; i<=val; i++) {
      position[val] = i;
      if(q->generateAll(at, dom, val, rid+1) != q->type) return !q->type;
      }
    return q->type;
    }

  bool check(int id, int val1, int val2) {
    if(val1 == val2) return id == ID_LEQ;
    else {
      int v1 = position[val1], v2 = position[val2];
      for(int k=val1+1; k<=val2; k++) if(position[k] <= v1) v1++;
      for(int k=val2+1; k<=val1; k++) if(position[k] <= v2) v2++;
      if(id == ID_LEQ) return v1 <= v2;
      else return v1 > v2;
      }
    }

  int checkterm(int id, int val1, int val2) { 
    if(check(ID_LEQ, val1, val2) ? (id==ID_MIN) : (id==ID_MAX)) 
      return val1;
    else return val2;
    }
  };

subrelptr RelOrder::genSub() {
  return std::make_shared<SubOrder>(this);
  }

rbool RelOrder::binform(int id, const term& a, const term& b) {
  if(a.p->getDom() != b.p->getDom()) return ffalse;
  if(eqterm(a, b)) return id == ID_LEQ ? ftrue : ffalse;
  return makebinary(this, id, a, b);
  };


// === reals ===

struct SubOrderedField : SubRelation {

  SubOrderedField(Relation *r) : SubRelation(r) {}
  
  bool exists(FormulaQ* q, int at, std::shared_ptr<Subdomain>& dom, int val, int rid) {
    throw unsolvable_exception();
    }

  bool check(int id, int val1, int val2) {
    throw unsolvable_exception();
    }

  int checkterm(int id, int val1, int val2) { 
    throw unsolvable_exception();
    }
  };

subrelptr RelOrderedField::genSub() {
  return std::make_shared<SubOrderedField>(this);
  }

std::vector<int> constvalues;

std::ostream& RelOrderedField::displaytermbin(std::ostream& os, const TermBinary* t) {
    int id = t->op;
    if(outlan == SMT || outlan == SMT_INC) {
      if(id == ID_MIN) return os << "(min " << " " << t->left << " " << t->right << ")";
      if(id == ID_MAX) return os << "(max " << " " << t->left << " " << t->right << ")";
      if(id == ID_PLUS) return os << "(+ " << t->left << " " << t->right << ")";
      if(id == ID_TIMES) return os << "(* " << t->left << " " << t->right << ")";
      if(id == ID_MINUS) return os << "(- " << t->left << " " << t->right << ")";
      if(id == ID_DIVIDE) return os << "(/ " << t->left << " " << t->right << ")";
      }
    os << "(" << t->left;
    if(id == ID_MIN) os << opmin;
    if(id == ID_MAX) os << opmax;
    if(id == ID_PLUS) os << opplus;
    if(id == ID_MINUS) os << opminus;
    if(id == ID_TIMES) os << optimes;
    if(id == ID_DIVIDE) os << opdivide;
    return os << t->right << ")";
    }

std::ostream& RelOrderedField::displayconst(std::ostream& os, const TermConst* t) {
    int id = t->op;
    return os << constvalues[id];
    }

term RelOrderedField::constant(Domain *d, int i) {
  int k = constvalues.size();
  constvalues.push_back(i);
  return term(std::make_shared<TermConst> (this,k,d));
  }

// === binary relation ===

struct SubBinary : SubRelation {
  std::vector<int> edges;
  loopmode lm;
  symmode sm;
  SubBinary(Relation *r, loopmode l, symmode s) : SubRelation(r), lm(l), sm(s) {}
  
  bool exists(FormulaQ* q, int at, std::shared_ptr<Subdomain>& dom, int val, int rid) {
    if(int(edges.size()) < val+1) edges.resize(val+1);
    
    int tochoose = val;
    
    if(lm == lmPossibleLoops) tochoose ++;
    if(sm == smAsymmetric) tochoose += val;
    
    for(int i=0; i<(1<<tochoose); i++) {
      edges[val] = i;
      if(q->generateAll(at, dom, val, rid+1) != q->type) return !q->type;
      }
    return q->type;
    }

  bool check(int id, int val1, int val2) {
    // printf("val1=%d val2=%d\n", val1, val2);
    if(val1 == val2) {
      if(lm == lmNoLoops) return id == ID_NOBINARY;
      if(lm == lmAllLoops) return id = ID_BINARY;
      }
    if(val1 > val2) {
      std::swap(val1, val2);
      if(sm == smAsymmetric) { val1 += val2; if(lm == lmPossibleLoops) val1++; }
      if(sm == smAntisymmetric) id ^= (ID_NOBINARY ^ ID_BINARY);
      }
      
    // printf("-> val1=%d val2=%d\n", val1, val2);
    bool b = edges[val2] & (1<<val1);
    if(id == ID_NOBINARY) return !b;
    return b;
    }
  };

subrelptr RelBinary::genSub() {
  return std::make_shared<SubBinary> (this, lm, sm);
  }

rbool RelBinary::binform(int id, const term& a, const term& b) {
  if(a.p->getDom() != b.p->getDom()) return ffalse;
  if(eqterm(a,b) && lm != lmPossibleLoops) return ((id == ID_BINARY) ^ (lm == lmNoLoops)) ? ftrue : ffalse;
  return makebinary(this,id,a,b);
  };

// === partition of the universe in 'nogroups' groups ===

#define ID_UN_EQ -1
#define ID_UN_NEQ -2

struct SubUnary : SubRelation {
  std::vector<int> groups;
  int nogroups;
  SubUnary(Relation *r, int n) : SubRelation(r), nogroups(n) {}
  
  bool exists(FormulaQ* q, int at, std::shared_ptr<Subdomain>& dom, int val, int rid) {
    if(groups.size() < val+1) groups.resize(val+1);
    
    for(int i=0; i<nogroups; i++) {
      groups[val] = i;
      if(q->generateAll(at, dom, val, rid+1) != q->type) return !q->type;
      }
    return q->type;
    }

bool check(int id, int val1, int val2) {
  if(id == ID_UN_EQ) return groups[val1] == groups[val2];
  if(id == ID_UN_NEQ) return groups[val1] != groups[val2];
  // printf("val1=%d val2=%d\n", val1, val2);
  if(val1 != val2) throw loisexception();
  return ((1<<groups[val1]) & id);
  }
};

subrelptr RelUnary::genSub() {
  return std::make_shared<SubUnary> (this, nogroups);
  }

std::ostream& RelUnary::display(std::ostream& os, int id, const term& t1, const term& t2) {
  os << oprel;

  if(id == ID_UN_EQ)
    return os << "[" << t1 << sym.eq << t2 << "]";

  if(id == ID_UN_NEQ)
    return os << "[" << t1 << sym.neq << t2 << "]";

  for(int i=0; i<nogroups; i++) if(id & (1<<i)) os << i;
  return os << "(" << t1 << ")";
  }

rbool RelUnary::binform(int id, const term& a, const term& b) {
  if(id >= 0) {
    int mask = ((1<<nogroups)-1);
    // note: no domain check -> they can actually be in different domains
    if((id & mask) == 0) return ffalse;
    if((id & mask) == mask) return ftrue;
    }
  if(id == ID_UN_EQ && eqterm(a,b)) return ftrue;
  if(id == ID_UN_NEQ && eqterm(a,b)) return ffalse;
  return makebinary(this,id,a,b);
  };

int RelUnary::negate(int id) { 
  if(id == ID_UN_EQ) return ID_UN_NEQ;
  if(id == ID_UN_NEQ) return ID_UN_EQ;
  return (~id) & 0x7FFFFFFF;
  }

// === infinite random tree ===

#define ID_ANCESTOR_OR_EQUAL 1
#define ID_NOT_ANCESTOR_OR_EQUAL 2

struct SubRelTree : SubRelation {
  std::vector<int> parent;
  SubRelTree(Relation *r) : SubRelation(r) {}
  
  bool exists(FormulaQ* q, int at, std::shared_ptr<Subdomain>& dom, int val, int rid) {

    if(parent.size() < val+2) parent.resize(val+2);
  
    // create a new tree
    if(val == 0) {
      parent[0] = -1;
      if(q->generateAll(at, dom, val, rid+1) != q->type) return !q->type;
      }
    
    // already created
    if(val > 0 && parent[val-1] == val) 
      return q->generateAll(at, dom, val, rid+1);
    
    // create a descendant
    for(int u=0; u<val; u++) {
      parent[val] = u;
      if(q->generateAll(at, dom, val, rid+1) != q->type) return !q->type;
      }
    
    // create an intermediate ancestor
    for(int u=0; u<val; u++) {
      parent[val] = parent[u]; parent[u] = val;
      bool b = q->generateAll(at, dom, val, rid+1);
      parent[u] = parent[val];
      if(b != q->type) return !q->type;
      }
    
    // create another child of an intermediate ancestor of u
    for(int u=0; u<val; u++) {
      dom->nextval++;
      parent[val+1] = parent[u]; parent[u] = val+1; parent[val] = val+1;
      bool b = q->generateAll(at, dom, val, rid+1);
      parent[u] = parent[val+1]; dom->nextval--;
      if(b != q->type) return !q->type;
      }
    
    return q->type;
    }

  bool check(int id, int val1, int val2) {
    while(true) {
      if(val2 == -1) return id == ID_NOT_ANCESTOR_OR_EQUAL;
      if(val2 == val1) return id == ID_ANCESTOR_OR_EQUAL;
      val2 = parent[val2];
      }
    }

  int checkterm(int id, int val1, int val2) { 
    int l1 = 0; int up1 = val1;
    while(up1 != -1) up1 = parent[up1], l1++;

    int l2 = 0; int up2 = val2;
    while(up2 != -1) up2 = parent[up2], l2++;
    
    up1 = val1; up2 = val2;
    while(l2 > l1) l2--, up2 = parent[up2];
    while(l1 > l2) l1--, up1 = parent[up1];
    
    while(up1 != up2) up1 = parent[up1], up2 = parent[up2];

    return up1;
    }
  };

rbool RelTree::binform(int id, const term& a, const term& b) {
  if(a.p->getDom() != b.p->getDom()) return ffalse;
  if(eqterm(a,b)) return id == ID_ANCESTOR_OR_EQUAL ? ftrue : ffalse;
  return makebinary(this,id,a,b);
  };


subrelptr RelTree::genSub() {
  return std::make_shared<SubRelTree> (this);
  }

}
