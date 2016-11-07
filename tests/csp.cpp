// Locally Finite Constraint Satisfaction Problems
// Bartek Klin, Eryk Kopczyski, Joanna Ochremiak, Szymon Toruczyk 
  
#include <stdlib.h>
#include <sys/time.h>

#include <set>

#include <map>
// #include <gecode/driver.hh>
// #include <gecode/int.hh>
// #include <gecode/minimodel.hh>
//
// using namespace Gecode;


#include "../include/loisextra.h"

#include <iostream>
using namespace std;
using namespace lois;



enum symmetry {EQ, LIN};

long long getVa() {
  struct timeval tval;
  gettimeofday(&tval, NULL);
  return tval.tv_sec * 1000000 + tval.tv_usec;
  }

long long lasttime;

void showTimeElapsed() {
  long long curtime = getVa();
  
  if(curtime - lasttime < 100000) 
    cout << "time elapsed: " << (curtime-lasttime) << " \u03bcs" << endl;

  else
    printf("time elapsed: %.6f s\n", (curtime-lasttime)/1000000.0);

  lasttime = curtime;
  }



void listAllChoices(lset& orbitids, std::vector<vptr>& variables, elem i, int a, int b, int& orbitcount, const symmetry &sym=LIN) {
  b++;
  if(a == b) { a++; b=0; }
  if(a == (int) variables.size()) {
    
    lbool neworbit = true;
    for(auto el: orbitids) 
      If(neworbit && (as<elpair>(el).first == i)) neworbit = false;

    If(neworbit) {
        orbitcount++;
        orbitids += elpair(i, orbitcount);
    }
    
    return;
    }
  
  term va(variables[a]);
  term vb(variables[b]);
  
  if (sym==LIN) {  
    If(va < vb) 
      listAllChoices(orbitids, variables, i, a, b, orbitcount, sym);

    If(va > vb)
      listAllChoices(orbitids, variables, i, a, b, orbitcount, sym);

    If(va == vb)
      listAllChoices(orbitids, variables, i, a, b, orbitcount, sym);
  }
  else if (sym==EQ) {
    If(va != vb)
      listAllChoices(orbitids, variables, i, a, b, orbitcount, sym);

    If(va == vb)
      listAllChoices(orbitids, variables, i, a, b, orbitcount, sym);
    }
  else throw exception();
}
  

/*  void addConstraint(lset& orbitids, const char* relname, const eltuple& inrel, 
    int inrel_pos = 0, vector<int> orbitsofar = vector<int> ()) {
    if(inrel_pos == inrel.size()) {
      std::cout << relname;
      for(int i=0; i<inrel_pos; i++) std::cout << " x" << orbitsofar[i];
      std::cout << std::endl;
      return;
      }
  
    for(auto o: orbitids) {
      elpair p = as<elpair> (o);
      If(p.first == inrel[inrel_pos]) {
        auto orbitsofar2 = orbitsofar;
        orbitsofar2.push_back(as<int> (p.second));
        addConstraint(orbitids, relname, inrel, inrel_pos+1, orbitsofar2);
        }
      }
    }
  */
  

lset computeOrbits(const lset& V, const symmetry& sym=LIN) {
  lset orbitids = newSet();
  int orbitcount = 0;
  for(elem v: V) {
    std::vector<vptr> variables;

    contextptr e = currentcontext;
    while(e != emptycontext) {
      if(!e) throw env_exception();
      for(auto w: e->var) variables.push_back(w);
      e = e->parent;
      }
    

    listAllChoices(orbitids, variables, v, 1, -1, orbitcount,sym);
        
    }
    return orbitids;
  }
  

  template<class T> void powersCollect(const lset& X, T beg, T end, int n, eltuple& v, lset& res)
  //computes a tuple of powers  X^k, for k ranging from beg to end using iterator
  {    
    if(beg==end) return;
    
    if(*beg == n) {
      res += v;
      beg++;
    }
    
    for(elem a: X) {
      v.push_back(a);
      n++;
      powersCollect(X, beg, end, n, v, res);
      v.pop_back();
      n--;
    }
  }    
    
  template <class Container> lset powers(const lset& X, const Container & powers) 
  //computes the union of k-th powers of X, for k in the vector powers
  {
    lset res = newSet();
    eltuple v = eltuple();

    powersCollect(X, powers.begin(), powers.end(), 0, v, res);
    
    return res;    
  }
  



  
//typedef std::pair<elem, eltuple> Constraint;
typedef lpair<elem,eltuple> Constraint;

struct ConstraintGraph 
  //a constraint graph is a set of ``variables'' and a set of ``constraints'',
  //where each constraint is a pair (id, t) where id is and  and t is a tuple of variables.
{
  lset variables; 
  lsetof<Constraint> constraints; 
  
  ConstraintGraph() {
    variables = newSet();
    constraints = newSet();
  }
  
  friend ostream& operator<<(ostream& os, const ConstraintGraph& dt);  
    

  int constraintSize(Constraint c)
  {
    return (as<eltuple>(c.second)).size();    
  }
  
  std::set<int> constraintSizes() {
    std::set<int> s;
    
    for (Constraint c: constraints)
      s.insert(constraintSize(c));
    
    return s;
  }
  
  
  //takes an n-tuple of k-tuples and returns a k-tuple of n-tuples
  eltuple transpose(const eltuple& tup) {
    eltuple res;
        
    for (int i=0; true; i++) {
      eltuple coord;
      for (int j=0; j<tup.size(); j++) {
        if (as<eltuple>(tup[j]).size()==i) {return res;}        
        coord.push_back(as<eltuple>(tup[j])[i]);
      }
      res.push_back(coord);      
    }
    return res;
  }
  
  ConstraintGraph unarizedConstraintGraph()
      //Computes the ConstraintGraph whose variables are n-tuples of variables of I for n∈lengths
    //(for n=1 tuples are unpacked)
      //and whose constraints of length n, for n∈lengths are replaced by unary contraints,
    //and with binary constraints for the graphs of the n projections of n-tuples
    {
      ConstraintGraph R;
      
      lsetof<eltuple> newV = newSet();
      
      set<int> arities = constraintSizes();
      arities.erase(1);
      arities.insert(3);
      newV = powers(variables, arities);
      
      for (Constraint c: constraints) {
        if (constraintSize(c)==1)
          R.constraints += c;
        else 
          R.constraints += Constraint(c.first, eltuple({c.second}));
      }       
       
      for (auto t: newV) {
        int n = as<eltuple>(t).size();
        for (int k=0; k<n ; k++)
          R.constraints += Constraint(elpair(n,k), eltuple({t, t[k]}));
      }
      
      R.variables =  variables | newV;
 
      return R;
    }
  
    ConstraintGraph( const ConstraintGraph &obj) //copy constructor
    {
      variables = obj.variables;
      constraints = obj.constraints;
    } 
  
  void addConstants() {
    for (auto x: variables)
      constraints += Constraint(elpair(elem("="),x),eltuple({x}));
  }
  
  //returns the set of constraint symbols which appear
  lset getConstraintSymbols()   
  {    
    lset ret;
    for (auto c: constraints) {
      lbool found = false;
      for (auto el:ret)
        If(el == c.first)
          found &= true;
      If (!found)
        ret += c.first;
    }  
    return ret;
  }
  
  
  //return the set of pairs (relation symbol s, set of tuples in constraint with symbol s)
  lsetof<lpair<elem,lset>> getRelations()
  {    
    lset res;
    lset symbols = getConstraintSymbols();
        
    for (auto sym: symbols) {
      lset curr = newSet();
      
      for (auto c: constraints)
        If(c.first == sym)
          curr += c.second;
      
      res += elpair(sym,curr);
    }
      
    return res;
  }
  
  //construct the cartesian power constraint graph, whose variables are n-tuples of variables,
  //and constraints are obtained by taking coordinatewise AND of n-tuples of constraints
  ConstraintGraph cartesianPower(const int & n) 
  {
    ConstraintGraph res;
    
    res.variables = ::cartesianPower(variables, n);
    res.constraints = newSet();
    
    lset ass = getRelations();
    //it would be nicer to have auto instead of lset but this results in error (cf. issue #7)    
    
    for (auto x: ass) {
      lset rel = asSet(as<elpair>(x).second);
      lset pow = ::cartesianPower(rel, n);

      for (auto tup: pow) {
        res.constraints += Constraint(as<elpair>(x).first, transpose(as<eltuple>(tup)));
      }
    }
    return res;
  }
  
  //replaces each vertex v by a pair (tag, v)
  void markVariables(elem tag)
  {
    variables = SETOF(elpair(tag,v), v:variables, true);
    
    lset new_constraints = newSet();
    
    for (auto c:constraints) {
      eltuple t;
      
      for (int i=0; i < c.second.size(); i++) {
        t.push_back(elpair(tag, c.second[i]));
      }
      
      new_constraints += Constraint(c.first, t);
    }
    
    constraints = new_constraints;    
  }
    
  ConstraintGraph squashedConstraintGraph()
  {
    return unarizedConstraintGraph().quotientConstraintGraph();    
  }
  
  ConstraintGraph quotientConstraintGraph(const symmetry& sym=LIN) 
    //Returns the finite ConstraintGraph obtained by quotienting I
    //by the equivalence relation which identifies two variables if they 
    //have the same atomic type with respect to "<"
  {   
    struct QuotientConstraintGraph {
      lsetof<Constraint> constraints = newSet();
  
      lset variable_orbits;
      
      void addConstraint(const elem& relname, const eltuple& inrel, 
        int inrel_pos = 0, vector<int> orbitsofar = vector<int> ()) {
        if(inrel_pos == inrel.size()) {
          eltuple tuple = vector<elem> ();
    
          for(int i=0; i<inrel_pos; i++)  
            tuple.push_back(elem(orbitsofar[i]));
    
          constraints += Constraint(elem(relname), tuple);
          return;
          }
  
        for(auto o: variable_orbits) {
          elpair p = as<elpair> (o);
          If(p.first == inrel[inrel_pos]) {
            auto orbitsofar2 = orbitsofar;
            orbitsofar2.push_back(as<int> (p.second));
            addConstraint(relname, inrel, inrel_pos+1, orbitsofar2);
            }
          }
        }
    } Q;
    
    Q.variable_orbits = computeOrbits(variables, sym);
    
    for(auto c: constraints)
      Q.addConstraint(c.first, as<eltuple> (c.second));

    int i=1;
    lset orbs = newSet();
  
    for (auto o: Q.variable_orbits)
      orbs += i++;
    
    ConstraintGraph R;
    
    R.variables = orbs;
    R.constraints = Q.constraints;
    
    return R;
  }
};

ostream& operator<<(ostream& os, const ConstraintGraph& I) 
{
  os
    << "The set of variables is:" << std::endl 
    << I.variables << std::endl << std::endl;

  os
    << "The set of constraints is:" << std::endl 
    << I.constraints << std::endl << std::endl;
  
  return os;
}


typedef pair<int, std::vector<int>> IdentityTerm;
typedef pair<IdentityTerm, IdentityTerm> IdentityTT;
typedef pair<IdentityTerm, int> IdentityTV;

eltuple mapTuple(std::vector<int> recipe, eltuple source) {
  eltuple res;
  
  for (int i=0; i<recipe.size(); i++)  
    res.push_back(source[recipe[i]]);
  
  return res;
}

//copies elements of source to target, without duplicates
void copyAsSet(const  std::vector<int> & source, std::vector<int> & target)
{
  for (int i=0; i<source.size(); i++) {
    bool found=false;
    int j;
    for (j=0; (!found) && (j<target.size()); j++) 
      if (target[j]==source[i]) 
        found=true;
    if (!found) 
      target.push_back(source[i]);
  }
}

//maps each element of source to its index at which it occurs in map
std::vector<int> applyMap(const std::vector<int> & source, std::vector<int> & map)
{  
  std::vector<int> res;
  for (int i=0; i<source.size(); i++) {
    for (int j=0; ; j++)
      if (map[j]==source[i]) {
        res.push_back(j);
        break;
      }
  }
  return res;
}
  


class CSPInstance {
public:
  ConstraintGraph source;
  ConstraintGraph target;

  CSPInstance(ConstraintGraph I, ConstraintGraph T) 
  {
    source = I;
    target = T;    
  }
  

//  constructs CSP instance corresponding to existence of terms satisfying h1 equations
CSPInstance(ConstraintGraph I, std::set<IdentityTT> identitiesTT, std::set<IdentityTV> identitiesTV)
{
  std::set<IdentityTerm> allTerms;

  for (auto i: identitiesTT) {
    allTerms.insert(i.first);
    allTerms.insert(i.second);
  }
  
  for (auto i:identitiesTV) {
    allTerms.insert(i.first);
  }

  std::map<int,int> arities;  //a map mapping a function symbol to its arity
  
  for (auto t:allTerms) {
    if (arities.find(t.first)==arities.end()) {
      arities[t.first]=t.second.size();
    }
    else if (arities[t.first]!=t.second.size())    
      throw exception();
  }    
  
  lset vars = newSet();
  lset cons = newSet();
  
  for (auto& x: arities) {
    ConstraintGraph cart = I.cartesianPower(x.second);
  cart.markVariables(elem(x.first));

    vars |= cart.variables;
    cons |= cart.constraints;
  }
    
  
  for (auto eq: identitiesTT) {
    IdentityTerm lhs = eq.first;
    IdentityTerm rhs = eq.second;
    
    int f = lhs.first;
    int g = rhs.first;
    
    auto args_f = lhs.second;
    auto args_g = rhs.second;
    
    
    vector<int> all_args;
    
    copyAsSet(args_f,all_args);
    copyAsSet(args_g,all_args);
        
    vector<int> args_f_rel=applyMap(args_f,all_args);
    vector<int> args_g_rel=applyMap(args_g,all_args);
    
    
    lset cart = cartesianPower(I.variables, all_args.size());
    
    for (auto t: cart) {
      auto ftup=elpair(elem(f), mapTuple(args_f_rel, as<eltuple>(t)));
      auto gtup=elpair(elem(g), mapTuple(args_g_rel, as<eltuple>(t)));
      
      If (ftup!=gtup)
        cons += Constraint(elem("="),eltuple({ftup,gtup}));
    }
  }
  
  for (auto eq: identitiesTV) {
    IdentityTerm lhs = eq.first;
    int x = eq.second;
    int f = lhs.first;
    auto args_f = lhs.second;
    
    vector<int> all_args;
    
    copyAsSet(args_f,all_args);
    vector<int> args_f_rel=applyMap(args_f,all_args);
    
    int n = arities[f];
    
    lset cart = cartesianPower(I.variables, all_args.size());
    
    for (auto t:cart)
      cons += Constraint(elpair(elem("="),as<eltuple>(t)[x]),
        eltuple({elpair(elem(f), mapTuple(args_f_rel, as<eltuple>(t)))}));  
  }
  
  source.variables = vars;
  source.constraints = cons;
  
  target = ConstraintGraph(I);
  target.addConstants();
}

  //test the existence of a homomorphism from source to target
  bool hasSolution() {
    return true;
  }    
};


int main() {
  lasttime = getVa();
  initLois();
  // pushSolverDiagnostic("checking: ");

//  solver = solverBasic() || solverIncremental("cvc4 --lang smt --incremental");
//  solver = solverBasic() || solverIncremental("./z3 -smt2 -in");
//  solver = solverBasic() || solverIncremental("./cvc4-new --lang smt --incremental");
    
    
  Domain dA("Atoms");
  lset A = dA.getSet();

  sym.neq = "≠";
      
  ConstraintGraph I; //the CSP ConstraintGraph
  
  // I.variables = SETOF (newSet(a,b), a:A, b:A, a!=b);
  // I.constraints = SETOF (Constraint(elem("NEQ"), eltuple({newSet(a,b), newSet(b,c)})), a:A, b:A, c:A, a!=c);
  I.variables = SETOF (a, a:A, true);
  I.constraints = SETOF (Constraint(elem("Dom"), eltuple({a})), a:A, true);
    
  std::cout << "THE INSTANCE: " << endl << I;
    
  ConstraintGraph Pow = I.cartesianPower(2);
  
  std::cout << "THE SQUARE:" << endl << Pow;
  
  
  CSPInstance CSP = CSPInstance(I,{IdentityTT(IdentityTerm(0,{0,1}),IdentityTerm(0,{1,0}))}, {IdentityTV(IdentityTerm(0,{0,0}),0)});
  
  std::cout << "THE TERM INSTANCE SOURCE:" << endl << CSP.source;
  std::cout << "THE TERM INSTANCE TARGET:" << endl << CSP.target;
  


  
/*  
  ConstraintGraph U = I.unarizedConstraintGraph();
    
  std::cout 
    << "The set of variables of the unarized ConstraintGraph is:" << std::endl 
    << U.variables << std::endl << std::endl;
  
  std::cout 
    << "The set of constraints of the unarized ConstraintGraph is:" << std::endl 
    << U.constraints << std::endl << std::endl;
  

  // std::cout << "orbits: " << computeOrbits(U.variables, EQ) << endl;
  
  ConstraintGraph S = U.quotientConstraintGraph(EQ);
  
  std::cout 
    << "The set of variables of the squashed ConstraintGraph is:" << std::endl 
    << S.variables << std::endl << std::endl;
  
  std::cout 
    << "The set of constraints of the squashed ConstraintGraph is:" << std::endl 
    << S.constraints << std::endl << std::endl;*/
  
  return 0;
  }

  