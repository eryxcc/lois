// Locally Finite Constraint Satisfaction Problems
// Bartek Klin, Eryk Kopczyski, Joanna Ochremiak, Szymon Toruczyk 
  
#include <stdlib.h>
#include <sys/time.h>

#include <set>

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
typedef elpair Constraint;

struct ConstraintGraph 
  //a constraint graph is a set of ``variables'' and a set of ``constraints'',
  //where each constraint is a pair (id, t) where id is and  and t is a tuple of variables.
{
  lset variables = newSet();
  lsetof<Constraint> constraints = newSet();
  
/*  ConstraintGraph powerConstraintGraph(int n)
    //Computes the power ConstraintGraph I^n, whose variables are n-tuples of variables of I,
    //and whose constraints are obtained by n-ary conjunction
  {
    eltuple tup_vars = new vector<lset>();
    eltuple tup_constraints = new vector<lset>();
    
    for (int i=0; i<n; i++)
      tup_vars.push_back(variables);
    
    
    
    ConstraintGraph res;
    res.variables = cartesian(tup);
    
    res.constraints = 
    
  }*/
  

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
  
  
    //constructs CSP instance corresponding to existence of terms satisfying h1 equations
  // CSPInstance TermInstance(set<string> fsym, set<string> vsym,
  // set<pair<pair<string,vect<string>>,string> fv,
  // set<pair<pair<string,vect<string>>,pair<pair<string,vect<string>>> ff
  // )
  // {
  //
  //
  // }
    
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

struct CSPInstance = {
  ConstraintGraph Instance;
  ConstraintGraph Template;
  
  CSPInstance(ConstraintGraph I, ConstraintGraph T) 
  {
    Instance = I;
    Template = T;    
  }
  
  //test the existence of a homomorphism from Instance to Template
  bool hasSolution() {
    return true;
  }    
}


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
  
  I.variables = SETOF (newSet(a,b), a:A, b:A, a!=b);
/*  for(elem a: A) for(elem b: A)
    If(a != b)
      I.variables += newSet(a, b);*/



  
  for(elem a: A) for(elem b: A)  for(elem c: A)
    If(a != c) 
    {
      I.constraints += Constraint(1, eltuple ({newSet(a,b), newSet(b,c)}));
    }
  
/*    set<int> myset = I.constraintSizes();
    myset.insert(1);
    cout << powers(I.variables, myset);*/
    
    
/*  std::cout 
    << "The set of variables is:" << std::endl 
    << I.variables << std::endl << std::endl;
  
  std::cout 
    << "The set of constraints is:" << std::endl 
    << I.constraints << std::endl << std::endl;
  
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

  