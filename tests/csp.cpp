// Locally Finite Constraint Satisfaction Problems
// Bartek Klin, Eryk Kopczyski, Joanna Ochremiak, Szymon Toruczyk 
  
#include <stdlib.h>
#include <sys/time.h>

#include "../include/loisextra.h"

#include <iostream>
using namespace std;
using namespace lois;

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

void listAllChoices(lset& orbitids, std::vector<vptr>& variables, elem i, int a, int b, int& orbitcount) {
  b++;
  if(a == b) { a++; b=0; }
  if(a == (int) variables.size()) {
    orbitcount++;
    orbitids += elpair(i, orbitcount);
    return;
    }
  
  term va(variables[a]);
  term vb(variables[b]);
  
  If(va < vb) 
    listAllChoices(orbitids, variables, i, a, b, orbitcount);

  If(va > vb)
    listAllChoices(orbitids, variables, i, a, b, orbitcount);

  If(va == vb)
    listAllChoices(orbitids, variables, i, a, b, orbitcount);
  }

void addRelation(lset& orbitids, const char* relname, const eltuple& inrel, 
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
      addRelation(orbitids, relname, inrel, inrel_pos+1, orbitsofar2);
      }
    }
  }

void computeOrbits(const lset& Indices, lset& orbitids, int& orbitcount) {
  orbitids = emptyset;
  orbitcount = 0;
  for(elem i: Indices) {
    std::vector<vptr> variables;

    contextptr e = currentcontext;
    while(e != emptycontext) {
      if(!e) throw env_exception();
      for(auto w: e->var) variables.push_back(w);
      e = e->parent;
      }
    
    listAllChoices(orbitids, variables, i, 1, -1, orbitcount);
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
  
  lset Indices;
  for(elem a: A) for(elem b: A)
    // If(a != b)
      Indices += elpair(a, b);

  std::cout 
    << "The set of variable indices is:" << std::endl 
    << Indices << std::endl << std::endl;
  
  lset orbitids;
  int orbitcount;  
  computeOrbits(Indices, orbitids, orbitcount);
  
  std::cout 
    << "Orbit IDs for each variable index: " << std::endl 
    << orbitids << std::endl << std::endl;

  std::cout << "Number of orbits: " << orbitcount << std::endl;
  
  for(elem a: A) for(elem b: A) 
//    If(a != b) <- put or remove, to make it satisfiable or not
    {
      addRelation(orbitids, "NEQ", eltuple ({elpair(a,b), elpair(b,a)}));
      }

  return 0;
  }
