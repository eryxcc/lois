// automatically test the LOIS library.

// use SPASS?
// #define TEST_SPASS

#include <stdlib.h>
#include <sys/time.h>

#include "../include/loisextra.h"

#include <iostream>
using namespace std;
using namespace lois;

// an auxiliary function for showTimeElapsed.

long long getVa() {
  struct timeval tval;
  gettimeofday(&tval, NULL);
  return tval.tv_sec * 1000000 + tval.tv_usec;
  }

long long lasttime;

// measure the time since lasttime, which is updated by 
// showTimeElapsed() itself.

void showTimeElapsed() {
  long long curtime = getVa();
  
  if(curtime - lasttime < 100000) 
    cout << "time elapsed: " << (curtime-lasttime) << " \u03bcs" << endl;

  else
    printf("time elapsed: %.6f s\n", (curtime-lasttime)/1000000.0);

  lasttime = curtime;
  }

// assert whether the given rbool is: (ALWAYS) always satisfied, 
// (NEVER) never satisfied, or (MAYBE) sometimes satisfied and sometimes not.

#define MAYBE 3
#define ALWAYS 1
#define NEVER 2

string modenames[4] = {"inconsistent", "always", "never", "maybe"};

void testSat(int mode, rbool phi) {
  int cnt = 0;
  If(phi) cnt ++;
  If(!phi) cnt += 2;
  if(cnt == mode) cout << "passed: ";
  else cout << "failed: ";
  cout << modenames[cnt] << " " << phi << endl;
  if(cnt != mode) exit(1);
  }

// tree test.

// table used to build the tree structure
// strindex[a][b] = index of LCA of elements #a and #b

int strindex[16][16];

// strindex is calculated, draw the tree
// (not used in the autotest, we simply count the trees)

void drawTree(vector<term>& v, int from) {
  cout << v[from];
  int children = 0;

  for(size_t i=0; i<v.size(); i++) if(i>from && strindex[i][i] == from)
    cout << "=" << v[i];
    
  for(size_t i=0; i<v.size(); i++) 
    if(strindex[from][i] == from && strindex[i][i] != from && strindex[from][from] != i) {
      bool direct = true;
      for(size_t j=0; j<v.size(); j++) 
        if(i!=j && j != from && strindex[from][j] == from && strindex[i][j] == j)
          direct = false;
      if(direct) {
        if(!children) cout << "["; else cout << ", ";
        children++;
        drawTree(v, i);
        }
      }

  if(children) cout << "]";
  }

int numstruct = 0;

// build the strindex table 

void analyzeStructure(vector<term>& v, RelTree& tree, const lset& A, size_t a=0, size_t b=0) {

  if(b>a) a=a+1, b=0;
  
  if(a == v.size()) {
    int i = 0;
    for(int j=1; j<a; j++) i = strindex[i][j];
    // cout << "possible tree structure: "; drawTree(v, i); cout << endl;
    numstruct++;
    return;
    }

  term lca = tree.lca(v[a], v[b]);

  {
 
    rbool psi = ftrue;
    for(size_t i=0; i<v.size(); i++) {
      If(psi && v[i] == lca) {
        strindex[a][b] = strindex[b][a] = i;
        analyzeStructure(v, tree, A, a, b+1);
        }
      psi = psi && (v[i] != lca);
      }

    If(psi) {
      strindex[a][b] = strindex[b][a] = v.size();
      v.push_back(lca);
      analyzeStructure(v, tree, A, a, b+1);
      v.pop_back();
      }
    }

  }

int countStructures(vector<term>& v, RelTree& tree, const lset& A, rbool phi, int expected) {
  numstruct = 0;
  If(phi) analyzeStructure(v, tree, A);
  if(numstruct == expected)
    std::cout << "passed: found " << numstruct << " trees for " << phi << std::endl;
  else {
    std::cout << "failed: " << phi << " (" << numstruct << ")" << std::endl;
    exit(1);
    }
  showTimeElapsed(); cout<<endl;
  return numstruct;
  }

void testTree(const lset& A) {
  
  RelTree tree(sym.arrow, sym.notarrow, sym.min);
  
  lsetof<term> At(A);

  for(term x: At) for(term y: At) for(term z: At) for(term w: At) {
    x.asVar()->name = "X";
    y.asVar()->name = "Y";
    z.asVar()->name = "Z";
    w.asVar()->name = "W";
    vector<term> v;
    v.push_back(x);
    v.push_back(y);
    v.push_back(z);
    v.push_back(w);
    
    countStructures(v, tree, A, ftrue, 416);
    countStructures(v, tree, A, x != y && x != z && x != w && y != z && y != w && z != w, 262);

    countStructures(v, tree, A, z==w, 32);
    countStructures(v, tree, A, x != y && x != z && y != z && z == w, 22);
    }
  }

// Test some basic properties of orders. In particular, we check 
// whether a subset of A has exactly one supremum, unless empty or
// unbounded.

lelem max(lset X) { 
  lset answer = newSet();
  for(elem x: X) If(FORALL(y, X, x >= y)) answer += x;
  return extract(answer);
  }

lelem min(lset X) { 
  lset answer = newSet();
  for(elem x: X) If(FORALL(y, X, x <= y)) answer += x;
  return extract(answer);
  }

lelem supremum(lset X, lset domain) { 
  return min(FILTER(m, domain, FORALL(x, X, m >= x)));
  }

void testOrder(lset A) {

  for(elem a: A) for(elem b: A) for(elem c: A) {
    rbool phi = (a<b);
    cout << phi << endl;

    testSat(NEVER, (a<b) && (b<c) && (c<a));
    testSat(MAYBE, a<b && b<c); 
    testSat(MAYBE, c<b && b<a);
    testSat(MAYBE, b<a && a<c); 
    testSat(MAYBE, c<a && a<b);
    testSat(MAYBE, a<c && c<b);
    testSat(MAYBE, b<c && c<a);

    lset three = newSet({a,b,c});
    lelem mx = max(three);
    cout << "max(" << three << ") = " << mx << endl;

    lelem sup = supremum(three, A);
    cout << "supremum = " << sup << endl;
    testSat(ALWAYS, cardinality(newSet(sup)) == 1);
    testSat(ALWAYS, sup == mx);
    }

  for(elem a: A) for(elem b: A) If(a<b) {
    lelem sup = supremum(newSet(a,b), A);
    cout << "sup(a,b) = " << sup << endl;
    testSat(ALWAYS, cardinality(newSet(sup)) == 1);
    }
  
  for(elem a: A) for(elem b: A) If(a<b) {
    lset interval = FILTER(z, A, a<z && z<b);
    lelem sup = supremum(interval, A);
    cout << "sup(interval) = " << sup << endl;
    testSat(ALWAYS, cardinality(newSet(sup)) == 1);
//  testSat(ALWAYS, card(sup) == 1);
    }

  lelem sup = supremum(A, A);

  cout << "sup(A) = " << sup << endl;
  testSat(ALWAYS, cardinality(newSet(sup)) == 0);
  }

// BFS on the random bipartite graph, as advertised in the paper.

void testRandomBipartite(lset A) {
  string R = "R";

  RelBinary edge(sym.edge, sym.noedge, lmNoLoops, smSymmetric);
  
  RelUnary partition(R, 2);

  // for some reason 'edge' is passed as const to lambda
  RelBinary *ebi = &edge; 

  for(elem x: A) for(elem y: A) {
    testSat(NEVER, !EXISTS(z, A, (*ebi)(x,z) && (*ebi)(y,z)));
    }

  testSat(ALWAYS, FORALL(x, A, !(*ebi)(x,x)) );

  lsetof<elpair> E;

  for(elem x: A) for(elem y: A) 
    If(edge(x,y) && !partition.together(x,y))
      E += elpair(x,y);

  cout << "E = " << E << endl;

  for(elem a: A) {
   
    lnum<int> iterations = 0;
    
    lset R = newSet(a);
    lset P = R;
    
    While(R != A) {
      iterations++;
      lset NP;
      cout << "R= " << R << endl;
      for(elem x: P) 
        for(elpair ed: E) 
          If(ed.first == x && !memberof(ed.second, R))
            R += ed.second, NP += ed.second;
      cout << "R= " << R << endl;
      P = NP;
      }
  
    cout << "iterations = " << iterations << endl;
    testSat(ALWAYS, iterations == 3);
    }

  }

// test whether the 'queue' semantics of the 'for' loop works, as
// advertised in the paper.

void testQueue() {
  lsetof<int> N;
  N += 0;
  for(int n: N) if(n < 10) N += (n+1);
  cout << N << endl;
  testSat(ALWAYS, cardinality(N) == 11);
  }

// test whether the "-=" operator works in the natural way, as
// advertised in the paper.

void testRemoval(lset A) {
  lset X, Y;
  
  for(elem a: A) for(elem b: A) {
    X += elpair(a, b);
    If(a != b) Y += elpair(a, b);
    }
  
  for(elem a: A) X -= elpair(a, a);
  
  testSat(ALWAYS, X == Y);
  }

// test whether an improper assignment throws an exception, as
// mentioned in the paper.

void testAssignment(lset A) {
  lbool phi;
  
  try {
    for(elem a: A) for(elem b: A) phi = (a==b);
  } catch(assignment_exception e) {
    printf("passed: assignment_exception\n");
    }
  }

// calculate the powerset.
lset orbitpowerset(lset X) {
  lset result;
  for(elem x: X) {
    lset orbit = getorbit(x, {});
    lset Y = X &~ orbit;
    lset P1 = orbitpowerset(Y);
    for(elem p: P1) {
      result += p;
      result += (asSet(p) | orbit);
      }
    return result;
    }
  return newSet(emptyset);
  }

void testPowerset(lset A) {
  cout << "orbit powerset of " << A << " is:" << endl << orbitpowerset(A) << endl;
  }
  

void testFullyPseudoparallel(lset A) {
  lset X;
  for(elem a: A) X += a;
  for(elem a: A) for(elem b: A) If(a != b) X += elpair(a,b);
  X += 42;
  cout << "X is: " << X << endl;
  cout << "singleton set is: " << getsingletonset(X) << endl;
  
  for(lelem e: fullypseudoparallel(X)) {
    cout << "e = " << e << endl;
    }
  }

int main() {
  lasttime = getVa();
  initLois();
  // pushSolverDiagnostic("checking: ");
  
#ifdef TEST_SPASS
  useDefaultSolver();
  pushSolverSPASS();
  pushSolverBasic();
  
  sym.edge = "edge";
  sym.greater = "gt";
  
#endif

  Domain dA("A");
  lset A = dA.getSet(); 

  testFullyPseudoparallel(A);
  showTimeElapsed(); cout<<endl;

  testAssignment(A);
  showTimeElapsed(); cout<<endl;
  
  testRandomBipartite(A);
  showTimeElapsed(); cout<<endl;
  
  testOrder(A);
  showTimeElapsed(); cout<<endl;

  testQueue();
  showTimeElapsed(); cout<<endl;

  testRemoval(A);
  showTimeElapsed(); cout<<endl;
  
  testPowerset(A * A);
  showTimeElapsed(); cout<<endl;

  testTree(A);
  
  return 0;
  }
