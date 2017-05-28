// measure various solvers.

#include <stdlib.h>
#include <sys/time.h>
#include "../include/loisextra.h"
using namespace lois;
using namespace std;

namespace lois {
extern string lastsolver;
}

// testSat tool.

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

Domain *getDomain(const lsetof<term>& A) {
  for(term x:A) return x.asVar()->getDom();
  }

// Test some basic properties of orders. In particular, we check 
// whether a subset of A has exactly one supremum, unless empty or
// unbounded.

lelemof<term> max(const lsetof<term>& X) { 
  lsetof<term> answer;
  for(term x: X) If(FORALL(y, X, x >= y)) answer += x;
  return extract(answer);
  }

lelemof<term> min(const lsetof<term>& X) { 
  lsetof<term> answer;
  for(term x: X) If(FORALL(y, X, x <= y)) answer += x;
  return extract(answer);
  }

lelemof<term> supremum(const lsetof<term>& X, const lsetof<term>& domain) { 
  return min(FILTER(m, domain, FORALL(x, X, m >= x), term));
  }

void testOrder(const lsetof<term>& A) {

  for(term a: A) for(term b: A) for(term c: A) {
    rbool phi = (a<b);
    cout << phi << endl;

    testSat(NEVER, (a<b) && (b<c) && (c<a));
    testSat(MAYBE, a<b && b<c); 
    testSat(MAYBE, c<b && b<a);
    testSat(MAYBE, b<a && a<c); 
    testSat(MAYBE, c<a && a<b);
    testSat(MAYBE, a<c && c<b);
    testSat(MAYBE, b<c && c<a);

    lsetof<term> three = newSet<term>({a,b,c});
    lelemof<term> mx = max(three);
    // cout << "max(" << three << ") = " << mx << endl;

    lelemof<term> sup = supremum(three, A);
    // cout << "supremum = " << sup << endl;
    testSat(ALWAYS, cardinality(newSet(sup)) == 1);
    testSat(ALWAYS, sup == mx);
    }

  for(term a: A) for(term b: A) If(a<b) {
    lelemof<term> sup = supremum(newSet(a,b), A);
    // cout << "sup(a,b) = " << sup << endl;
    testSat(ALWAYS, cardinality(newSet(sup)) == 1);
    }
  
  for(term a: A) for(term b: A) If(a<b) {
    lsetof<term> interval = FILTER(z, A, a<z && z<b, term);
    lelemof<term> sup = supremum(interval, A);
    // cout << "sup(interval) = " << sup << endl;
    testSat(ALWAYS, cardinality(newSet(sup)) == 1);
//  testSat(ALWAYS, card(sup) == 1);
    }

  lelemof<term> sup = supremum(A, A);
  
  // cout << "sup(A) = " << sup << endl;
  testSat(ALWAYS, cardinality(newSet(sup)) == 0);
  }

// test whether the 'queue' semantics of the 'for' loop works, as
// advertised in the paper.

void testQueue(const lsetof<term>& A) {
  lsetof<int> N;
  N += 0;
  for(int n: N) if(n < 10) N += (n+1);
  cout << N << endl;
  testSat(ALWAYS, cardinality(N) == 11);
  }

// test whether the "-=" operator works in the natural way, as
// advertised in the paper.

void testRemoval(const lsetof<term>& A) {
  lsetof<lpair<term,term>> X, Y;
  
  for(term a: A) for(term b: A) {
    X += make_lpair(a, b);
    If(a != b) Y += make_lpair(a, b);
    }
  
  for(term a: A) X -= make_lpair(a, a);
  
  testSat(ALWAYS, X == Y);
  }

// test whether an improper assignment throws an exception, as
// mentioned in the paper.

void testAssignment(const lsetof<term>& A) {
  lbool phi;
  
  try {
    for(term a: A) for(term b: A) phi = (a==b);
  } catch(assignment_exception e) {
    printf("passed: assignment_exception\n");
    }
  }

// reachability test

template<class A> lsetof<A> reach(const lsetof<lpair<A,A>>& E, const lsetof<A>& S) {
  lsetof<A> R = S;
  lsetof<A> P;
  While (P!=R) {
    P = R;
    for (auto e: E) 
      If (memberof(e.first, R))
        R += e.second;
  }
  return R;
}

void testReachable(const lsetof<term>& A) {

  typedef lpair<term, term> vertex;

  lsetof<vertex> Pairs;
     
  for(term a: A) for(term b: A) 
    Pairs += make_lpair(a,b);

  lsetof<lpair<vertex, vertex>> E;
  
  for(auto x: Pairs) 
    for (auto y: Pairs)
      If ((x.second == y.first)
        && ((y.second != x.first)))
        E += make_lpair(x,y);
  
  cout << E << endl;

  for(auto a: A) for(auto b: A) If(a!=b) {

    auto s = make_lpair(a,b);
    auto t = make_lpair(b,a);
    
    If (memberof(t, reach(E, newSet(s))))
      cout << "Reached" << endl;
    }
  }

// minimalization test
// not yet updated for static typing

/*
void minimalize1(lset Q, lset F, rsetof<eltuple> delta, lset alph) {
  cout << "Q = " << Q << endl;
  cout << "F = " << F << endl;
  cout << "delta = " << delta << endl;
  
  lset NF = Q &~ F;
  cout << "NF = " << NF << endl << endl;
  
  lset eq = (F * F) | (NF * NF);
  
  int it = 0;
  
  lbool cont = true;
  
  While(cont) {
    cout << "eq = " << eq << endl << endl;
    
    lset neq;
    
    for(term q1: Q) for(term q2: Q) If(memberof(elpair(q1,q2), eq)) {
      lbool b = true;
      for(term x: alph) for(eltuple t1: delta) for(eltuple t2: delta)
        If(t1[0] == q1 && t1[1] == x && t2[0] == q2 && t2[1] == x)
          If(!memberof(elpair(t1[2], t2[2]), eq))
            b = ffalse;
      If(b) neq += elpair(q1, q2);
      }
    
    If(eq == neq) cont = ffalse;
    eq = neq; it++;
    }

  cout << "number of iterations: " << it << endl << endl;
  
  lset classes;
  
  for(term a: Q) {
    lset t;    
    for(term b: Q) If(memberof(elpair(a,b), eq)) 
      t += b;
    classes += t;
    }
  
  cout << "unoptimized classes: " << classes << endl;
  classes = optimize(classes);

  cout << "classes: " << classes << endl;  
  }

void minimalize2(lset Q, lset F, rsetof<eltuple> delta, lset alph) {
  cout << "Q = " << Q << endl;
  cout << "F = " << F << endl;
  cout << "delta = " << delta << endl;
  
  lset classes;
  classes += F;
  classes += Q &~ F;
  
  int it = 0;
  
  lbool cont = true;
  
  While(cont) {
    lset nclasses;
    cout << "classes = " << classes << endl << endl;
    cont = ffalse;
    
    for(term EC: classes) {
      for(term q1: EC) {
        lset EC1;
        for(term q2: EC) {
          lbool b = true;
          for(eltuple t1: delta) for(eltuple t2: delta) 
            If(t1[0] == q1 && t2[0] == q2 && t1[1] == t2[1])
            If(EXISTS(c, classes, memberof(t1[2], asSet(c)) && !memberof(t2[2], asSet(c))))
              b &= ffalse;
          Ife(b) EC1 += q2;
          else cont = ftrue;
          }
        nclasses += EC1;
        }
      }

    classes = optimize(nclasses); it++;
    }

  cout << "number of iterations: " << it << endl << endl;
  
  cout << "classes: " << classes << endl;  
  }

// as in the SAT paper
void minimalize3(lset Q, lset F, rsetof<eltuple> delta, lset alph) {
  
  lsetof<elpair> E;
  for(term p: Q) for(term q: Q) for(term a: alph) 
    for(term pbis: Q) for(term qbis: Q) 
       If(memberof(eltuple({p,a,pbis}), delta) && memberof(eltuple({q,a,qbis}), delta))
         E += elpair(elpair(pbis,qbis), elpair(p,q));

  lset S = (F * (Q&~F)) | ((Q&~F) * F);
  lset equiv = Q*Q &~ reach(E, S);
  
  lset classes;
  
  for(term q: Q) {
    lset qclass;
    for(term p: Q) If(memberof(elpair(p,q), equiv)) qclass += p;
    classes += qclass;
    }
  
  cout << "classes: " << classes << endl;  
  }

void mtestA(lsetof<term>& A, int id) {
  lset Q;
  lset F;
  lsetof<eltuple> delta;
  
  lset R;

//   for (term x:A) for (term y:A) {
//     elem u=elpair(x,y);
//     If (x==y)
// 		R+= elpair(x,y);
// 	If (x!=y)
//       R+=eltuple({x,x,y});
// //    R += u;
//   }
//   cout << "R= " << R << endl;
    

  lset alph = A;

  Q += 0;
  for(term a: A) Q += eltuple({a});
  for(term a: A) for(term b: A) Q += eltuple({a, b});
  for(term a: A) for(term b: A) for(term c: A) Q += eltuple({a, b, c});
  Q += 4;
  
  for(term a: A) for(term b: A) for(term c: A) 
    If(a==b || a==c || b==c) F += eltuple({a,b,c});
    
  for(term x: A) {
    delta += eltuple({elof(0),x,elof(eltuple({x}))});

    for(term a: A)
      delta += eltuple({elof(eltuple({a})),x,elof(eltuple({a,x}))});

    for(term a: A) for(term b: A)
      delta += eltuple({elof(eltuple({a,b})),x,elof(eltuple({a,b,x}))});

    for(term a: A) for(term b: A) for(term c: A)
      delta += eltuple({elof(eltuple({a,b,c})),x,elof(4)});

    delta += eltuple({elof(4),x,elof(4)});
    }

  if(id == 1) minimalize1(Q, F, delta, alph);
  if(id == 2) minimalize2(Q, F, delta, alph);
  if(id == 3) minimalize3(Q, F, delta, alph);
  }

void mtestB(lsetof<term>& A, int id) {
  RelInt real(sym.greater, sym.leq, sym.max, sym.min, sym.plus, sym.times, sym.minus, sym.divide);
  pushorder po(&real);
  Domain *d = getDomain(A);

  term r0 = real.constant(d, 0);
  term r3 = real.constant(d, 3);

  lset Nat = FILTER(x, A, x >= r0, term);
  
  lset alph = Nat; alph += '$';
  lset Q = Nat * Nat;
  lset F = newSet(r3) * Nat;
  lsetof<eltuple> delta;
  for(term p: Q) {
    term m = as<term>(as<elpair> (p).first);
    term n = as<term>(as<elpair> (p).second);
    delta += eltuple({p, '$', elpair(m, r0)});
    for(term x: Nat)
      delta += eltuple({p, x, elpair(real.max(m, real.plus(n, as<term>(x))), real.plus(n, as<term>(x)))});
    }
  
  if(id == 1) minimalize1(Q, F, delta, alph);
  if(id == 2) minimalize2(Q, F, delta, alph);
  if(id == 3) minimalize3(Q, F, delta, alph);
  }

void testMinimizeA1(lset A) { mtestA(A, 1); }
                 
void testMinimizeA2(lset A) { mtestA(A, 2); }

void testMinimizeA3(lset A) { mtestA(A, 3); }

void testMinimizeB1(lset A) { mtestB(A, 1); }

void testMinimizeB2(lset A) { mtestB(A, 2); }

void testMinimizeB3(lset A) { mtestB(A, 3); }
*/

// testReal

void testInt(lsetof<term>& A) { 
  string symg = "gt";
  RelInt rint(symg, sym.leq, sym.max, sym.min, sym.plus, sym.times, sym.minus, sym.divide);
  Domain *D = getDomain(A);

  RelOrder *ord = mainOrder;
  using namespace orderedfield_ops;
  mainOrder = mainField = &rint; mainDomain = getDomain(A);

  term c1 = rint.constant(mainDomain, 1);
  term c2 = rint.constant(mainDomain, 2);
  term c5 = rint.constant(mainDomain, 5);
  term c8 = rint.constant(mainDomain, 8);
  term c9 = rint.constant(mainDomain, 9);

  lsetof<term> odd;
  for(term x: A) odd += c1 + c2 * x;
  lsetof<term> z;
  for(term x: odd) If(c5 < x) z += x + c1;
  cout << z << endl;
  If(memberof(c8, z)) cout << "eight in the set" << endl;
  If(memberof(c9, z)) cout << "nine in the set" << endl;

  mainOrder = ord;
  }

void testReal(lsetof<term>& A) { 
  string symg = "gt";
  RelReal real(symg, sym.leq, sym.max, sym.min, sym.plus, sym.times, sym.minus, sym.divide);
  pushorder po(&real);
  Domain *d = getDomain(A);

  RelOrder *ord = mainOrder;
  using namespace orderedfield_ops;
  mainOrder = mainField = &real; mainDomain = getDomain(A);

  term c2 = real.constant(mainDomain, 2);
  term c4 = real.constant(mainDomain, 4);
  term c5 = real.constant(mainDomain, 5);

  for(term x: A) {
    lbool b = false;

    b |= (x >= c4 && x <= c5 - c2);

    cout << "b = " << b << endl;

    If(b) cout << "b can be true" << endl;
    If(!b) cout << "b can be false" << endl;
    }

  mainOrder = ord;
  }

// packings on an interval

lsetof<lsetof<term>> packings(lsetof<term>& F, RelReal& real, term radius) {
  lsetof<lsetof<term>> done;
  lsetof<lsetof<term>> active = newSet(newSet<term>({}));

  for (auto P: active) {    
    cout << "Considering " << P << " where " << emptycontext << endl;
    lsetof<term> N;
        
    for (term x: P) for(term y: F)
      If(y > real.minus(x, radius) && y < real.plus(x, radius))
        N += y;
    
    Ife (F == N)
      done += P;
    else
      for (term x: F) If(!memberof(x, N)) {
        lsetof<term> newP = P | newSet<term>(x);
        active += newP;
        }

    }

  return done;
  } 

void testPacking(lsetof<term>& A) {
  RelReal real(sym.greater, sym.leq, sym.max, sym.min, sym.plus, sym.times, sym.minus, sym.divide);
  pushorder po(&real);
  Domain *d = getDomain(A);

  lsetof<term> F;

  RelOrder *ord = mainOrder;
  using namespace orderedfield_ops;
  mainOrder = mainField = &real; mainDomain = getDomain(A);

  term c0 = real.constant(mainDomain, 0);
  term c1 = real.constant(mainDomain, 1);
  term c5 = real.constant(mainDomain, 5);

  for(term x: A) If(x>c0 && x<c5) F += x;

  cout << "F = " << F << endl;

  auto pack = packings(F, real, c1);

  cout << "packings = " << pack << endl;

  for(int i=0; i<9; i++)
    for(auto P: pack) {
      If(cardinality(P) == i)
        cout << "There is a packing of size " << i << endl;
      }

  mainOrder = ord;
  }

// packings of circles

// Requires NRA, and z3 does not know whether it is possible to fit two points in
// distance at least 3 in a circle of radius 1 anyway...

typedef lpair<term,term> point;

term getX(point p) { return p.first; }
term getY(point p) { return p.second; }

lsetof<lsetof<point>> circlepackings(lsetof<point> F, RelReal& real, term mindist) {
  lsetof<lsetof<point>> done;
  lsetof<lsetof<point>> active;
  active += newSet<point>({});

  for (auto P: active) {
    cout << "Considering " << P << " where " << emptycontext << endl;
    lsetof<point> N;
        
    for (point eP: P) for(point eF: F) {            

      term dX = real.minus(getX(eP), getX(eF));
      term dY = real.minus(getY(eP), getY(eF));

      term distsq = real.plus(real.times(dX, dX), real.times(dY, dY));

      If(distsq < real.times(mindist,mindist)) N += eF;
      }

    Ife (F == N)
      done += P;
    
    for (point x: F) If(!memberof(x, N)) {
      lsetof<point> newP = P;
      newP += x;
      active += newP;
      }

    }

  return done;
  } 

void testCirclePacking(const lsetof<term>& A) {
  RelReal real(sym.greater, sym.leq, sym.max, sym.min, sym.plus, sym.times, sym.minus, sym.divide);
  pushorder po(&real);
  Domain *d = getDomain(A);

  lsetof<point> F;

  term mindist = real.constant(d, 3);
  term bigradius = real.constant(d, 1);

  for(term x: A) for(term y: A) {
    If(real.plus(real.times(x, x), real.times(y, y)) <= real.times(bigradius, bigradius))
      F += make_lpair(x, y);
    }

  cout << "F = " << F << endl;

  auto pack = circlepackings(F, real, mindist);

  cout << "circle packings = " << pack << endl;

  for(int i=0; i<19; i++)
    for(auto P: pack) {
      If(cardinality(P) == i)
        cout << "There is a circle packing of size " << i << endl;
      }
  }

// tools

long long getVa() {
  struct timeval tval;
  gettimeofday(&tval, NULL);
  return tval.tv_sec * 1000000 + tval.tv_usec;
  }

typedef void loistest(const lsetof<term>& A);

std::string whichSolver;

FILE *outtable, *outparse;

long long elapsed[9];
int queries[9];

solveptr fallback;

const char *answers[3] = {"unsat", "unknown", "sat"};

struct SolverSol : Solver {
  solveptr inner;

  virtual int solve(rbool phi) {
    long long ti = getVa();
    int t = inner->solve(phi);
    ti = getVa() - ti;
    int t0 = t;
    if(t == 1)
      t = fallback->solve(phi);
    queries[t0*3+t]++;
    elapsed[t0*3+t]+=ti;
    fprintf(outparse, "answer(%s)\n", answers[t]);
    return t;
    }
  
  virtual int solveEnv() {
    long long ti = getVa();
    int t = inner->solveEnv();
    ti = getVa() - ti;
    int t0 = t;
    if(t == 1)
      t = fallback->solveEnv();
    queries[t0*3+t]++;
    elapsed[t0*3+t]+=ti;
    fprintf(outparse, "answer(%s)\n", answers[t]);
    return t;
    }
  
  SolverSol(solveptr x) { 
    inner = x; }
  };

// time limit per query, in milliseconds
#define TLIM "30000"

// total time limit, in milliseconds
#define TOTALLIM(x) // x "180000"

void testSolvers(const char* name, loistest test) {
  if(name == NULL) {
    fprintf(outtable, "%-20s", "solver used");
    }
  else {
    fprintf(outtable, "%-20s", name);
    fflush(outtable);
    }
  
  fallback = 
    solverExhaustive(1000000000, false) ||
    solverCrash();
  
  string cmdline;

  for(int i=0; i<5; i++) {
    switch(i) {
      case 0: 
        whichSolver = "internal";
        solver = make_shared<SolverSol> (solverNamed("internal", solverExhaustive(1000000000, false)));
        break;
      
      case 1:
        whichSolver = "Z3";
        cmdline = "z3-*/bin/z3 -smt2 -in -t:" TLIM TOTALLIM(" -T:");
        if(name) solver = 
          solverBasic() || make_shared<SolverSol> (solverIncremental(cmdline.c_str()));
        break;
          
      case 2:
        whichSolver = "CVC4";
        // cmdline = "cvc4 --lang smt --incremental --tlimit-per=" TLIM TOTALLIM(" --tlimit=");
        cmdline = "./cvc4-2016-* --lang smt --incremental --tlimit-per=" TLIM TOTALLIM(" --tlimit=");
        if(name) solver = 
          solverBasic() || make_shared<SolverSol> (solverIncremental(cmdline));
        break;
      
      case 3:
        whichSolver = "CVC4*";
        cmdline = "cvc4 --lang smt --incremental --finite-model-find --tlimit-per=" TLIM TOTALLIM(" -tlimit=");
        if(name) solver = 
          solverBasic() || make_shared<SolverSol> (
            solverIncremental(cmdline)
            );
        break;

      case 4:
        whichSolver = "SPASS";
        cmdline = "SPASS -TimeLimit=30";
        if(name) solver = 
          solverBasic() || make_shared<SolverSol> (
            solverSPASS(cmdline + " query.spass > query.out")
            );
        break;
      }
    
    if(name == NULL) {
      /* if(i == 0) fprintf(outtable,"#q ");
      fprintf(outtable, " | %-10s", whichSolver.c_str());
      fflush(outtable); */
      fprintf(outparse, "testSolver(\"%s\", \"%s\");\n", whichSolver.c_str(), cmdline.c_str());
      fflush(outparse);
      continue;
      }
    
    else {
      fprintf(outparse, "solverfile(\"%s\", \"\%s\", \"%s\");\n",
        whichSolver.c_str(), name, lastsolver.c_str());
      fflush(outparse);
      }

    printf("TEST %s, SOLVER %s\n", name, whichSolver.c_str());
    fflush(stdout);
    Domain dA("A");
    auto A = dA.getSet();

    for(int i=0; i<9; i++) elapsed[i] = queries[i] = 0;

    try {  
      test(A);
    } catch(unsolvable_exception e) {
      fprintf(outparse, "unsolvable()\n"); fflush(outparse);
    }

    if(i == 0) {
//      fprintf(outtable,"%3d", queries);
      // fprintf(outparse, "queries(\"%s\", %d);\n", name, queries);
      // fflush(outparse);
      }


//    fprintf(outtable, " |%7d/%3d", int(elapsed/1000), int(bugcount));
//    fflush(outtable);

//     fprintf(outparse, "result(\"%s\", \"%s\", %Ld, %d);\n", name, whichSolver.c_str(), elapsed, bugcount);
//     fflush(outparse);
    for(int i=0; i<9; i++) 
       fprintf(outparse, "result(\"%s\", \"%s\", %d, %Ld, %d);\n", name, whichSolver.c_str(), i, elapsed[i], queries[i]);
    fflush(outparse);

//     printf("ELAPSED = %d, BUG = %d\n", int(elapsed), int(bugcount));
    solver = fallback;
    }
  fprintf(outtable, "\n");
  }

int main() {
  outtable = fopen("out/soltable.txt", "wt");
  outparse = fopen("out/solparse.txt", "wt");
  initLois();

  testSolvers(NULL, testOrder);
  
  smtLogic = "LIA";
  // testSolvers("int", testInt);

  /*
  testSolvers("real", testReal);
  testSolvers("packing", testPacking);

  smtLogic = "NRA";
  testSolvers("circle packing", testCirclePacking);

  smtLogic = "LIA";
  testSolvers("minimizeB1", testMinimizeB1);
  testSolvers("minimizeB2", testMinimizeB2);
  testSolvers("minimizeB3", testMinimizeB3); */

  smtLogic = "LRA";
  testSolvers("order", testOrder);
  testSolvers("reachable", testReachable);

/*
  testSolvers("minimizeA3", testMinimizeA3);
  testSolvers("minimizeA1", testMinimizeA1);
  testSolvers("minimizeA2", testMinimizeA2);
*/

  fclose(outtable);
  fclose(outparse);
  return 0;
  }
