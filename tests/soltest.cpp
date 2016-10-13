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

Domain *getDomain(rset A) {
  for(elem x:A) return as<term>(x).asVar()->getDom();
  }

// Test some basic properties of orders. In particular, we check 
// whether a subset of A has exactly one supremum, unless empty or
// unbounded.

lelem max(rset X) { 
  lset answer = newSet();
  for(elem x: X) If(FORALL(y, X, x >= y)) answer += x;
  return extract(answer);
  }

lelem min(rset X) { 
  lset answer = newSet();
  for(elem x: X) If(FORALL(y, X, x <= y)) answer += x;
  return extract(answer);
  }

lelem supremum(rset X, rset domain) { 
  return min(FILTER(m, domain, FORALL(x, X, m >= x)));
  }

void testOrder(rset A) {

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

    rset three = newSet({a,b,c});
    lelem mx = max(three);
    // cout << "max(" << three << ") = " << mx << endl;

    lelem sup = supremum(three, A);
    // cout << "supremum = " << sup << endl;
    testSat(ALWAYS, cardinality(newSet(sup)) == 1);
    testSat(ALWAYS, sup == mx);
    }

  for(elem a: A) for(elem b: A) If(a<b) {
    lelem sup = supremum(newSet(a,b), A);
    // cout << "sup(a,b) = " << sup << endl;
    testSat(ALWAYS, cardinality(newSet(sup)) == 1);
    }
  
  for(elem a: A) for(elem b: A) If(a<b) {
    rset interval = FILTER(z, A, a<z && z<b);
    lelem sup = supremum(interval, A);
    // cout << "sup(interval) = " << sup << endl;
    testSat(ALWAYS, cardinality(newSet(sup)) == 1);
//  testSat(ALWAYS, card(sup) == 1);
    }

  lelem sup = supremum(A, A);
  
  // cout << "sup(A) = " << sup << endl;
  testSat(ALWAYS, cardinality(newSet(sup)) == 0);
  }

// test whether the 'queue' semantics of the 'for' loop works, as
// advertised in the paper.

void testQueue(rset A) {
  lsetof<int> N;
  N += 0;
  for(int n: N) if(n < 10) N += (n+1);
  cout << N << endl;
  testSat(ALWAYS, cardinality(N) == 11);
  }

// test whether the "-=" operator works in the natural way, as
// advertised in the paper.

void testRemoval(rset A) {
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

void testAssignment(rset A) {
  lbool phi;
  
  try {
    for(elem a: A) for(elem b: A) phi = (a==b);
  } catch(assignment_exception e) {
    printf("passed: assignment_exception\n");
    }
  }

// reachability test

rset reach(rsetof<elpair> E, rset S) {
  lset R = S;
  lset P;
  While (P!=R) {
    P = R;
    for (elpair e: E) 
      If (memberof(e.first, R))
        R += e.second;
  }
  return R;
}

void testReachable(rset A) {

  lsetof<elpair> Pairs;
     
  for(elem a: A) for(elem b: A) 
    Pairs += elpair(a,b);

  lsetof<elpair> E;
  
  for(elpair x: Pairs) 
    for (elpair y: Pairs)
      If ((x.second == y.first)
        && ((y.second != x.first)))
        E += elpair(x,y);
  
  cout << E << endl;

  for(auto a: A) for(auto b: A) If(a!=b) {

    elpair s = elpair(a,b);
    elpair t = elpair(b,a);
    
    If (memberof(t, reach(E, newSet(s))))
      cout << "Reached" << endl;
    }
  }

// minimalization test

void minimalize1(rset Q, rset F, rsetof<eltuple> delta, rset alph) {
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
    
    for(elem q1: Q) for(elem q2: Q) If(memberof(elpair(q1,q2), eq)) {
      lbool b = true;
      for(elem x: alph) for(eltuple t1: delta) for(eltuple t2: delta)
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
  
  for(elem a: Q) {
    lset t;    
    for(elem b: Q) If(memberof(elpair(a,b), eq)) 
      t += b;
    classes += t;
    }
  
  cout << "unoptimized classes: " << classes << endl;
  classes = optimize(classes);

  cout << "classes: " << classes << endl;  
  }

void minimalize2(rset Q, rset F, rsetof<eltuple> delta, rset alph) {
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
    
    for(elem EC: classes) {
      for(elem q1: EC) {
        lset EC1;
        for(elem q2: EC) {
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
void minimalize3(rset Q, rset F, rsetof<eltuple> delta, rset alph) {
  
  lsetof<elpair> E;
  for(elem p: Q) for(elem q: Q) for(elem a: alph) 
    for(elem pbis: Q) for(elem qbis: Q) 
       If(memberof(eltuple({p,a,pbis}), delta) && memberof(eltuple({q,a,qbis}), delta))
         E += elpair(elpair(pbis,qbis), elpair(p,q));

  lset S = (F * (Q&~F)) | ((Q&~F) * F);
  lset equiv = Q*Q &~ reach(E, S);
  
  lset classes;
  
  for(elem q: Q) {
    lset qclass;
    for(elem p: Q) If(memberof(elpair(p,q), equiv)) qclass += p;
    classes += qclass;
    }
  
  cout << "classes: " << classes << endl;  
  }

void mtestA(rset A, int id) {
  lset Q;
  lset F;
  lsetof<eltuple> delta;
  
  lset R;

//   for (elem x:A) for (elem y:A) {
//     elem u=elpair(x,y);
//     If (x==y)
// 		R+= elpair(x,y);
// 	If (x!=y)
//       R+=eltuple({x,x,y});
// //    R += u;
//   }
//   cout << "R= " << R << endl;
    

  rset alph = A;

  Q += 0;
  for(elem a: A) Q += eltuple({a});
  for(elem a: A) for(elem b: A) Q += eltuple({a, b});
  for(elem a: A) for(elem b: A) for(elem c: A) Q += eltuple({a, b, c});
  Q += 4;
  
  for(elem a: A) for(elem b: A) for(elem c: A) 
    If(a==b || a==c || b==c) F += eltuple({a,b,c});
    
  for(elem x: A) {
    delta += eltuple({elof(0),x,elof(eltuple({x}))});

    for(elem a: A)
      delta += eltuple({elof(eltuple({a})),x,elof(eltuple({a,x}))});

    for(elem a: A) for(elem b: A)
      delta += eltuple({elof(eltuple({a,b})),x,elof(eltuple({a,b,x}))});

    for(elem a: A) for(elem b: A) for(elem c: A)
      delta += eltuple({elof(eltuple({a,b,c})),x,elof(4)});

    delta += eltuple({elof(4),x,elof(4)});
    }

  /* Q += 0; Q += 1; Q += 2; Q += 3; F += 0; F += 2; F += 1;
  for(elem x: A) delta += eltuple({elof(0), x, elof(1)});
  for(elem x: A) delta += eltuple({elof(1), x, elof(2)});
  for(elem x: A) delta += eltuple({elof(2), x, elof(3)});
  for(elem x: A) delta += eltuple({elof(3), x, elof(0)}); */

  if(id == 1) minimalize1(Q, F, delta, alph);
  if(id == 2) minimalize2(Q, F, delta, alph);
  if(id == 3) minimalize3(Q, F, delta, alph);
  }

void mtestB(rset A, int id) {
  RelInt real(sym.greater, sym.leq, sym.max, sym.min, sym.plus, sym.times, sym.minus, sym.divide);
  pushorder po(&real);
  Domain *d = getDomain(A);

  term r0 = real.constant(d, 0);
  term r3 = real.constant(d, 3);

  lset Nat = FILTER(x, A, x >= r0);
  
  lset alph = Nat; alph += '$';
  lset Q = Nat * Nat;
  lset F = newSet(r3) * Nat;
  lsetof<eltuple> delta;
  for(elem p: Q) {
    term m = as<term>(as<elpair> (p).first);
    term n = as<term>(as<elpair> (p).second);
    delta += eltuple({p, '$', elpair(m, r0)});
    for(elem x: Nat)
      delta += eltuple({p, x, elpair(real.max(m, real.plus(n, as<term>(x))), real.plus(n, as<term>(x)))});
    }
  
  if(id == 1) minimalize1(Q, F, delta, alph);
  if(id == 2) minimalize2(Q, F, delta, alph);
  if(id == 3) minimalize3(Q, F, delta, alph);
  }

void testMinimizeA1(rset A) { mtestA(A, 1); }
                 
void testMinimizeA2(rset A) { mtestA(A, 2); }

void testMinimizeA3(rset A) { mtestA(A, 3); }

void testMinimizeB1(rset A) { mtestB(A, 1); }

void testMinimizeB2(rset A) { mtestB(A, 2); }

void testMinimizeB3(rset A) { mtestB(A, 3); }
// testReal

void testInt(rset A) { 
  string symg = "gt";
  RelInt rint(symg, sym.leq, sym.max, sym.min, sym.plus, sym.times, sym.minus, sym.divide);
  Domain *D = getDomain(A);

  RelOrder *ord = mainOrder;
  using namespace orderedfield_ops;
  mainOrder = mainField = &rint; mainDomain = getDomain(A);

  lset odd;
  for(elem x: A) odd += 1 + 2 * x;
  lset z;
  for(elem x: odd) If(5 < x) z += x + 1;
  cout << z << endl;
  If(memberof(8, z)) cout << "eight in the set" << endl;
  If(memberof(9, z)) cout << "nine in the set" << endl;

  mainOrder = ord;
  }

void testReal(rset A) { 
  string symg = "gt";
  RelReal real(symg, sym.leq, sym.max, sym.min, sym.plus, sym.times, sym.minus, sym.divide);
  pushorder po(&real);
  Domain *d = getDomain(A);

  RelOrder *ord = mainOrder;
  using namespace orderedfield_ops;
  mainOrder = mainField = &real; mainDomain = getDomain(A);

  for(elem x: A) {
    lbool b = false;

    b |= (x >= 4 && x <= elem(5) - 2);

    cout << "b = " << b << endl;

    If(b) cout << "b can be true" << endl;
    If(!b) cout << "b can be false" << endl;
    }

  mainOrder = ord;
  }

// packings on an interval

rset packings(rset F, RelReal& real, term radius) {
  lset done;
  lset active = newSet(emptyset);

  for (auto P: active) {    
    cout << "Considering " << P << " where " << emptycontext << endl;
    lset N;
        
    for (elem x: P) for(elem y: F)
      If(y > real.minus(as<term>(x), radius) && y < real.plus(as<term>(x), radius))
        N += y;
    
    Ife (F == N)
      done += P;
    else
      for (elem x: F) If(!memberof(x, N)) {
        lset newP = (asSet(P) | newSet(x));
        active += newP;
        }

    }

  return done;
  } 

void testPacking(rset A) {
  RelReal real(sym.greater, sym.leq, sym.max, sym.min, sym.plus, sym.times, sym.minus, sym.divide);
  pushorder po(&real);
  Domain *d = getDomain(A);

  lset F;

  RelOrder *ord = mainOrder;
  using namespace orderedfield_ops;
  mainOrder = mainField = &real; mainDomain = getDomain(A);

  for(elem x: A) If(x>0 && x<5) F += x;

  cout << "F = " << F << endl;

  lset pack = packings(F, real, real.constant(d, 1));

  cout << "packings = " << pack << endl;

  for(int i=0; i<9; i++)
    for(auto P: pack) {
      If(cardinality(asSet(P)) == i)
        cout << "There is a packing of size " << i << endl;
      }

  mainOrder = ord;
  }

// packings of circles

// Requires NRA, and z3 does not know whether it is possible to fit two points in
// distance at least 3 in a circle of radius 1 anyway...

typedef elpair point;

term getX(elpair p) { return as<term> (p.first); }
term getY(elpair p) { return as<term> (p.second); }

rsetof<rsetof<point>> circlepackings(rsetof<point> F, RelReal& real, term mindist) {
  lsetof<rsetof<point>> done;
  lsetof<rsetof<point>> active;
  active += emptyset;

  for (rsetof<point> P: active) {
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

void testCirclePacking(rset A) {
  RelReal real(sym.greater, sym.leq, sym.max, sym.min, sym.plus, sym.times, sym.minus, sym.divide);
  pushorder po(&real);
  Domain *d = getDomain(A);

  lsetof<point> F;

  term mindist = real.constant(d, 3);
  term bigradius = real.constant(d, 1);

  for(elem x: A) for(elem y: A) {
    term tx = as<term> (x);
    term ty = as<term> (y);
    If(real.plus(real.times(tx, tx), real.times(ty, ty)) <= real.times(bigradius, bigradius))
      F += make_pair(tx, ty);
    }

  cout << "F = " << F << endl;

  lset pack = circlepackings(F, real, mindist);

  cout << "circle packings = " << pack << endl;

  for(int i=0; i<19; i++)
    for(auto P: pack) {
      If(cardinality(asSet(P)) == i)
        cout << "There is a circle packing of size " << i << endl;
      }
  }

// tools

long long getVa() {
  struct timeval tval;
  gettimeofday(&tval, NULL);
  return tval.tv_sec * 1000000 + tval.tv_usec;
  }

typedef void loistest(rset A);

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
    lset A = dA.getSet();

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

  testSolvers("minimizeA3", testMinimizeA3);
  testSolvers("minimizeA1", testMinimizeA1);
  testSolvers("minimizeA2", testMinimizeA2);


  fclose(outtable);
  fclose(outparse);
  return 0;
  }
