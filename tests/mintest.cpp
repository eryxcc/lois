// automatically test the LOIS library.

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

/*
inline lset& operator += (lset& x, const lset& y) {
  return x += elem(y);
  } */

/*
inline lset& operator += (lset& x, const rset y) {
  return x += elem(y);
  } */

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

int main() {
  lasttime = getVa();
  initLois();
  // pushSolverDiagnostic("checking: ");

  sym.neq = "≠";

// this will not output anything
// and thus work more effectively
#define solverNamed(x,y) y

// avoid indentation to create a smaller incremental0.smt
//  iindent.i = iunindent.i = 0;

//solver = 
//  solverBasic() ||
//  solverCompare({
//      solverNamed("internal", solverExhaustive(500, false) || solverExhaustive(2000000000, true)),  
//      solverNamed("z3", solverIncremental("z3-*/bin/z3 -smt2 -in")),
      // solverNamed("crash", solverIncremental("echo crash")),
        // solverNamed("cvc4", solverIncremental("cvc4 --lang smt --incremental --finite-model-find")),
     // solverNamed("cvc4", solverIncremental("cvc4 --lang smt --incremental"))
//      }) ||
//  solverCrash();

  Domain dA("Real");
  rset A = dA.getSet();
  
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

  // minimalize1(Q, F, delta, alph);
  showTimeElapsed();

  // minimalize2(Q, F, delta, alph);
  // showTimeElapsed();

  return 0;
  }
