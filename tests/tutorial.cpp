// Note that this example is parsed by
// parsetutorial.cpp (not included in the package) in order to include
// its fragments and output into the paper easily.

// This file has the following parts:
// - the "Tutorial" from the paper.
// - more examples from the paper. Output is also included in some cases.
// - example LOIS functions from the paper (no output included).

#define SECTION(x) cout << "SECTION " << x << endl;
#define FUNCTION(x) 


FUNCTION("initialization")

#include "../include/loisextra.h"

using namespace std;
using namespace lois;


FUNCTION("min")
lsetof<term> min(lsetof<term> X) { 
  lsetof<term> answer;
  for(term& x: X) If (FORALL(y, X, x <= y)) answer += x;
  return answer;
}

FUNCTION("interval")
lsetof<term> interval(term a, 
  term b, lsetof<term> A) {  
  return FILTER(x, A, 
    (a < x) && (x < b), term);
}
  
FUNCTION("sup")
lsetof<term> supremum(lsetof<term> X, lsetof<term> domain) { 
  return min(FILTER(m, domain, FORALL(x, X, m >= x), term));
}

FUNCTION("reachability")
template<class T> lsetof<T> reach 
  (lsetof<lpair<T,T>> E, 
     lsetof<T> S) {
  lsetof<T> R = S;
  lsetof<T> P;
  While (P!=R) {
    P = R;
    for (auto& e: E) 
      If (memberof(e.first, R))
        R += e.second;
  }
  return R;
}
FUNCTION("automaton")
template<class T, class U> struct automaton {
  lsetof<T> states;
  lsetof<lpair<T, lpair<U, T>>> transitions;
  lsetof<T> initial;
  lsetof<T> final;
  lsetof<U> alphabet;  
};
FUNCTION("emptiness")
template<class T, class U> rbool isEmpty(automaton<T,U> A) {
  lsetof<lpair<T,T>> E;
  for (auto& t: A.transitions)
    E += make_lpair(t.first,t.second.second);

  return (A.final&reach(E, A.initial))!=lsetof<T>();
}


FUNCTION("main_old")
int main() {
SECTION("main")
initLois();
sym.useLaTeX();
//sym.useUnicode();
SECTION("domain")
Domain dA("\\mathbb{A}");  
lsetof<term> A = dA.getSet();
cout << "A = ";
cout << A;
cout << endl;
SECTION("resetX")
lsetof<term> X = A;
// X += 5;

cout << "X = " << X << endl;
SECTION("pairs")
lsetof<lpair<term,term>> Pairs;
     
for(term a: A) for(term b: A) 
    Pairs += make_lpair(a,b);

cout << "Pairs = " << Pairs << endl;
SECTION("linear")  
lbool ok = true;

for (term a:A) for (term b:A) for (term c:A)
  If (((a<=b) && (b<=c) 
    && !(a<=c))) ok = false;

If (ok) 
  cout << "Transitive" << endl;
SECTION("more")
// --- more snippets follow ---
  
SECTION("newset")
  // X = lsetof<term>();
  cout << "X = " << X << endl;

SECTION("const")
  /* for(int i=1; i<=3; i++) 
    X += i; */
  // cout << "X = " << X << endl;

SECTION("loop")
  lsetof<lpair<term,term>> Y;
  for(term a: X) for(term b: X) If (a != b)
    Y += make_lpair(a, b);
  cout << "Y = " << Y << endl;

// SECTION("declareatom")
//
//   declareatom a1(&dA, "a_1"), a2(&dA, "a_2");
//   axiom ax1(a1 != a2);
//   lset A1, A2;
//   for(elem a: A) If (a != a1) A1 += a;
//   for(elem a: A) If (a != a2) A2 += a;
//   cout << (A1 == A2) << endl;

SECTION("control")

  lsetof<lsetof<lsetof<term>>> C; 
  for(term a: A) {
    lsetof<lsetof<term>> D;
    for(term b: A) {
      If (a != b) {
        D += newSet(b);
        }
      }
    cout << "D= " << D << emptycontext << endl;
    C += D;
    }
  cout << "C=" << C << endl;

// SECTION("pseudonum")
//
//   lset E;
//   E |= newSet(a1, a2);
//
//   lnum<int> x = 0;
//   for(elem e: optimize(E)) x++;
//   cout << "|E| = " << x << endl;

SECTION("to_many_Zs")
  {
SECTION("whypiecewise")

/*
  lset Z;
  for (elem x:A) for (elem y:A) {
    lelem u = x;
    If (x != y)
      u = newSet(elpair(x,y));
    Z += u; 
  } */

SECTION("to_many_Zs_end")
  }

{
SECTION("semantics1")
lsetof<term> X;
for (term x:A) 
  X += x;

SECTION("tech1")
} {

SECTION("semantics2")
lsetof<lsetof<term>> X;
for (term x:A) {
  lsetof<term> Y;
  for (term y: A) 
    If(x!=y)
      Y += y;
  X += Y;
  }

SECTION("tech2")
}

/*
SECTION("assignmentex")
lsetof<term> Z;
for (term x:A) for (term y:A) 
{
  lsetof<term> X = newSet(x);
  If (x != y)
    X = newSet(make_lpair(x,y));
  Z += X;
}
*/
SECTION("simpleloop")
  for (auto a:A) 
    auto X = newSet(a);

SECTION("j1")
  { // too many Y's

SECTION("setup")
/*
X = newSet<term>();
contextptr C = currentcontext;
for (elem a:A) for (elem b:A) {
  for (elem c:A) 
    If (a==b) X+=eltuple({a,b,c});
  If (a!=b) {
    cout << "X=" << X << endl;
    cout << "C=" << C << endl;
    }
  }
    */

SECTION("simpleloop2")
lsetof<term> Y;

for (term x:X) Y += x;

If (X == Y) cout << "equal" << endl;
If (X != Y) cout << "not equal" << endl;
SECTION("test_equal")  
  If (X == Y) cout << "equal" << endl;
  If (X != Y) cout << "not equal" << endl;


SECTION("j2")
  }
 
// SECTION{"petri"}
//   bool coverable(lset P, lset V, lelem I, lelem F)
//     //test coverability of a target vector F from an initial vector I
//       //in a Petri net P -- a finite set of vectors in Z^k
//     //V is the space N^k
//   {
//       lset R;
//     lset P;
//
//     for (elem v: V)
//       If (as<eltuple> (v) >= F)
//         R += v;
//
//     while (P != R) {
//       P = R;
//       for (elem v: V) {
//         If (EXISTS(w, A, ismember(eltuple(v,w,  > y))
//
//       }
//     }
//
//
//
//
//
//   }

 
SECTION("treetest")
  RelTree tree(sym.arrow, sym.notarrow, sym.min);

  for(term a: A) for(term b : A) If (tree.anceq(a, b)) {
    lsetof<term> I;
    for(term c: A) If (tree.anceq(a,c) && tree.anceq(c,b)) 
      I += c;
    
    lbool ok = true;
   
    for(term x: I) for(term y: I) for(term z: I) {
      If (!(tree.anceq(x,y) || tree.anceq(y,x))) 
        ok = false;
      If (tree.anceq(x,y) && tree.anceq(y,z) && !tree.anceq(x,z)) 
        ok = false;
      }

    Ife(ok) cout << "total order" << endl;
    else cout << "not a total order" << endl;
    }
             

  

SECTION("count")
  int i=0;
  for (auto a:Pairs) i++;
  cout << i << endl;

             
SECTION("set2")
  lsetof<lsetof<term>> UnorderedPairs;
       
  for(term a: A) for(term b: A) 
    If (a!=b)
      UnorderedPairs += newSet(a,b);

  cout << UnorderedPairs << endl;

SECTION("set3x")
  lsetof<lsetof<term>> coPairs;

  for(auto x: UnorderedPairs)
    coPairs += A &~ x;
  cout << coPairs << endl;
    
SECTION("set3")
lsetof<lsetof<term>> Intvs;
   
for(auto& p: Pairs) 
  If (p.first<p.second) 
    Intvs += interval(p.first, p.second, A);

  cout << "Intvs = " << Intvs << endl;
SECTION("sup_test")
  lbool sup = true;

  for(auto i: Intvs) 
    If (supremum(i,A) == lsetof<term>()) sup = false;
  
  If (sup) cout << "Yes" << endl;

SECTION("graph")
lsetof<lpair<lpair<term,term>, lpair<term,term>>> E;

for(auto x: Pairs) 
  for (auto y: Pairs)
    If ((x.second == y.first)
      && ((y.second != x.first)))
      E += make_lpair(x,y);

cout << E << endl;
SECTION("declareatom")
declareatom u(&dA, "u"),
  v(&dA, "v");
axiom ax(u != v);
  
SECTION("test_reach")
lbool reached = false;

for (term a: A) 
  for (term b:A) {
    auto s = make_lpair(a,b);
    auto t = make_lpair(b,a);
    If (memberof(t,reach(E, newSet(s))))
      reached = true;
  }

If (reached) cout << "Reached" << endl;
  
SECTION("displaycontext")

  contextptr c = currentcontext;

  for(auto a: A) for(auto b: A) {

    a.asVar()->name = "a";
    b.asVar()->name = "b";

    If (a<b)
      cout << c << endl;
    }

SECTION("lnumtest")

  for(auto a: A) for(auto b: A) {

    a.asVar()->name = "a";
    b.asVar()->name = "b";

    lelemof<int> i = 0;

    for(auto x: A) 
      If (x == a || x == b)
        i++;

    cout << i << endl;
    }

SECTION("quantifiers")

  cout << FORALL(x, A, EXISTS(y, A, x > y)) << endl;

SECTION("mapfilter")

  auto mf =MAP(x, A, FILTER(y, A, x>y, term), lsetof<term>);

  cout << mf << endl;

SECTION("finish")
  return 0;
  }

//--- more examples from the paper: functions


FUNCTION("cardinalitybad")

template<class T> int cardinalityBad(lsetof<T>& X) {
  int result = 0;
  for(auto x: X) result++;
  return result;
}

/*
FUNCTION("memberof")

rbool ismemberof(elem a, lset B) {
  lbool phi(ffalse);
  for(elem b: B) 
    If (a == b) 
      phi |= ftrue;
  return phi;
}

  rbool issubseteq(lset A, lset B) {
    lbool phi(ftrue);
    for(elem a: A) 
      If (!memberof(a, B)) 
        phi &= ffalse;
    return phi;
  }

FUNCTION("card")

  lsetof<int> card(lset X) {
    lsetof<int> answer;
    for(elem a: X) 
      If (answer.isEmpty()) {
        lsetof<int> butone = card(X &~ newSet(a));
        for(int x: butone) answer += (x+1);
        }
    If (answer.isEmpty()) answer += 0;
    return answer;
  }

FUNCTION("max")

  lelem max(lset X) { 
    lset answer;
    for(elem x: X) If (FORALL(y, X, x >= y)) answer += x;
    return extract(answer);
  }

FUNCTION("maxbad")

  lelem maxBad(lset X) {
    lelem ans;
    for(elem x: X) 
      If (ans.isUndefined() || x > ans)
        ans = x;
    return ans;
  }

  */