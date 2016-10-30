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
lset min(lset X) { 
  lset answer;
  for(elem x: X) If (FORALL(y, X, x <= y)) answer += x;
  return answer;
}

FUNCTION("interval")
lset interval(term a, 
  term b, lset A) {  
  return FILTER(x, A, 
    (a < x) && (x < b));
}
  
FUNCTION("sup")
lset supremum(lset X, lset domain) { 
  return min(FILTER(m, domain, FORALL(x, X, m >= x)));
}

FUNCTION("reachability")
lset reach 
  (lsetof<elpair> E, 
     lset S) {
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
FUNCTION("automaton")
struct automaton {
  lset states;
  lsetof<eltuple> transitions;
  lset initial;
  lset final;
  lset alphabet;  
};
FUNCTION("emptiness")
rbool isEmpty(automaton A) {
  lsetof<elpair> E;
  for (eltuple t: A.transitions)
    E += elpair(t[0],t[2]);

  return (A.final&reach(E, A.initial))!=newSet();
}


FUNCTION("main_old")
int main() {
SECTION("main")
initLois();
sym.useLaTeX();
SECTION("domain")
Domain dA("\\mathbb{A}");  
lset A = dA.getSet();
cout << "A = " << A << endl;
SECTION("resetX")
lset X = A;
X += 5;

cout << "X = " << X << endl;
SECTION("pairs")
lsetof<elpair> Pairs;
     
for(elem a: A) for(elem b: A) 
    Pairs += elpair(a,b);

cout << "Pairs = " 
  << Pairs << endl;
SECTION("linear")  
lbool ok = true;

for (elem a:A) for (elem b:A) for (elem c:A)
  If (((a<=b) && (b<=c) 
    && !(a<=c))) ok = false;

If (ok) 
  cout << "Transitive" << endl;
SECTION("more")
// --- more snippets follow ---
  
SECTION("newset")
  X = newSet();
  cout << "X = " << X << endl;

SECTION("const")
  for(int i=1; i<=3; i++) 
    X += i;
  cout << "X = " << X << endl;

SECTION("loop")
  lset Y;
  for(elem a: X) for(elem b: X) If (a != b)
    Y += elpair(a, b);
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

  lset C; 
  for(elem a: A) {
    lset D;
    for(elem b: A) 
      If (a != b) 
        D += newSet(b);
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

  lset Z;
  for (elem x:A) for (elem y:A) {
    lelem u = x;
    If (x != y)
      u = newSet(elpair(x,y));
    Z += u;
  }

SECTION("to_many_Zs_end")
  }

{
SECTION("semantics1")
lset X;
for (elem x:A) 
  X += x;

SECTION("tech1")
} {

SECTION("semantics2")
lset X;
for (elem x:A) {
  lset Y;
  for (elem y: A) 
    If(x!=y)
      Y += y;
  X += Y;
  }

SECTION("tech2")
}

SECTION("assignmentex")
lset Z;
for (elem x:A) for (elem y:A) 
{
  lset X = newSet(x);
  If (x != y)
    X = newSet(elpair(x,y));
  Z += X;
}
SECTION("simpleloop")
  for (elem a:A) 
    lset X = newSet(a);

SECTION("j1")
  { // too many Y's

SECTION("setup")
X = newSet();
contextptr C = currentcontext;
for (elem a:A) for (elem b:A) {
  for (elem c:A) 
    If (a==b) X+=eltuple({a,b,c});
  If (a!=b) {
    cout << "X=" << X << endl;
    cout << "C=" << C << endl;
SECTION("simpleloop2")
lset Y;

for (elem x:X) Y += x;

If (X == Y) cout << "equal" << endl;
If (X != Y) cout << "not equal" << endl;
SECTION("test_equal")  
  If (X == Y) cout << "equal" << endl;
  If (X != Y) cout << "not equal" << endl;
}
}

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

  for(elem a: A) for(elem b : A) If (tree.anceq(a, b)) {
    lset I;
    for(elem c: A) If (tree.anceq(a,c) && tree.anceq(c,b)) 
      I += c;
    
    lbool ok = true;

    for(elem x: I) for(elem y: I) for(elem z: I) {
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
  for (elem a:Pairs) i++;
  cout << i << endl;

         
SECTION("set2")
  lset UnorderedPairs;
       
  for(elem a: A) for(elem b: A) 
    If (a!=b)
      UnorderedPairs += newSet(a,b);

SECTION("set3x")
  lset coPairs;

  for(elem x: UnorderedPairs)
    coPairs += A &~ asSet(x);
  cout << coPairs << endl;

SECTION("set3")
lset Intvs;
   
for(elpair p: Pairs) 
  If (p.first<p.second) 
    Intvs += interval(
       as<term>(p.first),
       as<term>(p.second),A);

cout << "Intvs = " 
  << Intvs << endl;  
SECTION("sup_test")
  lbool sup = true;

  for(elem i: Intvs) 
    If (supremum(asSet(i),A) == newSet()) sup = false;
  
  If (sup) cout << "Yes" << endl;

SECTION("graph")
lsetof<elpair> E;

for(elpair x: Pairs) 
  for (elpair y: Pairs)
    If ((x.second == y.first)
      && ((y.second != x.first)))
      E += elpair(x,y);

cout << E << endl;
SECTION("declareatom")
declareatom u(&dA, "u"),
  v(&dA, "v");
axiom ax(u != v);
  
SECTION("test_reach")
lbool reached = false;

for (elem a: A) 
  for (elem b:A) {
    elpair s = elpair(a,b);
    elpair t = elpair(b,a);
    If (memberof(t,reach(E, newSet(s))))
      reached = true;
}

If (reached) cout << "Reached" << endl;
  
SECTION("displaycontext")

  contextptr c = currentcontext;

  for(auto a: A) for(auto b: A) {

    as<term>(a).asVar()->name = "a";
    as<term>(b).asVar()->name = "b";

    If (a<b)
      cout << c << endl;
    }

SECTION("lnumtest")

  for(auto a: A) for(auto b: A) {

    as<term>(a).asVar()->name = "a";
    as<term>(b).asVar()->name = "b";

    lnum<int> i = 0;

    for(auto x: A) 
      If (x == a || x == b)
        i++;

    cout << i << endl;
    }

SECTION("quantifiers")

  cout << FORALL(x, A, EXISTS(y, A, x > y)) << endl;

SECTION("mapfilter")

  cout << MAP(x, A, FILTER(y, A, x>y)) << endl;

SECTION("finish")
  return 0;
  }

//--- more examples from the paper: functions


FUNCTION("cardinalitybad")

int cardinalityBad(elem X) {
  int result = 0;
  for(elem x: X) result++;
  return result;
}

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

