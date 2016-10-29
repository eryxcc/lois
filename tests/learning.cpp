// Learning nominal automata (INCOMPLETE!)
// Based on the paper by Joshua Moerman, Matteo Sammartino, Alexandra Silva, 
// Bartek Klin, Micha Szynwelski

// See https://arxiv.org/pdf/1607.06268.pdf

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

typedef vector<elem> word;

struct transition {
  elem src;
  elem symbol;
  elem tgt;
  transition(elem _src, elem _symbol, elem _tgt) { src=_src; symbol=_symbol; tgt=_tgt; }
  };

inline std::ostream& operator << (std::ostream& os, transition t) { 
  return os << "(" << t.src << "," << t.symbol << "," << t.tgt << ")"; }

inline bool isused(vptr v, transition p) { return isused(v, p.src) || isused(v, p.symbol) || isused(v, p.tgt); }
inline transition alpha(transition p, vptr v1, vptr v2) { 
  return transition(alpha(p.src, v1, v2), alpha(p.symbol, v1, v2), alpha(p.tgt, v1, v2));
  }

inline rbool operator == (transition a, transition b) { 
  return a.src == b.src && a.symbol == b.symbol && a.tgt == b.tgt; }
inline rbool operator != (transition a, transition b) { 
  return a.src != b.src || a.symbol != b.symbol && a.tgt == b.tgt; }
inline rbool operator <  (transition a, transition b) { 
  lbool ret;
  Ife(a.src != b.src)
    ret = a.src < b.src;
  else Ife(a.symbol != b.symbol)
    ret = a.symbol < b.symbol;
  else ret = a.tgt < b.tgt;
  }
inline rbool operator <= (transition a, transition b) { 
  lbool ret;
  Ife(a.src != b.src)
    ret = a.src < b.src;
  else Ife(a.symbol != b.symbol)
    ret = a.symbol < b.symbol;
  else ret = a.tgt <= b.tgt;
  }

struct dfa {
  lset Q;
  lsetof<transition> delta;
  lset F;
  lelem I;
  };

lset sigma;

lbool wordinlanguage_aux(word w, const dfa& a, int pos, lelem state) {
  if(pos == w.size()) return memberof(state, a.F);
  elem c = w[pos];
  
  lbool yes = false;

  for(transition& t: a.delta) 
    If(t.src == state && t.symbol == w[pos]) 
      yes |= wordinlanguage_aux(w, a, pos+1, t.tgt);
    
  return yes;
  }

lbool wordinlanguage(word w, const dfa& a) {
  return wordinlanguage_aux(w, a, 0, a.I);
  }

word concat(word x, word y) {
  word res;
  for(elem a: x) res.push_back(a);
  for(elem b: y) res.push_back(b);
  return res;
  }

word concat(elem a, word y) {
  word res;
  res.push_back(a);
  for(elem b: y) res.push_back(b);
  return res;
  }

word concat(word x, elem b) {
  word res;
  for(elem a: x) res.push_back(a);
  res.push_back(b);
  return res;
  }

word concat(word x, word y, word z) {
  word res;
  for(elem a: x) res.push_back(a);
  for(elem b: y) res.push_back(b);
  for(elem c: z) res.push_back(c);
  return res;
  }

word concat(word x, elem b, word z) {
  word res;
  for(elem a: x) res.push_back(a);
  res.push_back(b);
  for(elem c: z) res.push_back(c);
  return res;
  }

lbool operator ^ (lbool x, lbool y) {
  return (x&&!y) || (y&&!x);
  }

void printDFA(const dfa& L) {
  std::cout << "Q = " << L.Q << std::endl;
  std::cout << "I = " << L.I << std::endl;
  std::cout << "F = " << L.F << std::endl;
  std::cout << "δ = " << L.delta << std::endl;
  }
  
void learning(const dfa& L) {

  std::cout << "DFA to learn:" << std::endl;
  printDFA(L);
  std::cout << std::endl;

  lsetof<word> S, E;
  S += word();
  E += word();
  
  again: ;
  
  lbool changed = true;

  While(changed) {

    std::cout << "Observation table:" << std::endl;
    std::cout << "S = " << S << std::endl;
    std::cout << "E = " << E << std::endl;
    std::cout << std::endl;  

    changed = false;

    for(auto s: S) for(auto a: sigma) {
      lbool ok2 = false;
      for(auto t: S) {
        lbool ok = true;
        for(auto e: E) 
          If(wordinlanguage(concat((s),a,(e)), L) ^ wordinlanguage(concat((t),(e)), L))
            ok = false;
        ok2 |= ok;
        }
      If(!ok2) {  
        std::cout << 
          "Not closed. Adding to S: " << concat(s,a) << " for " << emptycontext
          << std::endl << std::endl;
        S += concat((s),a); changed = true; 
        }
      }
    
    for(auto s: S) for(auto t: S) {
      lbool ok = true;
      for(auto e: E) If(wordinlanguage(concat((s),e), L) ^ wordinlanguage(concat((t),e), L))
        ok = false;
      If(ok) for(auto a: sigma) for(auto e: E)
        If(wordinlanguage(concat((s),a,(e)), L) ^ wordinlanguage(concat((t),a,(e)), L))
          If(!memberof(concat(a,e), E)) {
            std::cout << 
              "Not consistent. Adding to E: " << concat(a,e) << " for " << emptycontext
              << std::endl << std::endl;
            E += concat(a,e);
            changed = true;
            }
      }
    
    }
  
  dfa Learned;

  for(auto s: S) {
    lset state;
    for(auto e: E) 
      If(wordinlanguage(concat(s,e), L))
        state += e;
      
    If(!memberof(state, Learned.Q)) {
      Learned.Q += state;
      If(s == word()) Learned.I = state;
      If(wordinlanguage(s, L))
        Learned.F += state;

      for(auto a: sigma) {
        lset q2;
        for(auto e: E) If(wordinlanguage(concat(s,a,e), L))
          q2 += e;
        Learned.delta += transition(state, a, q2);
        }
      }
    }
  
  std::cout << "Guessed DFA:" << std::endl;
  printDFA(Learned);
  std::cout << std::endl;
  
  // (q1,q2) \in compare iff, after reading some word, the
  // DFA L is in state q1 and the DFA Learned
  // is in state q2
  lsetof<elpair> compare;
  
  // for each (q1,q2) in compare, witnesses contains ((q1,q2), w),
  // where w is the witness word
  lsetof<elpair> witnesses;

  for(elem e: newSet(L.I)) for(elem e2: newSet(Learned.I)) {
    auto initpair = elpair(e, e2);
    compare += initpair;
    witnesses += elpair(initpair, word());
    }
  
  for(auto witness: witnesses) {

    elem q1 = as<elpair> (witness.first).first;
    elem q2 = as<elpair> (witness.first).second;
    
    word w = as<word> (witness.second);
    
    lbool m1 = memberof(q1, L.F);
    lbool m2 = memberof(q2, Learned.F);

    If(m1 ^ m2) {
      std::cout << "Counterexample: " << w << " where " << emptycontext << std::endl << std::endl;
      while(true) {
        If(!memberof(w, S)) S += w;
        if(w.size() == 0) goto again;
        w.resize(w.size() - 1);
        }
      }
    
    for(auto a: sigma)
      for(transition& t1: L.delta) 
        If(t1.src == q1 && t1.symbol == a) 
      for(transition& t2: Learned.delta) 
        If(t2.src == q2 && t2.symbol == a) {
          elpair p2 = elpair(t1.tgt, t2.tgt);
          If(!memberof(p2, compare)) {
            compare += p2;
            witnesses += make_pair(p2, concat(w, a));
            }
          }
    }

  std::cout << "Learning successful!" << std::endl;
  }
  
int main() {
  lasttime = getVa();
  initLois();
  // pushSolverDiagnostic("checking: ");

  Domain dA("Atoms");
  lset A = dA.getSet();

  sigma = A;

  sym.neq = "≠";
  
  dfa target;
  
  // language L1 from the paper (repeated letter)
  
/*
  elem e0 = 0;
  elem e1 = 1;
  elem e2 = 2;
  target.Q += e0;
  target.Q += e1;
  target.Q += e2;
  for(auto a:A) target.Q += a;
  target.F += e1;
  target.I = e0;
  for(auto a:A) target.delta += transition(e0, a, a);
  for(auto a:A) target.delta += transition(a, a, e1);
  for(auto a:A) for(auto b: A) If(a != b) target.delta += transition(a, b, e2);
  for(auto a:A) target.delta += transition(e2, a, e2);
*/

  elem eini = 0;   // initial
  elem etrash = 1; // trash
  elem eaccept = 2; // accept
  target.Q += eini;
  target.Q += etrash;
  target.Q += eaccept;
  for(auto a: A) target.Q += a; // read one letter
  for(auto a: A) for(auto b: A) target.Q += elpair(a,b); // read two letters
  for(auto a: A) target.Q += elpair(a,eini); // read three letters, wait for 'a'
  target.F += eaccept;
  target.I = eini;
  
  // language L2 from the paper (repeated two letters: 'baba')

  for(auto a: A) target.delta += transition(eini, a, a);
  for(auto a: A) for(auto b: A)
    target.delta += transition(a, b, elpair(a,b));
  for(auto a: A) for(auto b: A)
    target.delta += transition(elpair(a,b), a, elpair(b,eini));
  for(auto a: A) for(auto b: A) for(auto c: A) If(a != c)
    target.delta += transition(elpair(a,b), c, etrash);
  for(auto a: A) 
    target.delta += transition(elpair(a,eini), a, eaccept);
  for(auto a: A) for(auto b: A) If(a != b)
    target.delta += transition(elpair(a,eini), b, etrash);
  for(auto a: A) target.delta += transition(etrash, a, etrash);
  for(auto a: A) target.delta += transition(eaccept, a, etrash);
  
  
  lset allwords;
  
  for(auto a:A) for(auto b: A) {
    word w;
    w.push_back(a);
    w.push_back(b);
    If(wordinlanguage(w, target)) allwords += w;
    }

  std::cout << "All words of length 2: " << allwords << std::endl << std::endl;
  
  allwords = newSet();
  
  for(auto a:A) for(auto b: A) for(auto c: A) for(auto d: A) {
    word w;
    w.push_back(a);
    w.push_back(b);
    w.push_back(c);
    w.push_back(d);
    If(wordinlanguage(w, target)) allwords += w;
    }

  std::cout << "All words of length 4: " << allwords << std::endl << std::endl;
  
  learning(target);
  
  showTimeElapsed();
  
  /* for(elem a: A) {
    lelem x;
    for(elem b: A) If(a==b) x = b;
    } */

  return 0;
  }
