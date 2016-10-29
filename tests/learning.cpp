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

struct automaton {
  lset Q;
  lsetof<transition> delta;
  lset F;
  lelem I;
  };

lset sigma;

lbool wordinlanguage_aux(word w, const automaton& a, int pos, lelem state) {
  if(pos == w.size()) return memberof(state, a.F);
  elem c = w[pos];
  
  lbool yes = false;

  for(transition& t: a.delta) 
    If(t.src == state && t.symbol == w[pos]) 
      yes |= wordinlanguage_aux(w, a, pos+1, t.tgt);
    
  return yes;
  }

lbool wordinlanguage(word w, const automaton& a) {
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

void printAutomaton(const automaton& L) {
  std::cout << "Q = " << L.Q << std::endl;
  std::cout << "I = " << L.I << std::endl;
  std::cout << "F = " << L.F << std::endl;
  std::cout << "δ = " << L.delta << std::endl;
  }
  
void learning(const automaton& L) {

  std::cout << "Automaton to learn:" << std::endl;
  printAutomaton(L);
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
  
  automaton Learned;

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
  
  std::cout << "Guessed automaton:" << std::endl;
  printAutomaton(Learned);
  std::cout << std::endl;
  
  lsetof<elpair> compare;
  lsetof<elpair> samples;

  // todo in LOIS : newSet for <lelem, lelem>
  for(elem e: newSet(L.I)) for(elem e2: newSet(Learned.I)) {
    auto initpair = elpair(e, e2);
    compare += initpair;
    samples += elpair(initpair, word());
    }
  
  for(auto p: compare) {
    
    for(auto samp: samples) If(samp.first == p) {

      word w = as<word> (samp.second);
      
      lbool m1 = memberof(p.first, L.F);
      lbool m2 = memberof(p.second, Learned.F);
  
      If((m1 && !m2) || (m2 && !m1)) {
        std::cout << "Counterexample: " << w << " where " << emptycontext << std::endl << std::endl;
        while(true) {
          If(!memberof(w, S)) S += w;
          if(w.size() == 0) goto again;
          w.resize(w.size() - 1);
          }
        }
      
      for(auto a: sigma)
        for(transition& t1: L.delta) 
          If(t1.src == p.first && t1.symbol == a) 
        for(transition& t2: Learned.delta) 
          If(t2.src == p.second && t2.symbol == a) {
            elpair p2 = elpair(t1.tgt, t2.tgt);
            If(!memberof(p2, compare)) {
              compare += p2;
              samples += make_pair(p2, concat(as<word> (samp.second), a));
              }
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
  sigma = dA.getSet();
  lset A = sigma;

  sym.neq = "≠";
  
  automaton target;
  elem e0 = 0;
  elem e1 = 1;
  elem e2 = 2;
  target.Q += e0;
  target.Q += e1;
  target.Q += e2;
  for(auto a:A) target.Q += a;
  target.F += e1;
  // for(auto a: A) target.F += a;
  target.I = e0;
  for(auto a:A) target.delta += transition(e0, a, a);
  for(auto a:A) target.delta += transition(a, a, e1);
  for(auto a:A) for(auto b: A) If(a != b) target.delta += transition(a, b, e2);
  for(auto a:A) target.delta += transition(e2, a, e2);
  
  lset allwords;
  
  for(auto a:A) for(auto b: A) {
    word w;
    w.push_back(a);
    w.push_back(b);
    If(wordinlanguage(w, target)) allwords += w;
    }

  std::cout << "All words of length 2: " << allwords << std::endl;
  
  learning(target);
  
  showTimeElapsed();
  
  /* for(elem a: A) {
    lelem x;
    for(elem b: A) If(a==b) x = b;
    } */

  return 0;
  }
