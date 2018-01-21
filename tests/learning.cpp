// Learning nominal automata (INCOMPLETE!)
// Based on the paper by Joshua Moerman, Matteo Sammartino, Alexandra Silva, 
// Bartek Klin, Micha Szynwelski

// flags

// use a tailored algorithm for membership checking, instead of using the automaton
// (works for the queue)
#define MEMBEROPT
               
// iterate membership instead of recursing
// #define MEMBER_ITERATIVE

// should we add counterexamples to E, rather than to S?
#define ADD_COLUMNS

// should we learn DFA or NFA?
// #define LEARN_NFA

// consider only words of length <= LENGTH_LIMIT when verifying
// (sometimes required for NFA languages)
#define LENGTH_LIMIT 4

// test on the queue language
// #define QUEUESIZE 3

// test on doubleword (currently allowing only 1 = [aa], 2 = [abab])
#define DOUBLEWORD 2

// test on stack language
// #define STACKSIZE 2

// test on nth last letter
// the word is in language iff the first letter equals (NTHLAST-1) last one
// #define NTHLAST 2

// find equal letters
// #define EQLANG

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

#include "../include/lois-weak.h"

typedef elem Symbol;
typedef elem State;

typedef lvector<Symbol> word;

ostream& operator << (ostream& os, const word& w) {
  for(auto s: w) os << s;
  return os;
  }

#include "../include/lois-automaton.h"
typedef automaton<State, Symbol> myautomaton;
typedef transition<State, Symbol> mytransition;

lbool wordinlanguage_aux(const word& w, const myautomaton& a, int pos, const lelemof<State>& state) {
//  std::cout << "in state = " << state << " in " << emptycontext << std::endl;
  if(pos == w.size()) return memberof(state, a.F);
  elem c = w[pos];
  
  lbool yes = false;

  for(auto& t: a.delta) 
    If(t.src == state && t.symbol == w[pos]) 
      yes |= wordinlanguage_aux(w, a, pos+1, t.tgt);
    
  return yes;
  }

lbool wordinlanguage(word w, const myautomaton& a) {
#ifdef MEMBEROPT
#ifdef QUEUESIZE
  word queu;
  lbool res = true;
  int qid = 0;
  for(elem c: w) {
    auto p = as<lpair<elem, term> > (c);
    Ife(p.first == elof(1)) { 
      if(queu.size() == qid + QUEUESIZE) return false;
      queu.push_back(elof(p.second));
      }
    else {
      if(qid == queu.size()) return false;
      res &= (queu[qid] == elof(p.second));
      qid++;
      }
    }
  return res;
#endif
#endif
//  std::cout << "Checking word: " << w  << " in " << emptycontext << std::endl;

#ifdef MEMBER_ITERATIVE
  lset state = a.I;
  lset nextstate;
  for(elem c: w) {
    nextstate = newSet();
    for(transition& t: a.delta) 
      If(memberof(t.src, state) && t.symbol == c) 
        nextstate += t.tgt;
    state = nextstate;
    }
  
  return (state & a.F) != emptyset;
#endif

  lbool ret = false;
  for(auto it: a.I) ret |= wordinlanguage_aux(w, a, 0, it);
  return ret;
  }

word concat(word x, word y) {
  word res;
  for(Symbol a: x) res.push_back(a);
  for(Symbol b: y) res.push_back(b);
  return res;
  }

word concat(Symbol a, word y) {
  word res;
  res.push_back(a);
  for(Symbol b: y) res.push_back(b);
  return res;
  }

word concat(word x, Symbol b) {
  word res;
  for(Symbol a: x) res.push_back(a);
  res.push_back(b);
  return res;
  }

word concat(word x, word y, word z) {
  word res;
  for(Symbol a: x) res.push_back(a);
  for(Symbol b: y) res.push_back(b);
  for(Symbol c: z) res.push_back(c);
  return res;
  }

word concat(word x, Symbol b, word z) {
  word res;
  for(Symbol a: x) res.push_back(a);
  res.push_back(b);
  for(Symbol c: z) res.push_back(c);
  return res;
  }

lbool operator ^ (lbool x, lbool y) {
  return (x&&!y) || (y&&!x);
  }

void printAutomaton(const myautomaton& L) {
  std::cout << "Q = " << L.Q << std::endl;
  std::cout << "I = " << L.I << std::endl;
  std::cout << "F = " << L.F << std::endl;
  std::cout << "δ = " << L.delta << std::endl;
  }

#ifndef LEARN_NFA
void learning(const myautomaton& L) {

  std::cout << "DFA to learn:" << std::endl;
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
    
    std::cout << "Checking closedness..." << std::endl;    

    for(auto s: S) for(auto a: L.alph) {
      lbool ok2 = false;
      for(auto t: S) {
        lbool ok = true;
        for(auto e: E) {
//        std::cout << "s="<<s<<" t="<<t<<" a="<<a<<" e="<<e << std::endl;
          If(wordinlanguage(concat((s),a,(e)), L) ^ wordinlanguage(concat((t),(e)), L))
            ok = false;
          }
        ok2 |= ok;
        }
      If(!ok2) {  
        std::cout << 
          "Not closed. Adding to S: " << concat(s,a) << " for " << emptycontext
          << std::endl << std::endl;
        S += concat((s),a); changed = true; 
        }
      }
    
    std::cout << "Checking consistency..." << std::endl;    
    
    for(auto s: S) for(auto t: S) If(s != t) {
      lbool ok = true;
      for(auto e: E) If(wordinlanguage(concat((s),e), L) ^ wordinlanguage(concat((t),e), L))
        ok = false;
      If(ok) for(auto a: L.alph) for(auto e: E)
        If(wordinlanguage(concat((s),a,(e)), L) ^ wordinlanguage(concat((t),a,(e)), L))
          /*If(!memberof(concat(a,e), E))*/ {
            std::cout << 
              "Not consistent. Adding to E: " << concat(a,e) << " for " << emptycontext
              << std::endl << std::endl;
            E += concat(a,e);
            changed = true;
            }
      }
    
    }
  
  myautomaton Learned;

  for(auto s: S) {
    lsetof<word> state;
    for(auto e: E) 
      If(wordinlanguage(concat(s,e), L))
        state += e;
      
    elem estate = elof(state);
      
    If(!memberof(estate, Learned.Q)) {
      Learned.Q += estate;
      If(s == word()) Learned.I = newSet<State> (estate);
      If(wordinlanguage(s, L))
        Learned.F += estate;

      for(auto a: L.alph) {
        lsetof<word> q2;
        for(auto e: E) If(wordinlanguage(concat(s,a,e), L))
          q2 += e;
        Learned.delta += mytransition(estate, a, elof(q2));
        }
      }
    }
  
  std::cout << "Guessed DFA:" << std::endl;
  printAutomaton(Learned);
  std::cout << std::endl;

  // (q1,q2) \in compare iff, after reading some word, the
  // DFA L is in state q1 and the DFA Learned
  // is in state q2
  lsetof<lpair<State, State>> compare;
  
  // for each (q1,q2) in compare, witnesses contains ((q1,q2), w),
  // where w is the witness word
  lsetof<lpair<lpair<State,State>, word> > witnesses;

  for(elem e: L.I) for(elem e2: Learned.I) {
    auto initpair = make_lpair(e, e2);
    compare += initpair;
    witnesses += make_lpair(initpair, word());
    }
  
  for(auto witness: witnesses) {

    elem q1 = (witness.first).first;
    elem q2 = (witness.first).second;
    
    word w = (witness.second);
    
    lbool m1 = memberof(q1, L.F);
    lbool m2 = memberof(q2, Learned.F);

    If(m1 ^ m2) {
      std::cout << "Counterexample: " << w << " where " << emptycontext << std::endl << std::endl;
      while(true) {

#ifdef ADD_COLUMNS
        If(!memberof(w, E)) E += w;
#else
        If(!memberof(w, S)) S += w;
#endif
        if(w.size() == 0) goto again;
      
#ifdef ADD_COLUMNS
        for(int i=0; i<int(w.size()-1); i++) w[i] = w[i+1];
#endif
        w.inner.resize(w.size() - 1);
        }
      }
    
    for(auto a: L.alph)
      for(mytransition& t1: L.delta) 
        If(t1.src == q1 && t1.symbol == a) 
      for(mytransition& t2: Learned.delta) 
        If(t2.src == q2 && t2.symbol == a) {
          auto p2 = make_lpair(t1.tgt, t2.tgt);
          If(!memberof(p2, compare)) {
            compare += p2;
            witnesses += make_lpair(p2, concat(w, a));
            }
          }
    }

  std::cout << "Learning successful!" << std::endl;
  }

#else
void learning(const myautomaton& L) {
  std::cout << "NFA to learn:" << std::endl;
  printAutomaton(L);
  std::cout << std::endl;

  lsetof<word> S, E;
  S += word();
  E += word();
  
  while(true) {
    again:
    std::cout << "S = " << S << std::endl;
    std::cout << "E = " << E << std::endl;
  
    lset uprows;
    
    for(word s: S) {
      lset ourrow;
      for(word e: E) If(wordinlanguage(concat(s,e), L))
        ourrow += e;
      
      If(!memberof(ourrow, uprows))
        uprows += ourrow;
      }
    
    std::cout << "uprows = " << uprows << std::endl;

    lset allrows = uprows;

    for(word s: S) for(elem a: sigma) {
      lset ourrow;
      for(word e: E) If(wordinlanguage(concat(s,a,e), L))
        ourrow += e;
      
      If(!memberof(ourrow, allrows))
        allrows += ourrow;
      }

    std::cout << "uprows = " << uprows << std::endl;
    std::cout << "allrows = " << allrows << std::endl;
    
    lset primerows;
    for(auto row: allrows) {
      lset aggregator;
      for(auto row2: allrows) 
        If(row2 != row && subseteq(asSet(row2), asSet(row)))
          aggregator |= asSet (row2);
      If(aggregator != row)
        primerows += row;
      }
    
    std::cout << "primerows = " << primerows << std::endl;

    for(word s: S) for(elem a: sigma) {
      lset ourrow;
      for(word e: E) If(wordinlanguage(concat(s,a,e), L))
        ourrow += e;
      
      If(memberof(ourrow, primerows) && !memberof(ourrow, uprows)) {
        std::cout << "Adding row: " << concat(s,a) << std::endl;
        S += concat(s,a);
        goto again;
        }
      }
    
    for(word s1: S) for(word s2: S) {
      lbool contained = true;
      for(word e: E)
        If(wordinlanguage(concat(s1,e),L) && !wordinlanguage(concat(s2,e),L))
          contained = false;
        
      If(contained) for(elem a: sigma) {
        for(word e: E)
          If(wordinlanguage(concat(s1,a,e),L) && !wordinlanguage(concat(s2,a,e),L)) {
          std::cout << "Adding column: " << concat(a,e) << std::endl;
          E += concat(a, e);
          goto again;
          }
        }
      }
        
    myautomaton learned;
    learned.Q = primerows & uprows;

    lset initrow;
    for(word e: E) If(wordinlanguage(e, L)) initrow += e;
        
    for(auto row: learned.Q) {
      If(subseteq(asSet (row), initrow)) learned.I += row;
      If(memberof(word(), asSet (row))) learned.F += row;
      }

    for(auto s: S) {
      lset rowS;
      for(word e: E) If(wordinlanguage(concat(s, e), L)) rowS += e;

      for(auto a: sigma) {
        lset rowsa;
        for(word e: E) If(wordinlanguage(concat(s, a, e), L)) rowsa += e;

        for(auto r: learned.Q) If(subseteq(asSet (r), rowsa))
          If(!memberof(transition(rowS, a, r), learned.delta))
            learned.delta += transition(rowS, a, r);
        }
      }
    
    std::cout << std::endl << "Guessed DFA:" << std::endl;
    printAutomaton(learned);
    std::cout << std::endl;

    lsetof<elpair> compare;
    lsetof<elpair> witnesses;
    
    compare += elpair(L.I, learned.I);
    witnesses += elpair(elpair(L.I, learned.I), word());
    
    for(auto wit: witnesses) {
      lset states1 = asSet(as<elpair> (wit.first).first);
      lset states2 = asSet(as<elpair> (wit.first).second);
      auto w = as<word> (wit.second);
      
      If(((states1 & L.F) != emptyset) ^ ((states2 & learned.F) != emptyset)) {
        std::cout << "Counterexample: " << w << " where " << emptycontext << std::endl << std::endl;
        while(true) {
  
  #ifdef ADD_COLUMNS
          If(!memberof(w, E)) E += w;
  #else
  #warn ADD_COLUMNS recommended for learnNFA
          If(!memberof(w, S)) S += w;
  #endif
          if(w.size() == 0) goto again;
        
  #ifdef ADD_COLUMNS
          for(int i=0; i<int(w.size()-1); i++) w[i] = w[i+1];
  #endif
          w.resize(w.size() - 1);
          }
        }

#ifdef LENGTH_LIMIT
      if(w.size() >= LENGTH_LIMIT) {
        std::cout << "Ignoring suffixes of: " << w << " where " << emptycontext << std::endl << std::endl;
        continue;
        }
#endif

      for(elem a: sigma) {
        lset nextstate1;
        lset nextstate2;
        for(auto t: L.delta) 
          If(t.symbol == a && memberof(t.src, states1))
            If(!memberof(t.tgt, nextstate1))
            nextstate1 += t.tgt;

        for(auto t: learned.delta) 
          If(t.symbol == a && memberof(t.src, states2))
            If(!memberof(t.tgt, nextstate2))
              nextstate2 += t.tgt;
            
        auto p = elpair(nextstate1, nextstate2);
        
        If(!memberof(p, compare)) {
          compare += p;
          witnesses += elpair(p, concat(w, a));
          }
        }
      }

    std::cout << "Learned correctly!" << std::endl;
    break;
    }  
  
  }
#endif

template<class T, class U> elem elpair(T x, U y) { return elof(make_lpair(x,y)); }

void buildStackAutomaton(myautomaton& target, lsetof<term>& A, elem& etrash, elem& epush, elem& epop, word w, int more) {
 
  auto ew = elof(w);
  target.Q += elof(w);
  target.F += elof(w);
  
  if(more) for(term a: A) {
    word w2 = concat(w, elof(a));
    auto ew2 = elof(w2);
    buildStackAutomaton(target, A, etrash, epush, epop, w2, more-1);
    target.delta += mytransition(ew, elof(make_lpair(epush, a)), ew2);
    target.delta += mytransition(ew2, elof(make_lpair(epop, a)), ew);
    for(term b: A) If(a != b)
      target.delta += mytransition(elof(w2), (elpair(epop, b)), etrash);
    }
  else for(term a: A) target.delta += mytransition(elof(w), (elpair(epush, a)), etrash);
  }

void buildQueueAutomaton(myautomaton& target, lsetof<term>& A, elem& etrash, elem& epush, elem& epop, word w, int more) {
 
  auto ew = elof(w);
  target.Q += ew;
  target.F += ew;
  
  if(w.size() == 0) {
    for(term a: A) target.delta += mytransition(ew, (elpair(epop, a)), etrash);
    }
  else {
    word w2;
    for(int i=1; i<(int) w.size(); i++) w2.push_back(w[i]);
    auto ew2 = elof(w2);
    target.delta += mytransition(ew, (elpair(epop, w[0])), ew2);
    for(term a: A) If(elof(a) != w[0]) target.delta += mytransition(ew, elof(elpair(epop, a)), etrash);
    }
  
  if(more) for(term a: A) {
    word w2 = concat(w, elof(a));
    auto ew2 = elof(w2);
    target.delta += mytransition(ew, elpair(epush, a), ew2);
    buildQueueAutomaton(target, A, etrash, epush, epop, w2, more-1);
    }
  else for(term a: A) target.delta += mytransition(ew, elpair(epush, a), etrash);
  }

void buildNthLast(myautomaton& target, lset& A, word w, int nth) {
  
  if(w.size() == nth) {    
    If(w[0] == w[1]) 
      target.F += elof(w);
    for(auto a: A) {
      word w2 = w;
      for(int i=1; i<nth-1; i++) w2[i] = w2[i+1];
      w2[nth-1] = a;
      target.delta += mytransition(elof(w), a, elof(w2));
      }
    }
  else {
    for(auto a: A) {
      word w2 = concat(w, a);
      target.delta += mytransition(elof(w), a, elof(w2));
      buildNthLast(target, A, w2, nth);
      }
    }    
  }

template<class T> lset be_lset(const lsetof<T>& s) {
  lset val;
  for(auto el: s) val += elof(el);
  return val;
  }
  
int main() {
  lasttime = getVa();
  initLois();
  // pushSolverDiagnostic("checking: ");

//  solver = solverBasic() || solverIncremental("cvc4 --lang smt --incremental");
//  solver = solverBasic() || solverIncremental("./z3 -smt2 -in");
//  solver = solverBasic() || solverIncremental("./cvc4-new --lang smt --incremental");

  Domain dA("Atoms");
  auto A = dA.getSet();

  lset sigma;
  for(auto a: A) sigma += elof(a);

  sym.neq = "≠";
  
  myautomaton target;
  target.alph = sigma;
  
  // language L1 from the paper (repeated letter)

#if DOUBLEWORD == 1
#define WHICH doubleword
  if(true) {
    elem e0 = elof(0);
    elem e1 = elof(1);
    elem e2 = elof(2);
    target.Q += e0;
    target.Q += e1;
    target.Q += e2;
    for(auto a:A) target.Q += a;
    target.F += e1;
    target.I += e0;
    for(auto a:A) target.delta += mytransition(e0, a, a);
    for(auto a:A) target.delta += mytransition(a, a, e1);
    for(auto a:A) for(auto b: A) If(a != b) target.delta += transition(a, b, e2);
    for(auto a:A) target.delta += mytransition(e2, a, e2);
    }                             
#endif

  // language L2 from the paper (repeated two letters: 'baba')

#if DOUBLEWORD == 2
#define WHICH doubleword
  if(true) {
    elem eini = elof(0);   // initial
    elem etrash = elof(1); // trash
    elem eaccept = elof(2); // accept
    target.Q += eini;
    target.Q += etrash;
    target.Q += eaccept;
    for(auto a: A) target.Q += elof(a); // read one letter
    for(auto a: A) for(auto b: A) target.Q += elpair(a,b); // read two letters
    for(auto a: A) target.Q += elpair(a,eini); // read three letters, wait for 'a'
    target.F += eaccept;
    target.I += eini;
    
    for(auto a: A) target.delta += mytransition(eini, elof(a), elof(a));
    for(auto a: A) for(auto b: A)
      target.delta += mytransition(elof(a), elof(b), elpair(a,b));
    for(auto a: A) for(auto b: A)
      target.delta += mytransition(elpair(a,b), elof(a), elpair(b,eini));
    for(auto a: A) for(auto b: A) for(auto c: A) If(a != c)
      target.delta += mytransition(elpair(a,b), elof(c), etrash);
    for(auto a: A) 
      target.delta += mytransition(elpair(a,eini), elof(a), eaccept);
    for(auto a: A) for(auto b: A) If(a != b)
      target.delta += mytransition(elpair(a,eini), elof(b), etrash);
    for(auto a: A) target.delta += mytransition(etrash, elof(a), etrash);
    for(auto a: A) target.delta += mytransition(eaccept, elof(a), etrash);
    }
#endif
  
#ifdef STACKSIZE
#define WHICH stacksize
  if(true) {
    int reps = 0;
    elem etrash = elof(0);
    elem epush = elof(1);
    elem epop = elof(2);
    target.alph = be_lset(newSet(epush, epop) * A);
    target.Q += etrash;
    target.I += elof(word());
    for(auto a: A) target.delta += mytransition(elof(word()), elpair(epop, elof(a)), etrash);
    
    buildStackAutomaton(target, A, etrash, epush, epop, word(), STACKSIZE);
    }
#endif
  
#ifdef QUEUESIZE
#define WHICH queuesize
  if(true) {
    int reps = 0;
    elem etrash = elof(0);
    elem epush = elof(1);
    elem epop = elof(2);
    target.alph = be_lset(newSet(epush, epop) * A);
    target.Q += etrash;
    target.I += elof(word());
    
    buildQueueAutomaton(target, A, etrash, epush, epop, word(), QUEUESIZE);
    }
#endif

#ifdef NTHLAST
#define WHICH nthlast
  if(true) {
    target.I += elof(word());
    buildNthLast(target, sigma, word(), NTHLAST);
    }
#endif

#ifdef EQLANG
#define WHICH eqlang

#ifndef LEARN_NFA
#error Learning NFA as DFA!
#endif

#ifndef LENGTH_LIMIT
#warn Will not work without LENGTH_LIMIT
#endif
  if(true) {
    printf("eq lang\n");
    elem einit = 0;
    elem eaccept = 1;
    target.Q += einit;
    target.Q += eaccept;
    target.I += einit;
    target.F += eaccept;
    for(elem a: A) {
      target.delta += transition(einit, a, a);
      target.delta += transition(einit, a, einit);
      target.delta += transition(a, a, eaccept);
      target.delta += transition(eaccept, a, eaccept);
      for(elem b: A) target.delta += transition(a, b, a);
      }
    }
#endif

  lsetof<word> allwords;

  for(auto a: target.alph) for(auto b: target.alph) {
    word w;
    w.push_back(a);
    w.push_back(b);
    If(wordinlanguage(w, target)) allwords += w;
    }

  std::cout << "All words of length 2: " << allwords << std::endl << std::endl;
  
  allwords = newSet<word>({});
  
  for(auto a: target.alph) for(auto b: target.alph) for(auto c: target.alph) for(auto d: target.alph) {
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
