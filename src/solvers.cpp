#include "../include/loisinternal.h"

#include <fstream>
#include <string.h>

namespace lois {

// 0 = not satisfiable
// 1 = unknown
// 2 = satisfiable

// just crash when trying to solve something
struct SolverCrash : Solver {
  virtual int solve(rbool phi) { throw unsolvable_exception(); }
  };

// just assume that everything is satisfiable (2), unknown (1), or unsat (0)
struct SolverDummy : Solver {
  int answer;
  virtual int solve(rbool phi) {
    return answer;
    }
  SolverDummy(int i) : answer(i) {}
  };

solveptr solver = std::make_shared<SolverCrash> ();

// BASIC SOLVER: can only solve True and False

struct SolverBasic : Solver {

  virtual int solve(rbool phi) {
    if(phi.isTrue()) return 2;
    if(phi.isFalse()) return 0;
    return 1;
    }
  
  virtual int solveEnv() {
    rbool phi = outenv();
    if(phi.isTrue()) return 2;
    if(phi.isFalse()) return 0;
    return 1;
    }
  
  SolverBasic() {}
  };

// EXHAUSTIVE SOLVER: exhaustive search

bool notForInternalSolver(rbool phi) {
  featurelist fl;
  phi->listFeatures(fl);
  for(auto r: fl.relations) if(!r->worksWith(L0)) {
    std::cout << "The internal solver cannot use " << r->getName() << std::endl;
    return 1;
    }

  // only one tree allowed
  RelTree *tree = NULL;
  for(auto r: fl.relations) {
    RelTree *o = dynamic_cast<RelTree*> (r);
    if(tree && o && o != tree) {
      std::cout << "Too many trees" << std::endl;
      return true;
      }
    if(o) tree = o;
    }
  
  return false;
  }

struct SolverExhaustive : Solver {
  int tries;
  bool verbose;

  virtual int solve(rbool phi) {
    if(notForInternalSolver(phi)) return 1;
    checkid++;
    phi->initDomains();
    triesleft = tries;
    // std::cout << "{exhaustive solving " << phi << "} ";
    try {
      bool b = phi->verify();
      if(triesleft > 0) {
        if(verbose) printf("answer: %d (%d values tried)\n", b, tries-triesleft);
        return b ? 2 : 0;
        }
      }
    catch(unsolvable_exception) {
      if(verbose) printf("unable to verify\n");
      return 1;
      }
    return 1;
    }

  SolverExhaustive(int t, bool v = false) : tries(t), verbose(v) {}
  };

// SMT SOLVER: call a SMT-lib solver for help

struct SolverSMT : Solver {
  const std::string cmd;

  virtual int solve(rbool phi) {

    featurelist fl;
    phi->listFeatures(fl);
    if(fl.domains.size() != 1) {
      std::cout << "SMT cannot use multiple domains yet" << std::endl;
      return 1;
      }
    for(auto r: fl.relations) if(!r->worksWith(SMT)) {
      std::cout << "SMT cannot use " << r->getName() << std::endl;
      return 1;
      }

    std::ofstream ofs;

    ofs.open("smtquery.smt", std::ofstream::out);
    
    // set 4 for nice indentation
    // indent = 4; 
    indent = -100;
    
    ofs << "(benchmark test\n  :extrasorts (A)\n  :logic " << smtLogic << "\n  :formula \n    ";
    
    OutputLanguage sav = outlan;
    outlan = SMT;
    ofs << phi << std::endl << "    )" << std::endl;
    outlan = sav;
    ofs.close();
    // system("cvc3 -lang smt < smtquery.smt > smtquery.out");  
    if(system(cmd.c_str()))
      printf("system error\n");
    FILE *f = fopen("smtquery.out", "rt");
    char buf[500];
    char *b = fgets(buf, 500, f);
    if(!b) { printf("ERROR\n"); exit(3); }
    fclose(f);
    std::string ans = buf;
    if(ans == "sat\n") return 2;
    if(ans == "unsat\n") return 0;
    return 1;
    }
  
  SolverSMT(std::string c) : cmd(c) { }
  };

// CVC SOLVER: call the CVC3 solver for help

struct SolverCVC : Solver {
  const std::string cmd;

  virtual int solve(rbool phi) {

    featurelist fl;
    phi->listFeatures(fl);
    if(fl.domains.size() != 1) {
      std::cout << "CVC cannot use multiple domains yet" << std::endl;
      return 1;
      }
    for(auto r: fl.relations) if(!r->worksWith(CVC3)) {
      std::cout << "CVC cannot use " << r->getName() << std::endl;
      return 1;
      }

    std::ofstream ofs;
    ofs.open("smtquery.cvc", std::ofstream::out);
    ofs << "A:TYPE;" << std::endl;
    OutputLanguage sav = outlan;
    outlan = CVC3;
    ofs << "CHECKSAT " << phi << ";" << std::endl;
    outlan = sav;
    ofs.close();
    if(system(cmd.c_str()))
      printf("system error\n");
    FILE *f = fopen("smtquery.out", "rt");
    char buf[500];
    char *b = fgets(buf, 500, f);
    if(!b) { printf("ERROR\n"); exit(3); }
    fclose(f);
    std::string ans = buf;
    printf("Answer: %s", buf);
    if(ans == "Satisfiable.\n") return 2;
    if(ans == "Unsatisfiable.\n") return 0;
    // for CVC4
    if(ans == "sat\n") return 2;
    if(ans == "unsat\n") return 0;
    return 1;
    }
  
  SolverCVC(std::string c = "cvc3 < smtquery.cvc > smtquery.out") : cmd(c) { }
  };

// SPASS SOLVER: call the SPASS solver for help

struct SolverSPASS : Solver {
  std::string cmd;

  virtual int solve(rbool phi) {
    std::ofstream ofs;

    featurelist fl;
    phi->listFeatures(fl);
    for(auto r: fl.relations) if(!r->worksWith(SPASS)) {
      ofs << "SPASS cannot use " << r->getName() << std::endl;
      return 1;
      }

    ofs.open("query.spass", std::ofstream::out);
    ofs << "begin_problem(LoisProblem)." << std::endl;
    ofs << std::endl;

    ofs << "list_of_descriptions." << std::endl;
    ofs << "name({*LOIS*})." << std::endl;
    ofs << "author({*LOIS*})." << std::endl;
    ofs << "status(unknown)." << std::endl;
    ofs << "description({*generated by LOIS*})." << std::endl;
    ofs << "end_of_list." << std::endl;
    ofs << std::endl;

    OutputLanguage sav = outlan;

    outlan = SPASS_SYMBOLS;
    ofs << "list_of_symbols." << std::endl;
    for(auto r: fl.relations) {
      RelOrder *o = dynamic_cast<RelOrder*> (r);
      if(o) ofs << "predicates[(" << spassescape(o->opgreater) << ", 2)]." << std::endl;
      }
    for(auto d: fl.domains) 
      ofs << "sorts[" << d->name << "]." << std::endl;
    ofs << "end_of_list." << std::endl;
    ofs << std::endl;

    outlan = SPASS_AXIOMS;
    ofs << "list_of_formulae(axioms)." << std::endl;
    for(auto r: fl.relations) {
      RelOrder *o = dynamic_cast<RelOrder*> (r);
      // to do: axiomatize min/max
      // to do: specify domains
      // (BUT: it does not work effectively with order even now)
#define G << spassescape(o->opgreater) <<
      ofs << "formula(forall([x,y,z],implies(and(" G "(x,y), " G "(y,z)), " G "(x,z))),1)." << std::endl;
      ofs << "formula(forall([x,y],not(and(" G "(x,y), equal(x,y)))),2)." << std::endl;
      ofs << "formula(forall([x,y],or(" G "(x,y), " G "(y,x), equal(x,y))),3)." << std::endl;
      ofs << "formula(forall([x],exists([y], " G "(y,x))),4)." << std::endl;
      ofs << "formula(forall([x],exists([y], " G "(x,y))),5)." << std::endl;
      ofs << "formula(forall([x,z], implies(" G "(x,z), exists([y], and(" G "(x,y), " G "(y,z))))), 6)." << std::endl;
#undef G
      }
    ofs << "end_of_list." << std::endl;
    ofs << std::endl;

    outlan = SPASS;
    ofs << "list_of_formulae(conjectures)." << std::endl;
    ofs << "formula(not(" << phi << "))." << std::endl;
    ofs << "end_of_list." << std::endl;
    ofs << std::endl;

    ofs << "end_problem." << std::endl;

    outlan = sav;
    ofs.close();
    if(system(cmd.c_str()))
      printf("system error\n");
    FILE *f = fopen("query.out", "rt");
    char buf[500];
    int ans = 1;
    while(fgets(buf, 500, f)) {
      if(strstr(buf, "Proof found")) ans=0;
      if(strstr(buf, "Completion found")) ans=2;
      }
    fclose(f);
    printf("SPASS answer: %d\n", ans);
    return ans;
    }
  
  SolverSPASS(std::string c = "SPASS query.spass > query.out") : cmd(c) { }
  };

// diagnostic

#include <sys/time.h>
long long getVa() {
  struct timeval tval;
  gettimeofday(&tval, NULL);
  return tval.tv_sec * 1000000 + tval.tv_usec;
  }

struct SolverVerbose : Solver {
  std::string msg;
  solveptr s;
  
  virtual int solve(rbool phi) {
    std::cout << ido << msg << phi << std::endl;
    return s->solve(phi);
    }
  
  virtual int solveEnv() {
    std::cout << ido << msg << emptycontext << std::endl;
    std::cout << ido << msg << "OUTENV=" << outenv() << std::endl;
    return s->solveEnv();
    }
  
  SolverVerbose(std::string m, solveptr _s) : msg(m), s(_s) { }
  };

struct SolverNamed : Solver {
  std::string name;
  solveptr s;
  
  virtual int solve(rbool phi) {
//  std::cout << "trying " << name << ": " << phi << std::endl; fflush(stdout);
    long long t = getVa();
    int r = s->solve(phi);
    std::cout << ido;
    printf("result of %s: %d, %10lld \u03bcs\n", name.c_str(), r, getVa() - t);
    return r;
    }
  
  virtual int solveEnv() {
//  std::cout << "trying " << name << ": " << emptycontext << std::endl; fflush(stdout);
    long long t = getVa();
    int r = s->solveEnv();
    std::cout << ido;
    printf("result of %s: %d, %10lld \u03bcs\n", name.c_str(), r, getVa() - t);
    return r;
    }
  
  SolverNamed(std::string _n, solveptr _s) : name(_n), s(_s) { }
  };

struct SolverCompare : Solver {
  std::vector<solveptr> solvers;

  virtual int solve(rbool phi) {
    int endr = 1;
    for(solveptr s: solvers) {
      int r = s->solve(phi);
      if(endr == 3) continue;
      else if(endr == 1) endr = r;
      else if(r != 1 && r != endr) endr = 3;
      }
    if(endr == 3) throw inconsistency_exception();
    return endr;
    }
  
  virtual int solveEnv() {
    int endr = 1;
    for(solveptr s: solvers) {
      int r = s->solveEnv();
      if(endr == 3) continue;
      else if(endr == 1) endr = r;
      else if(r != 1 && r != endr) endr = 3;
      }
    if(endr == 3) throw inconsistency_exception();
    return endr;
    }
  
  void addSolver(solveptr s) { solvers.push_back(s); }

  SolverCompare() { }
  };

// stack of solvers

struct SolverStack : Solver {
  solveptr s, fallback;

  virtual int solve(rbool phi) {
    int a = s->solve(phi);
    if(a != 1) return a;
    return fallback->solve(phi);
    }

  virtual int solveEnv() {
    int a = s->solveEnv();
    if(a != 1) return a;
    return fallback->solveEnv();
    }

  SolverStack(solveptr _s, solveptr _fb) : s(_s), fallback(_fb) { }
  };

// INCREMENTAL SMT SOLVER: call a SMT-lib solver for help,
// in interactive incremental mode

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>

#define BUFSIZE 20480

// return the ancestor of 'of' which is a child of 'a', or emptycontext otherwise
// useful to traverse the tree of contexts
contextptr findchildcontext(contextptr a, contextptr of) {
  while(of != emptycontext && of->parent != a) of = of->parent;
  return of;
  }

static int inccounter = 0;

bool wearechild = false;

std::string smtsort() {
  if(smtLogic == "LRA" || smtLogic == "NRA")
    return "Real";
  if(smtLogic == "LIA" || smtLogic == "NIA")
    return "Int";
  throw other_exception("logic unknown to LOIS");
  }

std::string smtLogic = "LRA";

std::string lastsolver;

struct SolverIncremental : Solver {

  std::ofstream ofs;

  int solver_input[2], solver_output[2], child, filepipe;
  
  char buf[60];

  bool crashed;
  
  contextptr last;
  
  int ourindent;
  
  bool isNotInSMTLogic(rbool phi) {
    featurelist fl;
    phi->listFeatures(fl);

//  for(auto d: fl.domains) std::cout << d->name << " ";
//  for(auto r: fl.relations) std::cout << r->getName() << " "; std::cout << std::endl;
    
    if(fl.domains.size() != 1) {
      std::cout << "SMT cannot use multiple domains yet" << std::endl;
      return true;
      }
    for(auto r: fl.relations) if(!r->worksWith(SMT_INC)) {
      std::cout << "SMT cannot use " << r->getName() << std::endl;
      return true;
      }
    
    // only one order allowed
    RelOrder *order = NULL;
    for(auto r: fl.relations) {
      RelOrder *o = dynamic_cast<RelOrder*> (r);
      if(o) {
        if(order && o != order) {
          std::cout << "Only one order relation allowed for " << smtLogic << std::endl;
          return true;
          }
        if(o) order = o;
        }
      }

    if(smtsort() == "Int") for(auto r: fl.relations) {
      RelReal *o = dynamic_cast<RelReal*> (r);
      if(o) {
        std::cout << smtLogic << " does not support RelReal" << std::endl;
        return true;
        }
      }
    
    if(smtsort() == "Real") for(auto r: fl.relations) {
      RelInt *o = dynamic_cast<RelInt*> (r);
      if(o) {
        std::cout << smtLogic << " does not support RelInt" << std::endl;
        return true;
        }
      }
    
    return false;
    }

  virtual int solve(rbool phi) { 
  
    if(isNotInSMTLogic(phi)) return 1;

    OutputLanguage sav = outlan;

//    std::string s = readline();
//    if(s != "") printf("{%s}\n", s.c_str());

    outlan = SMT_INC;

    ofs << "(push 1)" << std::endl;
    ofs << "(assert\n    " << phi << std::endl << "    )" << std::endl;
    ofs << "(check-sat)" << std::endl;
    ofs << "(pop 1)" << std::endl;

    outlan = sav;
    if(crashed) return 1;
    pump();
    while(!crashed) {
      std::string s = readline();
      if(s == "sat") return 2;
      if(s == "unsat") return 0;
      if(s == "unknown") return 1;
      }
    return 1;
    }
  
  virtual int solveEnv() {
#ifdef DO_NOT_USE_STACK
    return solve(outenv());
#endif

    if(isNotInSMTLogic(outenv())) return 1;

    OutputLanguage sav = outlan;
    outlan = SMT_INC;
    
    std::swap(indent, ourindent);
    
    while(last != currentcontext) {
      contextptr p = findchildcontext(last, currentcontext);
      if(p == emptycontext) {
        last = last->parent;
        ofs << "(pop 1)" << iunindent << ieol;
        }
      else {
        last = p;
        ofs << "(push 1)" << iindent << ieol;
        for(auto v: p->var)
          ofs << "(declare-const " << v << " " << smtsort() << ")" << ieol;
        if(!p->phi.isTrue())
          ofs << "(assert" << iindent << ieol << p->phi << ieol << ")" << iunindent << ieol;
        }
      }

    ofs << "(check-sat)" << ieol;

    std::swap(indent, ourindent);
    
    outlan = sav;
    if(crashed) return 1;
    pump();
    bool haderror = false;
    while(!crashed) {
      std::string s = readline();
      if(s == "sat") return haderror ? 1 : 2;
      if(s == "unknown") return 1;
      if(s == "unsat") return haderror ? 1 : 0;
      if(s.find("error") != std::string::npos) {
        haderror = true;
        printf("error reported: %s\n", s.c_str());
        }
      }
    return 1;
    }
  
  SolverIncremental(std::string c) { 
    // fprintf(stderr, "[create incremental %d]\n", inccounter);
    sprintf(buf, "incremental%d.smt", inccounter++);
    lastsolver = buf;

    crashed = false; last = emptycontext;
    // ofs << "(declare-sort A 0)\n";

    ourindent = 0;
    
    if(pipe(solver_input))
      throw other_exception("pipe input");
    if(pipe(solver_output))
      throw other_exception("pipe output");
    
    child = fork();
    if(child < 0) throw other_exception("fork");
    else if(child == 0) {
      wearechild = true;
      close(solver_input[1]);
      close(solver_output[0]);
      dup2(solver_output[1], 1);
      dup2(solver_input[0], 0);
      int err = system(c.c_str());
//    fprintf(stderr, "[kill %s]\n", buf);
      exit(0);
      }
    
    close(solver_input[0]);
    close(solver_output[1]);

    ofs.open(buf, std::ofstream::out);
    ofs << "(set-logic " << smtLogic << ")\n";
    ofs << "(define-fun min ((a "<<smtsort()<<") (b "<<smtsort()<<")) "<<smtsort()<<" (ite (<= a b) a b))\n";
    ofs << "(define-fun max ((a "<<smtsort()<<") (b "<<smtsort()<<")) "<<smtsort()<<" (ite (>= a b) a b))\n";
//  fprintf(stderr, "[print logic %s]\n", buf);

    filepipe = open(buf, O_RDONLY);
    }

  void pump() {
    char buf[BUFSIZE];
    int numbytes = 0;
    while(true) {
      ofs.flush();
      int bytes = read(filepipe, buf, BUFSIZE);
      if(bytes <= 0) break;
      numbytes += bytes;
      if(write(solver_input[1], buf, bytes) != bytes)
        throw other_exception("writing to solver");
      }
    if(numbytes)
      printf("Sent %dB to the solver\n", numbytes);
    }

  std::string readline() {
    std::string ans;
    while(true) {
      int status;
      int wpr = waitpid(child, &status, WNOHANG);
      if(wpr == child) {
        crashed = true;
        return "";
        }

      char buf;
      int bytes = read(solver_output[0], &buf, 1);
      if(bytes<0 && errno == EAGAIN) {
        // printf("AGAIN\n");
        continue;
        }
      if(bytes > 0) {
        if(buf == 10) return ans;
        else ans += buf;
        }
      if(bytes < 0) { crashed = true; return ""; }
      }
    }
  
  ~SolverIncremental() {
    if(wearechild || crashed) {
      close(solver_input[1]);
      close(solver_output[0]);
      return;
      }
    ofs << "(exit)\n" << std::endl;
//  fprintf(stderr, "[exit %s]\n", buf);
    ofs.close();
//  fprintf(stderr, "[close %s]\n", buf);
    pump();
    int status;
    waitpid(child, &status, 0);
    close(filepipe);
    close(solver_input[1]);
    close(solver_output[0]);
//  kill(child, SIGTERM);
/*  while(true) {
      std::string s = readline();
      if(s == "") break;
      printf("bye %s\n", s.c_str());
      } */
    }
  };

// create solvers

solveptr solverDummy(int i) { return std::make_shared<SolverDummy> (i); }
solveptr solverCrash() { return std::make_shared<SolverCrash> (); }

solveptr solverBasic() { 
  return std::make_shared<SolverBasic> (); 
  }
 
solveptr solverExhaustive(int t, bool v) { 
  return std::make_shared<SolverExhaustive> (t, v); 
  }

solveptr solverVerbose(std::string s, solveptr sv) { 
  return std::make_shared<SolverVerbose> (s, sv); 
  }

solveptr solverIncremental(std::string t) { 
  return std::make_shared<SolverIncremental> (t);
  }

solveptr solverStack(solveptr s, solveptr fallback) { 
  return std::make_shared<SolverStack> (s, fallback);
  }

solveptr solverSMT(std::string t) {
  return std::make_shared<SolverSMT> (t); 
  }

solveptr solverCVC() { 
  return std::make_shared<SolverCVC> ();
  }

solveptr solverCVC(std::string t) { 
  return std::make_shared<SolverCVC> (t); 
  }

solveptr solverSPASS(std::string t) { return std::make_shared<SolverSPASS> (t); }
solveptr solverSPASS() { return std::make_shared<SolverSPASS> (); }

solveptr solverCompare(std::initializer_list<solveptr> p) {
  auto a = std::make_shared<SolverCompare> ();
  for(auto x: p) a->addSolver(x);
  return a;
  }
  
solveptr solverNamed(std::string n, solveptr s) {
  return std::make_shared<SolverNamed> (n,s);
  }
  
extern int aggressive_simplification_tries;

void useDefaultSolver(int i, int j) {
  aggressive_simplification_tries = j;
  solver = solverExhaustive(j, true) || solverCrash();
  if(i > 0) solver = solverExhaustive(i, false) || solver;
  solver = solverBasic() || solver;
  }

std::string spassescape(symbol tt) {
  std::string res;
  std::string t = tt.asString();
  for(int i=0; i<int(t.size()); i++)
    if(isalpha(t[i]))
      res += t[i];
    else {
      int v = t[i];
      res += "_";
      if(v<0) res += 'N', v = -v;
      while(v>0) { res += 'a' + (v%26); v /= 26; }
      res += "_";
      }
  return res;
  }

}
