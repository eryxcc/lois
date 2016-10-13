#define PRES
#include "../include/loisextra.h"
#include <string.h>

using namespace lois;
using namespace std;

struct stackviewer {
    lset *Z;
    string name;
    contextptr e;
    stackviewer *v;
};

stackviewer *sv = NULL;

namespace lois {
    extern int vids;
}

struct addtostackviewer {
    stackviewer t;
    addtostackviewer(lset& s, string n) {
        t.Z = &s;
        t.name = n;
        t.e = currentcontext;
        t.v = sv;
        sv = &t;
    }
    ~addtostackviewer() { sv = t.v; }
};

struct stackviewere {
    elem Z;
    string name;
    contextptr e;
    stackviewere *v;
};

stackviewere *sve = NULL;

struct addtostackviewere {
    stackviewere t;
    addtostackviewere(elem x, string n) {
        t.Z = x;
        t.name = n;
        t.e = currentcontext;
        t.v = sve;
        sve = &t;
    }
    ~addtostackviewere() { sve = t.v; }
};

int pageid = 0;

void viewEnv(string s) {
    // system("clear");
    // 30 lines empty!
    printf("\\begin{frame}\n");
    printf("\\frametitle{Simulation of the stack}\n");

    printf("\\parbox[b][13em][b]{15em}{$\\begin{array}[b]{|c|}\n");


    contextptr e = currentcontext;
    stackviewer *csv = sv;
    stackviewere *csve = sve;
    while(true) {
        if(csv != NULL && csv->e == e) {
            // int svid = vids; vids = 10;

            cout << "\\hline \\makebox[12.5em]{\\color{local}$" << csv->name << " = " << *(csv->Z) << "$} \\\\" <<endl;

            // cout << "    LOCAL      " << csv->name << " = " << *(csv->Z) << endl;
            // vids = svid;
            csv = csv->v;
        }
        else if(csve != NULL && csve->e == e) {
            cout << "\\hline \\makebox[12.5em]{\\color{local}$" << csve->name << " = " << (csve->Z) << "$} \\\\" <<endl;

            // cout << "    LOCAL      " << csve->name << " = " << csve->Z << endl;
            csve = csve->v;
        }
        else if(e == NULL)
            break;
        else {
          if(!e->phi.isTrue())
            cout << "\\hline \\makebox[12.5em]{\\color{constraint}$" << e->phi << "$} \\\\" <<endl;
//          cout << "    CONSTRAINT " << e->phi << endl;
          for(auto v: e->var) 
//          cout << "    VARIABLE   " << v << sym.in << "A" << endl;
            cout << "\\hline \\makebox[12.5em]{\\color{variable}$" << v << sym.in << "\\bbA" << "$} \\\\" <<endl;
          e = e->parent;
        }
    }
//  cout << "    STACK      ########" << endl << endl;
    char buf[120];
    sprintf(buf, "talk/page%04d.cpp", pageid);
    FILE *g = fopen(buf, "wt"); // or screen!
    FILE *f = fopen("tests/stackdemo.cpp", "rt");
    bool on = false;
    bool here = false;
    while(!feof(f)) {
        char buf[512];
        fgets(buf, 512, f);
        if(buf[0] == '/' && strstr(buf, "END")) on = false;
        if(strstr(buf, "viewEnv"))
            here = strstr(buf, s.c_str());
        else if(strstr(buf, "addtostack"))
            ;
        else if(on) {
            if(here) fprintf(g, "/**/");
            else fprintf(g, "    ");
            here = false;
            fprintf(g, "%s", buf);
        }
        if(buf[0] == '/' && strstr(buf, "START")) on = true;
    }
  fclose(f); fclose(g);

  printf(
    "\\hline\\end{array}$}\\parbox[b][13em][b]{18em}{\\lstinputlisting{page%04d.cpp}}\n", pageid);
  printf("\\vskip 2em\n{\\color{constraint}constraints}, {\\color{local}local variables}, {\\color{variable}pseudoparalell variables}\n");
  printf("\\end{frame}\n");

  // "\\hline\\end{tabular}&\\lstinputlisting{page%04d.cpp}\\end{tabular}\\end{frame}\n",
    
  pageid++;
  // getchar
}

int main() {
    initLois();
#ifdef PRES
  sym.useLaTeX();
#endif

//Domain dA("A");
  Domain dA("\\bbA");
  lset A = dA.getSet();
    addtostackviewer a0(A, "A");

// START
  viewEnv("1");
  lset C;
  addtostackviewer a1(C, "C");
  viewEnv("2");
  for(elem x: A) {
    addtostackviewere a3(x, "x");
    viewEnv("3");
    lset D;
    addtostackviewer a2(D, "D");
    viewEnv("4");
    for(elem y: A) {
      addtostackviewere a4(y, "y");
      viewEnv("5");
      If(y != x) {
        viewEnv("0a");
        D += newSet(y);
        viewEnv("6");
        }
      viewEnv("7");
      }
    viewEnv("0b");
    C += D;
    viewEnv("8");
    }
  viewEnv("9");
// END
  return 0;
  }
