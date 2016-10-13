#include "../include/loisinternal.h"

namespace lois {

int indent = 0;

indenter iindent {2, false}, ieol {0, true}, iunindent {-2, false};

std::ostream& operator << (std::ostream& os, indenter a) { 
  indent += a.i;
  if(indent >= 0 && a.doit) {
    os << std::endl; for(int i=0; i<indent; i++) os << " ";
    }
  return os;
  }

OutputLanguage outlan = L0;

symbols sym;

void symbols::useASCII() {

  exists.set("ex ");
  forall.set("all ");
  _and.set(" && ");
  _or.set(" || ");
  eq.set(" == ");
  neq.set(" != ");
  in.set(" in ");
  
  emptyset.set("{}");
  ssunion.set("; ");
  leftbrace.set("{");
  rightbrace.set("}");
  sfepipe.set("|");
  sfecomma.set(", "); 

  leq.set(" <= ");
  geq.set(" >= ");
  edge.set(" ~ ");
  noedge.set(" !~ ");

  greater.set(" > ");
  less.set(" < ");
  min.set("^"); max.set("v");
  plus.set("+"); times.set("*"); minus.set("-"); divide.set("/");

  arrow.set(" ->> "); // " \u2192 ");
  notarrow.set(" />> "); // " \u219B ");
  pseudo.set("P");
  }

void symbols::useUnicode() {

  exists.set("\u2203 ");
  forall.set("\u2200 ");
  _and.set(" \u2227 ");
  _or.set(" \u2228 ");
  eq.set("=");
  neq.set(" \u2260 ");
  in.set("\u2208");
  
  emptyset.set("\u2205");
  ssunion.set("; "); // "} \u222a {";
  leftbrace.set("{");
  rightbrace.set("}");
  sfepipe.set("|");
  sfecomma.set(", ");  

  leq.set(" \u2264 ");
  geq.set(" \u2265 ");
  edge.set(" \u223c ");
  noedge.set(" \u2241 ");

  greater.set(" > ");
  less.set(" < ");
  min.set("\u2227");
  max.set("\u2228");
  plus.set("+"); times.set("*"); minus.set("-"); divide.set("/");

  arrow.set(" ->> "); // " \u2192 ");
  notarrow.set(" />> "); // " \u219B ");
  pseudo.set("\u2647");
  }

void symbols::useLaTeX() {

  exists.set("{\\exists}");
  forall.set("{\\forall}");
  _and.set(" {\\wedge} ");
  _or.set(" {\\vee} ");
  eq.set(".set(");
  neq.set(" {\\neq} ");
  in.set("{\\in}");
  
  emptyset.set("{\\emptyset}");
  ssunion.set(";\\ ");
  leftbrace.set("\\left\\{");
  rightbrace.set("\\right\\}");
  sfepipe.set("|");
  sfecomma.set(",\\ ");  

  leq.set(" {\\leq} ");
  geq.set(" {\\geq} ");
  edge.set(" {\\sim} ");
  noedge.set(" {\\not\\sim} ");
  greater.set(" > ");
  less.set(" < ");
  max.set("{\\vee}");
  min.set("{\\wedge}");
  plus.set("+"); times.set("{\\cdot}"); minus.set("-"); divide.set("/");

  arrow.set(" {\\succurlyeq} ");
  notarrow.set(" {\\not\\succurlyeq} ");
  pseudo.set("{\\pseudo}");
  }

RelOrder *mainOrder;

void initLois() { 
  sym.useUnicode();
  mainOrder = new RelOrder(sym.greater, sym.leq, sym.max, sym.min);
  useDefaultSolver(50000);
  emptyset = newSet();
  smtLogic = "LRA";
  }

}
