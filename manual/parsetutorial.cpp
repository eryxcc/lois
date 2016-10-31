// tutorial parser, by Eryk Kopczynski

// This program parses the source of the Tutorial (tests/tutorial.cpp)
// and its output (out/tutorial.txt) in order to create the snippet
// files (manual/snippets/...) that are later included by LaTeX.

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;

bool firstline = true, lastempty = true;

int main(int argc, char **argv) {

  ofstream outf;
  
  ifstream inf;
  inf.open(argv[1]);
  if(!inf.is_open()) exit(1);
  outf.open("manual/snippets/output-pre");

  while(!inf.eof()) {
    string l;
    getline(inf, l);
    if(l.substr(0, 8) == "SECTION ") {
      outf.close();
      outf.open("manual/snippets/output-" + l.substr(8));
      }
    
    else outf << l << endl;
    }
  
  inf.close();
  outf.close();

  inf.open(argv[2]);
  outf.open("manual/snippets/snippet-pre");

  while(!inf.eof()) {
    string l;
    getline(inf, l);
    if(l.substr(0, 9) == "SECTION(\"") {
      outf.close();
      l = l.substr(9);
      l = l.substr(0, l.find("\""));
      outf.open("manual/snippets/snippet-" + l);
      firstline = true; lastempty = false;
      }
    
    else if(l.substr(0, 10) == "FUNCTION(\"") {
      outf.close();
      l = l.substr(10);
      l = l.substr(0, l.find("\""));
      outf.open("manual/snippets/snippet-" + l);
      firstline = true;
      }
    
    else {
      for(size_t i=0; i<l.size(); i++) if(l[i] == '\n') l[i] = '*';
      if(firstline && l == "") lastempty = false;
      else if(l == "") lastempty = true;
      else { 
        if(lastempty) outf << endl, lastempty = false;
        outf << l << endl;
        }
      firstline = false;
      }
    }
  
  
  inf.close();
  outf.close();

  outf.open("manual/snippets/ok");
  outf << "ok" << endl;
  outf.close();

  }
