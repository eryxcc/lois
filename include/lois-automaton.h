namespace lois {

  template<class State, class Symbol> struct transition {
    State src;
    Symbol symbol;
    State tgt;
    transition(const State& _src, const Symbol& _symbol, const State& _tgt) { src=_src; symbol=_symbol; tgt=_tgt; }
    transition() {}
    };
  
  template<class State, class Symbol> transition<State,Symbol> make_transition(const State &q1, const Symbol& x, const State &q2) {
    return transition<State,Symbol> (q1,x,q2);
    }

  template<class State, class Symbol> std::ostream& operator << (std::ostream& os, transition<State,Symbol> t) { 
    return os << "(" << t.src << "," << t.symbol << "," << t.tgt << ")"; }

  template<class State, class Symbol> bool isused(vptr v, transition<State,Symbol> p) { 
    return isused(v, p.src) || isused(v, p.symbol) || isused(v, p.tgt); 
    }

  template<class State, class Symbol> transition<State,Symbol> substitute(const transition<State,Symbol>& p, const varsubstlist& l) { 
    return make_transition(substitute(p.src, l), substitute(p.symbol, l), substitute(p.tgt, l));
    }

  template<class State, class Symbol> rbool operator == (const transition<State,Symbol>& a, const transition<State,Symbol>& b) { 
    return a.src == b.src && a.symbol == b.symbol && a.tgt == b.tgt; }
  template<class State, class Symbol> rbool operator != (const transition<State,Symbol>& a, const transition<State,Symbol>& b) { 
    return a.src != b.src || a.symbol != b.symbol && a.tgt == b.tgt; }
  template<class State, class Symbol> rbool operator <  (const transition<State,Symbol>& a, const transition<State,Symbol>& b) { 
    lbool ret;
    Ife(a.src != b.src)
      ret = a.src < b.src;
    else Ife(a.symbol != b.symbol)
      ret = a.symbol < b.symbol;
    else ret = a.tgt < b.tgt;
    }
  template<class State, class Symbol> rbool operator <= (const transition<State,Symbol>& a, const transition<State,Symbol>& b) { 
    lbool ret;
    Ife(a.src != b.src)
      ret = a.src < b.src;
    else Ife(a.symbol != b.symbol)
      ret = a.symbol < b.symbol;
    else ret = a.tgt <= b.tgt;
    }

  
  template<class State, class Symbol> struct automaton {
    lsetof<Symbol> alph;
    lsetof<State> Q;
    lsetof<transition<State,Symbol>> delta;
    lsetof<State> F;
    lsetof<State> I;
    };

  
  };
