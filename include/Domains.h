#ifndef _DOMAINS_H_
#define _DOMAINS_H_

#include "AMReX_Vector.H"
#include "Domain.h"

class Domains {
private:
  amrex::Vector<Domain*> domainVector;
  int iSelected;

public:
  Domains() = default;
  ~Domains() { clear(); }

  void add_new_domain() { domainVector.push_back(new Domain()); }

  void clear() {
    for (auto& dom : domainVector) {
      delete dom;
    }
    domainVector.clear();
  }

  Domain& operator()(int i) {
    if (i >= domainVector.size())
      amrex::Abort("Error: index out of range for class Domains!");

    return *(domainVector[i]);
  };

  int size(){
      return domainVector.size();
  }

  void select(int i){iSelected = i; }

  int selected()const{return iSelected;}
};

#endif
