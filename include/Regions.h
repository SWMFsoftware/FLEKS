#ifndef _REGIONS_H_
#define _REGIONS_H_

#include <sstream>

#include "Shape.h"

class Regions {
public:
  Regions() = default;

  Regions(amrex::Vector<std::unique_ptr<Shape> >& in, std::string list) {
    define(in, list);
  }

  void define(amrex::Vector<std::unique_ptr<Shape> >& in, std::string list) {
    for (auto& shape : in) {
      shapes.push_back(shape.get());
    }

    includeList.clear();
    excludeList.clear();

    std::stringstream ss(list);
    std::string token;
    while (ss >> token) {
      if (token.find("+") != std::string::npos) {
        includeList.push_back(token);
      } else if (token.find("-") != std::string::npos) {
        excludeList.push_back(token);
      }
    }
  }

  // Does the 'shape' in the 'list'?
  bool contains(const Shape* shape,
                const amrex::Vector<std::string>& list) const {
    bool doContain = false;

    for (auto& l : list) {
      if (l == shape->get_name()) {
        doContain = true;
        break;
      }
    }

    return doContain;
  }

  // Is the 'shape' in the include list?
  bool is_included(const Shape* shape) const {
    return contains(shape, includeList);
  }

  // Is the 'shape' in the exclude list?
  bool is_excluded(const Shape* shape) const {
    return contains(shape, excludeList);
  }

  // Is the point 'xyz' inside the include list but outside the exclude list?
  bool is_inside(amrex::Real* xyz) const {
    bool isIncluded = false;
    for (auto& shape : shapes) {
      if (is_included(shape))
        isIncluded = shape->is_inside(xyz);

      if (isIncluded)
        break;
    }

    bool isExcluded = false;
    for (auto& shape : shapes) {
      if (is_excluded(shape))
        isExcluded = shape->is_inside(xyz);

      if (isExcluded)
        break;
    }

    return isIncluded && !isExcluded;
  }

private:
  amrex::Vector<Shape*> shapes;
  amrex::Vector<std::string> includeList;
  amrex::Vector<std::string> excludeList;
};

#endif