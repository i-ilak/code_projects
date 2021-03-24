#ifndef PENNA_VECTOR_H
#define PENNA_VECTOR_H

#include <vector>

template <class T>
class penna_vector : public std::vector<T> {

  using base_t = std::vector<T>;

public:

  using typename base_t::value_type;
  using typename base_t::iterator;
  using typename base_t::size_type;

  using base_t::begin;
  using base_t::end;
  using base_t::back;
  using base_t::pop_back;

  penna_vector(size_type s=0, value_type x=value_type()) : base_t(s, x) {}

  template <class IT>
  penna_vector(IT first, IT last) : base_t(first, last) {}

  template <class P>
  void remove_if(P pred) {
    iterator it = begin();
    while(it != end()) {
      if(pred(*it)) {
        *it = back(); // copy the last to this location
        pop_back(); // remove the last one
      }
      else
        ++it;
    }
  }

};

#endif
