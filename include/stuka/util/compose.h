//
// Created by Avinash Madavan on 2019-07-19.
//

#ifndef STUKA_UTIL_COMPOSE_H_
#define STUKA_UTIL_COMPOSE_H_

#include <utility>
#include <memory>

namespace stuka {

// Hack courtesy of https://stackoverflow.com/questions/43205257/iterate-over-c-variadic-template
template<class F>
struct compose_t {
  F f;

  template<class Lhs, class Rhs>
  inline friend auto operator+(compose_t<Lhs> lhs, compose_t<Rhs> rhs);

  template<class...Args>
  decltype(auto) operator()(Args &&...args) {
    return f(std::forward<Args>(args)...);
  }
};

template<class Lhs, class Rhs>
auto operator+(compose_t<Lhs> lhs, compose_t<Rhs> rhs) {
  auto r = [lhs = std::move(lhs).f, rhs = std::move(rhs).f](auto &&...args)
      -> decltype(auto) { return lhs(rhs(decltype(args)(args)...)); };
  return compose_t<decltype(r)>{std::move(r)};
}

template<class F>
compose_t<F> compose(F f) { return {std::forward<F>(f)}; }

template<class T>
auto toShared() {
  return [](auto &&...args) {
    return std::make_shared<T>(decltype(args)(args)...);
  };
}

}

#endif //STUKA_UTIL_COMPOSE_H_
