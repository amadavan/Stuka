//
// Created by Avinash Madavan on 12/1/18.
//

#include <string>
#include <iostream>

#include <stuka/stuka.h>

#include "example/simple_lp.h"
#include "example/simple_qp.h"
#include "example/borrelli_729.h"
#include "example/borrelli_731.h"

std::string print_header(const std::string &header, const char sep = '=') {
  size_t n_chars = header.length();
  size_t l_char = 40 - (size_t) ceil(n_chars/2.) - 1;
  size_t r_char = 40 - (size_t) floor(n_chars/2.) - 1;

  std::string tmp;
  for (size_t i = 0; i < l_char; ++i) tmp += sep;
  tmp += ' ' + header + ' ';
  for (size_t i = 0; i < r_char; ++i) tmp += sep;
  return tmp;
}

int main() {

  std::vector<stuka::example::LinearProgramExample*> ex_lp = {new stuka::example::SimpleLP()};
  std::vector<stuka::example::QuadraticProgramExample*> ex_qp = {new stuka::example::SimpleQP()};
  std::vector<stuka::example::DecomposedLinearProgramExample*> ex_dlp = {new stuka::example::Borrelli729(),
                                                                         new stuka::example::Borrelli731()};

  stuka::OptimizeState res;

  for (stuka::example::LinearProgramExample *ex : ex_lp) {
    std::cout << print_header(ex->name()) << std::endl;
    std::cout << print_header("Default", '-') << std::endl;
    res = stuka::util::linprog(ex->gen());
    std::cout << res.x.transpose() << "\t\tRuntime: " << res.runtime << std::endl;
  }

  for (stuka::example::QuadraticProgramExample *ex : ex_qp) {
    std::cout << print_header(ex->name()) << std::endl;
    std::cout << print_header("Default", '-') << std::endl;
    res = stuka::util::quadprog(ex->gen());
    std::cout << res.x.transpose() << "\t\tRuntime: " << res.runtime << std::endl;
  }

  for (stuka::example::DecomposedLinearProgramExample *ex : ex_dlp) {
    std::cout << print_header(ex->name()) << std::endl;
    std::cout << print_header("Default LP", '-') << std::endl;
    res = stuka::util::linprog(ex->full());
    std::cout << res.x.transpose() << "\t\tRuntime: " << res.runtime << std::endl;

    std::cout << print_header("Critical Region Exploration", '-') << std::endl;
    stuka::Options opts;
    opts.dlp_solver = stuka::CRE;
    res = stuka::util::linprog(ex->gen(), opts);
    std::cout << res.x.transpose() << "\t\tRuntime: " << res.runtime << std::endl;

    std::cout << print_header("Bender's Decomposition", '-') << std::endl;
    res = stuka::util::linprog(ex->gen());
    std::cout << res.x.transpose() << "\t\tRuntime: " << res.runtime << std::endl;

  }
}