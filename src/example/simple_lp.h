//
// Created by Avinash Madavan on 7/16/18.
//

#ifndef SWITCHBACK_EXAMPLE_SIMPLE_LP_H
#define SWITCHBACK_EXAMPLE_SIMPLE_LP_H

#include "example.h"

namespace stuka { namespace example {

    class SimpleLP : public LinearProgramExample {
    public:
        stuka::LP::LinearProgram gen() override;
        std::string name() override {return "Simple LP Example";}
    };

}}


#endif //SWITCHBACK_EXAMPLE_SIMPLE_LP_H
