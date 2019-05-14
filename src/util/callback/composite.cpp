//
// Created by Avinash Madavan on 2019-05-13.
//

#include <stuka/util/callback/composite.h>

stuka::util::callback::Composite::Composite(
    const std::vector<const std::shared_ptr<stuka::util::callback::BaseCallback>> &cbs_) : cbs_(cbs_) {}

void stuka::util::callback::Composite::initialize(const stuka::OptimizeState state) {
  for (std::vector<const std::shared_ptr<stuka::util::callback::BaseCallback>>::const_iterator it = cbs_.begin();
       it != cbs_.end(); ++it)
    (*it)->initialize(state);
}

void stuka::util::callback::Composite::callback(const stuka::OptimizeState state) {
  for (std::vector<const std::shared_ptr<stuka::util::callback::BaseCallback>>::const_iterator it = cbs_.begin();
       it != cbs_.end(); ++it)
    (*it)->callback(state);
}
