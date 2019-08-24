//
// Created by Avinash Madavan on 2019-05-13.
//

#include <stuka/util/callback/composite.h>

stuka::util::callback::Composite::Composite(
    const std::vector<std::shared_ptr<stuka::util::callback::BaseCallback>> &cbs) : cbs_(cbs) {}

void stuka::util::callback::Composite::initialize(const stuka::OptimizeState state) {
  for (std::vector<std::shared_ptr<stuka::util::callback::BaseCallback>>::const_iterator it = cbs_.begin();
       it != cbs_.end(); ++it)
    (*it)->initialize(state);
}

void stuka::util::callback::Composite::callback(const stuka::OptimizeState state) {
  for (std::vector<std::shared_ptr<stuka::util::callback::BaseCallback>>::const_iterator it = cbs_.begin();
       it != cbs_.end(); ++it)
    (*it)->callback(state);
}

void stuka::util::callback::Composite::finish(const stuka::OptimizeState state) {
  for (std::vector<std::shared_ptr<stuka::util::callback::BaseCallback>>::const_iterator it = cbs_.begin();
       it != cbs_.end(); ++it)
    (*it)->finish(state);
}
