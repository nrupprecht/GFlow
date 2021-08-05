#include <base/modifier.hpp>

using namespace GFlowSimulation;

Modifier::Modifier(GFlow *gflow)
    : Base(gflow), remove(false) {};

bool Modifier::getRemove() {
  return remove;
}

void Modifier::setRemove(bool r) {
  remove = r;
}
