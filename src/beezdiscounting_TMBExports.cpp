// beezdiscounting TMB model dispatcher
// Registers the MixedDiscounting TMB model for the beezdiscounting package.

#define TMB_LIB_INIT R_init_beezdiscounting
#include <TMB.hpp>
#include "MixedDiscounting.h"
#include "ChoiceDiscounting.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model);
  if (model == "MixedDiscounting") {
    return MixedDiscounting(this);
  } else if (model == "ChoiceDiscounting") {
    return ChoiceDiscounting(this);
  } else {
    error("Unknown model");
  }
  return Type(0);
}
