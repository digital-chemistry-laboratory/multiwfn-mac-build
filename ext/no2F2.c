#include <stdio.h>
#include <stdlib.h>

#include "2F2.h"

double hyp2F2(double a1, double a2, double b1, double b2, double z) {
  fprintf(stderr, "Fractional derivatives are not supported!");
  fprintf(stderr, "Please, recompile Multiwfn with its support!");
  exit(1);
  return 0.0;
}
