#include <stdlib.h>

#include "common.h"

/* wrap the random function so we can easily change the selected random generator */

void COV_init_rand() {
    srand48(555L);
}

/* Must produce a values in the [0, 1) range */
double COV_rand() {
    double ret;
    ret =  drand48();
    return ret;
}
