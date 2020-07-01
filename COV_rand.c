#include <stdlib.h>
#include <stdio.h>

#if defined(USE_GETRANDOM)
# include <linux/random.h>
# include <errno.h>
#else /* if defined(USE_GETRANDOM) */
# include <time.h>
#endif /* if defined(USE_GETRANDOM) */

#include "common.h"

long seed = 555L;

int use_fixed_seed = 0;

/* wrap the random function so we can easily change the selected random
  generator */

void COV_init_rand() {
#if defined(USE_GETRANDOM)
  unsigned short sv[3];
#else /* if defined(USE_GETRANDOM) */
  struct timespec t;
#endif /* if defined(USE_GETRANDOM) */
  int ret;

  if (use_fixed_seed) {
    printf("Using fixed seed %lu\n", seed);
    srand48(seed);
  } else {
#if defined(USE_GETRANDOM)
    ret = getrandom((void *)sv, sizeof(sv), 0);

    while (ret != sizeof(sv) && (errno == EAGAIN || errno == EINTR)) {
      ret = getrandom((void *)sv, sizeof(sv), 0);
    }
    seed48(sv);
#else /* if defined(USE_GETRANDOM) */
    ret  = clock_gettime(CLOCK_MONOTONIC, &t);
    seed = t.tv_sec * t.tv_nsec;
    srand48(seed);
    printf("Using seed %lu\n", seed);
#endif /* if defined(USE_GETRANDOM) */
  }
}

/* Must produce a values in the [0, 1) range */
double COV_rand() {
  double ret;

  ret =  drand48();
  return ret;
}
