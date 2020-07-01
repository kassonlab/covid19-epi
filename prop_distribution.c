#include <stdio.h>
#include <stdlib.h>

#include "common.h"

/* Return a random index value in the range [lb, ub] based on probablility
  distribution dist */
/* dist must consist of monotonically increasing probablility values */
long rand_distr_indx(double *dist, long lb, long ub)
{
  long   i, ret = -1;
  double r = COV_rand();

  for (i = ub; i >= lb; i--) {
    if (r > dist[i]) {
      ret = i + 1;
      break;
    }
  }

  if (ret < 0) {
    ret = lb;
  }
  return ret;
}

/* Generate a monotonically increasing vector of probablilities */
int generate_inc_distr_vec(double **dist_vec, double *prob, int n, char *name)
{
  int i, ret = 0;
  double eps = 1e-10;

  *dist_vec      = (double *)malloc(n * sizeof(double));
  (*dist_vec)[0] = prob[0];

  for (i = 1; i < n; i++) {
    (*dist_vec)[i] = (*dist_vec)[i - 1] + prob[i];
  }

  /* Check that the sum is close enough to 1, we're dealing with doubles... */
  if (((*dist_vec)[n - 1] > 1.0 + eps) || ((*dist_vec)[n - 1] < 1.0 - eps)) {
    fprintf(stderr,
            "propbablilities for %s doesn't sum to 1 within %f\n",
            name,
            eps);
    ret = 1;
  }

  /* Force upper value to 1, needed to make sure no random value ends up outside
    of the distribution */
  (*dist_vec)[n - 1] = 1.0;

  return ret;
}
