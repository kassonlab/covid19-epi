#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "prop_distribution.h"
#include "COV_rand.h"
#include "age_dist.h"

void age_dist (double * age, int population, FILE* stats, int * age_distrib) {

  int i; // Counter
  int ret, age_sz;

  int age_start[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
  double age_dist[] = {0.10721, 0.11553, 0.12439, 0.13475, 0.12578, 0.12704, 0.10782, 0.09952, 0.05796};
  double *age_dst;

  age_sz = (int)(sizeof(age_dist)/sizeof(age_dist[0]));
  ret = generate_inc_distr_vec(&age_dst, age_dist, age_sz, "age");
  if (ret) {
    fprintf(stderr, "Bailing out on age initialization\n");
    exit(1);
  }

  for (i=0; i<population; i++) {
    int indx;
    indx = rand_distr_indx(age_dst, 0, age_sz-1);
    age[i] = age_start[indx] + COV_rand()*10;
    age_distrib[(int)floor(age[i]/5)]++;
  }
  free(age_dst);


  /* Print age distribution */
  int age_dist_test[9]={0};
  for (i=0; i<population; i++) {
    age_dist_test[(int)floor(age[i]/10)]++;
  }

  fprintf(stats, "Testing age distribution\n");
  for (i=0; i<9; i++) {
    fprintf(stats, "age %i fraction %f actual %f \n", i, (double)age_dist_test[i]/population, age_dist[i]) ;
  }
  fprintf(stats, "\n\n");
  fflush(stats);
}
