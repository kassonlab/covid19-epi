/* infections.c
 * Code to manage infection initialization
 * Part of COVID-19 infectious spread model.
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "COV_rand.h"

// Percent per county taken from C19.se infections as of 2020/3/25.
// Initial infections calculated from population admitted to intensive care
// per day from 2020/3/14 to 2020/3/24.
static const double initial_per_county[21] =
{ 0.4234, 0.0404, 0.0336, 0.0843, 0.0257, 0.0079, 0.0071, 0.0020, 0.00475,
  0.0973, 0.0261, 0.1088, 0.0178, 0.0230, 0.0115, 0.0158, 0.0127,  0.0075,
  0.0233, 0.0131, 0.0139 };

/***** initialization array, based on ICU numbers ending on 31 March 2020, day
 0 is 3/21 ******/
static const int initialize_base[15] =
{ 3955, 4068, 5198, 3955, 3616, 4633, 4633, 4859, 5085, 3051, 1921, 2712,
  1469, 1695, 452 };

static void initialize_infections(int    *initial_infections,
                                  double *tau,
                                  int    *infected,
                                  int    *severe,
                                  int    *symptomatic,
                                  int    *county,
                                  int    *num_infect,
                                  int     num_counties,
                                  double  symptomatic_per,
                                  int     population,
                                  double  dt,
                                  double  t,
                                  double *lat_locale,
                                  double *lon_locale,
                                  int    *num_infect_county,
                                  int    *num_infect_age,
                                  double *age,
                                  int   **county_p,
                                  int    *county_size,
                                  int    *locale_HH,
                                  int    *HH,
                                  double *tmp_lat,
                                  double *tmp_lon);

static int read_file_ints(char *fname, int* res) {
  /* reads ints from a file, allocates memory, returns # read */
  int ctr, tmp;
  FILE *infile = fopen(fname, "r");
  if (infile == NULL) {
    fprintf(stderr, "Error opening file %s\n", fname);
    exit(-1);
  }
  for (ctr=0; fscanf(infile, "%d", &tmp) > 0; ctr++) {
    ;
  }
  rewind(infile);
  res = (int *)calloc(ctr, sizeof(int));
  for (int i=0; i < ctr; i++){
    fscanf(infile, "%d", &(res[i]));
  }
  fclose(infile);
  return ctr;
}

void place_initial_infections(char   *initial_infect_filename,
                              char   *initial_immune_filename,
                              FILE   *stats,
                              int    *initial_infections,
                              double *tau,
                              int    *infected,
                              int    *immune,
                              int    *severe,
                              int    *symptomatic,
                              int    *county,
                              int    *num_infect,
                              int     num_counties,
                              double  symptomatic_per,
                              int     population,
                              double  dt,
                              double *lat_locale,
                              double *lon_locale,
                              int    *num_infect_county,
                              int    *num_infect_age,
                              double *age,
                              int   **county_p,
                              int    *county_size,
                              int    *locale_HH,
                              int    *HH) {
  // Place and initialize infections
  // Infections are randomly placed based on number of initial infections
  // Includes infections from t=-11 to t=-1.
  int   *initialize = NULL;
  int   *init_immune = NULL;
  double tmp_t;
  int    i, j;
  double *tmp_lat_v = (double *)calloc(population, sizeof(double));
  double *tmp_lon_v = (double *)calloc(population, sizeof(double));
  int infect_len = 15;
  int immune_len = 15;

  if (initial_infect_filename != NULL) {
    infect_len = read_file_ints(initial_infect_filename, initialize);
    assert((infect_len == 15) || (infect_len == 15*21));
    /* initial infections should either be provided by time or by time & county. */
  } else {
    initialize = (int *)initialize_base;
  }
  if (initial_immune_filename != NULL) {
    immune_len = read_file_ints(initial_immune_filename, init_immune);
    assert((immune_len == 21));
    /* initial immune pop should equal to number of counties. */
  }
  fprintf(stats, "Initial Infections by county \n");
  fflush(stats);

  for (tmp_t = -14; tmp_t <= 0; tmp_t++) {
    int l = -tmp_t;
 
    for (j = 0; j < 21; j++) {
      if (infect_len == 15) {
        /* Initial infections provided as total per time. Distribute by county. */
        initial_infections[j] = initial_per_county[j] * initialize[l];
      } else {
        /* Initial infections provided as total per county per time. */
        initial_infections[j] = initialize[l*21 + j];
      }
      fprintf(stats,
              "time %f county %i initial_infections %i fraction_of_infections %f total_intialized %i\n",
              tmp_t,
              j,
              initial_infections[j],
              initial_per_county[j],
              initialize[l]);
    }
    fprintf(stats, "\n\n");
    fflush(stats);
  
    /* Randomly assign initial infections */
    initialize_infections(initial_infections,
                          tau,
                          infected,
                          severe,
                          symptomatic,
                          county,
                          num_infect,
                          num_counties,
                          symptomatic_per,
                          population,
                          dt,
                          tmp_t,
                          lat_locale,
                          lon_locale,
                          num_infect_county,
                          num_infect_age,
                          age,
                          county_p,
                          county_size,
                          locale_HH,
                          HH,
                          tmp_lat_v,
                          tmp_lon_v);
  }
  /* initialization of immune is much simpler.
   * Algorithm: pick a random person in county. If infected, draw again.
   */
  if (init_immune != NULL) {
    for (i = 0; i < num_counties; i++) {
      if (init_immune[i] > county_size[i]) {
        fprintf(stderr, "Warning: init immune %d larger than county size %d\n",
                init_immune[i], county_size[i]);
        init_immune[i] = county_size[i];
      }
      for (j = 0; j < init_immune[i]; j++) {
        int tmp_immune = (int)(COV_rand() * county_size[i]);
        while (infected[county_p[i][tmp_immune]] > 0) {
          /* retry as long as person is infected */
          tmp_immune = (int)(COV_rand() * county_size[i]);
        }
        immune[county_p[i][tmp_immune]] = 1;
      }
    }
    free(init_immune);
  }

  if (initial_infect_filename != NULL) {
    free(initialize);
  }
  free(tmp_lat_v);
  free(tmp_lon_v);
  fflush(stats);
  printf("All infections initialized\n");
  fflush(stdout);
}

void initialize_infections(int    *initial_infections,
                           double *tau,
                           int    *infected,
                           int    *severe,
                           int    *symptomatic,
                           int    *county,
                           int    *num_infect,
                           int     num_counties,
                           double  symptomatic_per,
                           int     population,
                           double  dt,
                           double  t,
                           double *lat_locale,
                           double *lon_locale,
                           int    *num_infect_county,
                           int    *num_infect_age,
                           double *age,
                           int   **county_p,
                           int    *county_size,
                           int    *locale_HH,
                           int    *HH,
                           double *tmp_lat,
                           double *tmp_lon) {
  int person_infected = 0, county_person_inf = 0;
  int tmp_infect      = 0;
  double min_diff     = 1000;
  double diff_lat_lon = 10;
  int    i, j;
  
  for (i = 0; i < num_counties; i++) {
    tmp_infect = 0;
    
    while ((tmp_infect < initial_infections[i])) {
      int tmp_j = 0;
      
      // Test up to population to see if we can find someone who fits a
      // perviously determined cluster.  If not, leave this loop and pick a
      // random person.
      county_person_inf = (int)(COV_rand() * county_size[i]);
      person_infected   = county_p[i][county_person_inf];
      
      while (((county[person_infected] != i) ||
              (infected[person_infected] !=
               0) || min_diff > 0.5) && tmp_j < county_size[i]) {
                /* pick first available person in the county to infect if the randomly
                 choosen aren't free to pick */
                county_person_inf++;
                
                if (county_person_inf >= county_size[i]) {
                  county_person_inf = 0;
                }
                person_infected = county_p[i][county_person_inf];
                
                if (county[person_infected] != i) {
                  fprintf(stderr,
                          "Error: random person is not in expected county %d, is in %d\n",
                          i,
                          county[person_infected]);
                }
                min_diff = 1000;
                
                if ((t < -13) || (num_infect_county[i] == 0)) {
                  min_diff = 0;
                } else {
                  for (j = 0; j < *num_infect; j++) {
                    diff_lat_lon =
                    (fabsf(tmp_lat[j] -
                           lat_locale[locale_HH[HH[person_infected]]]) +
                     fabsf(tmp_lon[j] - lon_locale[locale_HH[HH[person_infected]]]));
                    
                    if (diff_lat_lon < min_diff) {
                      min_diff = diff_lat_lon;
                    }
                  }
                }
                tmp_j++;
              }
      
      if (min_diff > 1) {
        while ((county[person_infected] != i) ||
               (infected[person_infected] != 0)) {
          county_person_inf++;
          
          if (county_person_inf >= county_size[i]) {
            county_person_inf = 0;
          }
          person_infected = county_p[i][county_person_inf];
          
          if (county[person_infected] != i) {
            fprintf(stderr,
                    "Error: random person is not in expected county %d, is in %d\n",
                    i,
                    county[person_infected]);
          }
        }
      }
      
      infected[person_infected] = 1;
      severe[person_infected]   = (int)(COV_rand() * 2);
      
      if (COV_rand() < symptomatic_per) {
        symptomatic[person_infected] = 1;
      }
      tau[person_infected] = t;
      tmp_lat[*num_infect] = lat_locale[locale_HH[HH[person_infected]]];
      tmp_lon[*num_infect] = lon_locale[locale_HH[HH[person_infected]]];
      *num_infect          = *num_infect + 1;
      num_infect_county[i]++;
      num_infect_age[(int)(floor(age[person_infected] / 5))]++;
      tmp_infect++;
    }
  }
}

