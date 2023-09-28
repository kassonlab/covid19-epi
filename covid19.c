/* COVID-19 infectious spread model
 * Original code 2020 by Jasmine Gardner.
 * Subsequent modifications by development team.
 * Hosted at https://github.com/kassonlab/covid19-epi
 * Original publication DOI 10.1101/2020.04.11.20062133
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


#ifdef _OPENMP
# include <omp.h>
#endif /* ifdef _OPENMP */

#include "common.h"
#include "interventions.h"
#include "infections.h"

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

/* supported population range */
#define POP_MIN 100000
#define POP_MAX 10098554

int typical_max_HH_sz = 7;

int sort_household_members = 0;

int full_kappa = 0;

/* Beta scaling factors balance the number of potentially infectious contacts
   that each infected individual has between household, workplace, and community.
   A uniform prior would state that 1/3 of the contacts are made in each of these
   contexts.  This is also the assumption used in doi:10.25561/77482.  Base values
   for beta parameters were obtained by Ferguson and colleagues based on SE Asia.
   To scale these parameters, beta scaling factors are applied to maintain the
   balance of contacts.  Work by Mossong et al., 2008 (doi:10.1371/journal.pmed.0050074)
   shows that the balance of contacts between home, workplace, and community is
   a) relatively similar to the uniform prior and b) most importantly, does 
   not significantly vary across 8 European countries despite population changing
   by a factor of >10 and average household size changing from ~2 to >2.8.
   These direct epidemiological data thus support the beta scaling approach used here.
   The concept applied is that each infected individual will have N contacts in a day;
   each of these contacts has a probability of infection calculated using the
   equations for household, workplace, and community transmission.  Because those equations
   are summed over all individual-individual pairs sharing household, workplace,
   or community, the beta scaling factors are required for normalization of the contact
   probabilities.  The R0_scaling parameter is then the one fitted parameter,
   procedure discussed and data given in the supplement to doi.org:10.1093/cid/ciaa864,
   that covers infectious probability of those contacts.

   Overall population is 8.4X smaller in Sweden, attributed as a 8.4X increase
       in community transmission factor.
    Household size is approximately 2X smaller in Sweden, 2 vs. 4, attributed as
       a 2X increase in household transmission factor.
    R0 Scaling factor used to set approximate doubling time.  R0_scale=2.2
       attributes to a doubling time of approximately 3 days while R0_scale=1.8
       attributes to a doubling time of approximately 5 days.
 */
double betac_scale         = 8.4, betah_scale = 2.0, betaw_scale = 1.0,
       R0_scale            = 1.8;
double population_fraction = 0; /* Fraction for adjustable interventions. */

#define pi 3.14159265358979323846

/* Adding in a method for considering demographic data */
void dem_data(double *dem, int population, FILE *dem_info, int *dem_distrib){
  int i; // counter
  int ret, dem_sz;
  
  int dem_start[] = {0, 1, 2, 3, 4 , 5};
  double dem_dist[] = {0.05, 0.1, 0.3, 0.45, 0.8, 0.9} //PUT NUMBERS HERE WHICH CORRESPOND WITH THE DISTRIBUTION OF THESE CATEGORIES THIS NEEDS TO BE A CDF TYPE OF SITUATION
  double *dem_dst;

  dem_sz = (int)(sizeof(dem_dist) / sizeof(dem_dist[0]));
  ret = generate_inc_distr_vec(&dem_dst, dem_dist, dem_sz, "demographic");

  if (ret) {
    fprintf(stderr, "Bailing out on demographic initialization\n");
    exit(1);
  }

  for (i = 0; i < population; i++) {
    // Generate a random number between 0 and 1
    double rand_num = (double)rand() / RAND_MAX;

    int category = 0;
      while (category < dem_sz && rand_num > dem_dist[category]){
        category++;
      }

      dem[i] = dem_start[category]

  free(dem_dst);
}
    

void age_dist(double *age, int population, FILE *stats, int *age_distrib) {
  int i; // Counter
  int ret, age_sz;

  int age_start[]   = { 0, 10, 20, 30, 40, 50, 60, 70, 80 };
  double age_dist[] =
  { 0.10721, 0.11553, 0.12439, 0.13475, 0.12578, 0.12704, 0.10782, 0.09952,
    0.05796 };
  double *age_dst;

  age_sz = (int)(sizeof(age_dist) / sizeof(age_dist[0]));
  ret    = generate_inc_distr_vec(&age_dst, age_dist, age_sz, "age");

  if (ret) {
    fprintf(stderr, "Bailing out on age initialization\n");
    exit(1);
  }

  for (i = 0; i < population; i++) {
    int indx;
    indx   = rand_distr_indx(age_dst, 0, age_sz - 1);
    age[i] = age_start[indx] + COV_rand() * 10;
    age_distrib[(int)floor(age[i] / 5)]++;
  }
  free(age_dst);


  /* Print age distribution */
  int age_dist_test[9] = { 0 };

  for (i = 0; i < population; i++) {
    age_dist_test[(int)floor(age[i] / 10)]++;
  }

    fprintf(stats, "Testing age distribution\n");

  for (i = 0; i < 9; i++) {
    fprintf(stats,
            "age %i fraction %f actual %f \n",
            i,
            (double)age_dist_test[i] / population,
            age_dist[i]);
  }
    fprintf(stats, "\n\n");
  fflush(stats);
}

int cmp_int(void const *a_, void const *b_) {
  int a = *(int *)a_;
  int b = *(int *)b_;

  if (a != b) {
    return (a < b) ? -1 : 1;
  }
  return 0;
}

/* Puts households in certain locations based on population density
 * distribution.  Only really correct for full population but works for
 * smaller populations.  Biases towards smaller households for smaller
 * populations.  Each household is also fed into a locality based on the
 * shortest distance to the center of that locality for purposes of
 * school and workplace choice.
 */
void household_lat_long(int     num_households,
                        int    *HH,
                        double *lat_city,
                        double *long_city,
                        int     num_cities,
                        int    *city,
                        int    *county,
                        int    *city_county,
                        int    *city_size,
                        int    *county_size,
                        int     population,
                        double *age,
                        int   **per_HH_members,
                        char  **county_names,
                        double *county_pop,
                        double  tot_pop_actual,
                        FILE   *stats,
                        int   **county_p,
                        int    *num_locale,
                        double *lat_locale,
                        double *lon_locale,
                        double *pop_density_init_num,
                        int   **locale_to_HH,
                        int    *locale_HH,
                        double  land_pop_total_density) {
  /* Initialize helpers for population density */
  double tot_pop_density = 0;

  /* Initialize helpers for household distribution */
  int *city_HH   = malloc(num_households * sizeof(int));
  int *county_HH = malloc(num_households * sizeof(int));
  int  county_num;
  int  city_num;

  double dist1;
  double min_dist = 1000000;
  int    placement;
  int    i, j;

  int tmp_county_count[21]   = { 0 };
  int tmp_county_density[21] = { 0 };

  /* Fill up households. */
  int HH_count    = 0; // Counter
  int HH_person   = 0; // Counter
  int max_HH_size = 0;


  /* Initially set all households to -1 */
  for (i = 0; i < population; i++)
  {
    HH[i] = -1;
  }

  /* Save list of households in each locale. */
  int   original_num_locale = *num_locale;
  int **county_list;
  county_list = (int **)calloc(*num_locale, sizeof(int *));
  int *locale_count;
  locale_count = (int *)calloc(*num_locale, sizeof(int));
  int *locale_HH_count;
  locale_HH_count = (int *)calloc(*num_locale, sizeof(int));


  city_num = -1;

  while ((HH_count < num_households) && (HH_count < *num_locale))
  {
    double tmp_lat, tmp_lon;
    min_dist = 1000000;

    tmp_lat = lat_locale[HH_count];
    tmp_lon = lon_locale[HH_count];

    /* Scale population density */
    pop_density_init_num[HH_count] = ceil(
      pop_density_init_num[HH_count] * population / land_pop_total_density);
    tot_pop_density               += pop_density_init_num[HH_count];

    if (HH_count + 1 < *num_locale)
    {
      // Conditional added to prevent read overrun
      if ((tot_pop_density > population) ||
          (ceil(pop_density_init_num[HH_count +
                                     1] * population / land_pop_total_density) <
           0.5))
      {
        *num_locale = HH_count + 1;
      }
    }

    // Determine city of each population square.  Use city data to determine
    // which schools students attend.  Workplaces are placed by county. //
    for (j = 0; j < num_cities; j++)
    {
      dist1 = distance(tmp_lat, tmp_lon, lat_city[j], long_city[j], 'K');

      if (dist1 < min_dist)
      {
        min_dist = dist1;
        city_num = j;
      }
    }

    if (city_num < 0)
    {
      fprintf(stderr, "Error in household_lat_long: city_num < 0\n");
      exit(0);
    }

    county_num = city_county[city_num];
    tmp_county_count[county_num]++;
    tmp_county_density[county_num] += pop_density_init_num[HH_count];

    /* Set up household */
    city_HH[HH_count]   = city_num;
    county_HH[HH_count] = county_num;
    arrput(county_list[HH_count],  HH_count);
    locale_HH[HH_count] = HH_count;
    arrput(locale_to_HH[HH_count], HH_count);
    locale_HH_count[HH_count] += 1;

    /* rooms_HH[HH_count] = number of rooms */



    /* Allocate an adult to the household */
    while (age[HH_person] < 20)
    {
      HH_person++;
    }
    HH[HH_person]     = HH_count;
    city[HH_person]   = city_HH[HH_count];
    county[HH_person] = county_HH[HH_count];
    county_size[county[HH_person]]++;
    arrput(county_p[county[HH_person]],   HH_person);
    city_size[city[HH_person]]++;
    arrput(per_HH_members[HH[HH_person]], HH_person);
    locale_count[HH_count] += 1;

    if (arrlen(per_HH_members[HH[HH_person]]) > max_HH_size)
    {
      max_HH_size = arrlen(per_HH_members[HH[HH_person]]);
    }
    HH_person++;
    HH_count++;
  }

  /* Uncomment to test population denity WRT county.
      // Test population based on method vs actual population density on county
         level.  Lines up fairly well. //
      for (i=0;i<21; i++) {
          printf("counties %i %s number of locales %i calculated density %f
             actual density %f\n", i, county_names[i], tmp_county_count[i],
             tmp_county_density[i]/(double)tot_pop_density,
             county_pop[i]/tot_pop_actual);
          fflush(stdout);
      }
   */


  placement = 0; // Keeps track of household placement after first loop of
                 // locales.

  /* First add adult head of household to the rest of the households. */
  while (HH_count < num_households)
  {
    while (age[HH_person] < 20)
    {
      HH_person++;
    }

    /* Place another household in locales with population density greater than
      2.2*(households in locale). */
    if (locale_HH_count[placement] * 2.2 > pop_density_init_num[placement])
    {
      placement++; /* continue with next locale */
    }
    city_HH[HH_count]   = city_HH[placement];
    county_HH[HH_count] = county_HH[placement];
    arrput(county_list[placement],  HH_count);
    locale_HH[HH_count] = placement;
    arrput(locale_to_HH[placement], HH_count);
    locale_HH_count[placement] += 1;

    /* Set up head of household. */
    HH[HH_person]     = HH_count;
    city[HH_person]   = city_HH[HH_count];
    county[HH_person] = county_HH[HH_count];
    county_size[county[HH_person]]++;
    arrput(county_p[county[HH_person]],   HH_person);
    city_size[city[HH_person]]++;
    locale_count[placement] += 1;
    arrput(per_HH_members[HH[HH_person]], HH_person);

    if (arrlen(per_HH_members[HH[HH_person]]) > max_HH_size)
    {
      max_HH_size = arrlen(per_HH_members[HH[HH_person]]);
    }
    HH_person++;
    HH_count++;
  }

  placement = 0;

  /* Distribute remaining people randomly.  This could be changed to a
    distribution to more realistically reflect household size in the future. */
  for (HH_person = 0; HH_person < population; HH_person++)
  {
    int tmp_HH, tmp_county_HH;

    if (HH[HH_person] == -1)
    {
      /* Place people in random households within locale until locale is full.
        */
      if (locale_count[placement] + 1  > pop_density_init_num[placement])
      {
        placement++; /* continue with next locale */
      }

      /* Pick a random household in the locale. */
      int only_once;
      only_once     = 0;
      tmp_county_HH = (int)(COV_rand() * locale_HH_count[placement]);
      tmp_HH        = county_list[placement][tmp_county_HH];

      while (arrlen(per_HH_members[tmp_HH]) + 1 > typical_max_HH_sz)
      {
        tmp_county_HH++;

        if (tmp_county_HH >= locale_HH_count[placement])
        {
          if (only_once)
          {
            if (placement + 1 >= *num_locale)
            {
              printf("Bailing out, no more locale to put people in\n");
              exit(0);
            }
            placement++;
            tmp_county_HH = (int)(COV_rand() * locale_HH_count[placement]);
            only_once     = 0;
          }
          else
          {
            /* Start over from the beginning of county_list[placement] in the
              search for a unfilled household */
            /* but only do that once */
            tmp_county_HH = 0;
            only_once     = 1;
          }
        }
        tmp_HH = county_list[placement][tmp_county_HH];
      }
      HH[HH_person] = tmp_HH;
      fflush(stdout);
      city[HH_person]   = city_HH[HH[HH_person]];
      county[HH_person] = county_HH[HH[HH_person]];
      county_size[county[HH_person]]++;
      arrput(county_p[county[HH_person]],   HH_person);
      city_size[city[HH_person]]++;
      locale_count[placement] += 1;
      arrput(per_HH_members[HH[HH_person]], HH_person);

      if (arrlen(per_HH_members[HH[HH_person]]) > max_HH_size)
      {
        max_HH_size = arrlen(per_HH_members[HH[HH_person]]);
      }
    }
  }

  printf("max_HH_size = %d\n", max_HH_size);
  fflush(stdout);

  if (sort_household_members)             /*sort_household_members only set to 0/1*/
  {
    for (size_t hh = 0; hh < num_households; ++hh)
    {
      int *set = per_HH_members[hh];

      if (set)
      {
        qsort(set, arrlen(set), sizeof(int), cmp_int);
      }
    }
  }


  /* Print household distribution */
  int *HH_dist_test;
  HH_dist_test = (int *)calloc(max_HH_size + 1, sizeof(int));
  int *HH_county_test;
  HH_county_test = (int *)calloc(21, sizeof(int));

  for (i = 0; i < num_households; i++)
  {
    HH_dist_test[arrlen(per_HH_members[i])]++;
    HH_county_test[county_HH[i]]++;
  }

    fprintf(stats, "Household distributions \n");

  for (i = 0; i <= max_HH_size; i++)
  {
    fprintf(stats,
            "household_size %i fraction_of_households %f num_households %i total_households %i \n",
            i,
            HH_dist_test[i] / (double)num_households,
            HH_dist_test[i],
            num_households);
    fflush(stats);
  }
  free(HH_dist_test);

    fprintf(stats, "\n\n County distribution \n");

  for (i = 0; i < 21; i++)
  {
    fprintf(stats,
            "%s county %i population %i fraction %f actual %f num_households %i \n",
            county_names[i],
            i,
            county_size[i],
            county_size[i] / (double)population,
            county_pop[i] / tot_pop_actual,
            HH_county_test[i]);
    fflush(stats);
  }
  fprintf(stats, "\n\n");
  fflush(stats);
  free(HH_county_test);

  for (i = 0; i < original_num_locale; i++)
  {
    arrfree(county_list[i]);
  }
  free(county_list);
  free(locale_count);
  free(locale_HH_count);

  free(city_HH);
  free(county_HH);
}

void city_lat_long(int    *num_cities,
                   double *lat_city,
                   double *long_city,
                   char  **cities,
                   int    *county,
                   char  **county_names,
                   int     num_county) {
  /* Get longitude and latitude of cities in Sweden from CSV */
  const char *cityfile = "cities_all.csv";
  FILE *fp             = fopen(cityfile, "r"); // Not sure about the validity of
                                               // this file.  Could use a better
                                               // source.

  /* Counters for parsing file. */
  char   tmp[199];                             // municipality name
  char   tmp1[199];                            // locality name
  char   tmp2[199];                            // county name
  double tmp_lat;                              // tmp latitude
  double tmp_lon;                              // tmp longitude
  int    i, j;

  if (fp == NULL)
  {
    fprintf(stderr, "Error opening file %s.  Exiting.\n", cityfile);
    exit(0);
  }

  for (i = 0; i < 2000; i++) {
    int got = fscanf(fp,
                     "%[^,\n]%*c%[^,\n]%*c%[^,\n]%*c%lf%*c%lf\n",
                     tmp1,
                     tmp,
                     tmp2,
                     &tmp_lat,
                     &tmp_lon);

    if (got != 5) break;  // wrong number of tokens - maybe end of file

    /* Save city name, longitude, latitude, and increase number of cities.
      Looking at unique municipalities.*/
    cities[*num_cities]    = tmp1;
    lat_city[*num_cities]  = tmp_lat;
    long_city[*num_cities] = tmp_lon;
    *num_cities            = *num_cities + 1;

    /* Save county that the city is in. */
    for (j = 0; j < num_county; j++) {
      if ((int)(strcmp(county_names[j], tmp2)) == 0) {
        county[i] = j;
        break;
      } else if (j == num_county - 1) {
        printf("no_county_match %s %s %i %i\n",
               county_names[j],
               tmp2,
               j,
               strcmp(county_names[j], tmp2));
      }
    }
  }

  fclose(fp);
}

void job_dist(int    *job_status,
              int   **job_status_city,
              double *age,
              int    *county,
              int    *city,
              int     population,
              int     num_cities,
              int     num_counties,
              int    *county_size,
              FILE   *stats) {
  int i, j;


  /* All currently randomly placed based on county.  Would like to do per town
    but each town needs inhabitants. See commented section below for per city
    distribution of schools */
  for (i = 0; i < population; i++) {
    if ((age[i] < 1) || (age[i] > 75)) {
      job_status[i] = 0;
      job_status_city[0][county[i]]++;
    } else if ((age[i] >= 1) && (age[i] < 3)) {
      if (COV_rand() < 0.7800) {
        job_status[i] = 1;
        job_status_city[1][city[i]]++;
      } else {
        job_status[i] = 0;
        job_status_city[0][city[i]]++;
      }
    } else if ((age[i] >= 3) && (age[i] < 6)) {
      if (COV_rand() < 0.9500) {
        job_status[i] = 1;
        job_status_city[1][city[i]]++;
      } else {
        job_status[i] = 0;
        job_status_city[0][city[i]]++;
      }
    } else if ((age[i] >= 6) && (age[i] < 15)) {
      job_status[i] = 2;
      job_status_city[2][city[i]]++;
    } else if ((age[i] >= 15) && (age[i] < 22)) {
      job_status[i] = 3;
      job_status_city[3][county[i]]++;
    } else if ((age[i] >= 22) && (age[i] <= 65)) {
      if (COV_rand() < 0.773) {
        // 17.25% of workforce is in healthcare from OECD 2017 statstics.
        //  Assume 1/4 of these are in hospitals.
        if (COV_rand() < 0.04325) {
          job_status[i] = 5;
          job_status_city[5][county[i]]++; // Workplace is based on county, not
                                           // city.
        } else {
          job_status[i] = 4;
          job_status_city[4][county[i]]++; // Workplace is based on county, not
                                           // city.
        }
      } else {
        job_status[i] = 0;
        job_status_city[0][county[i]]++;
      }
    } else if ((age[i] >= 65) && (age[i] <= 75)) {
      if (COV_rand() < 0.172) {
        job_status[i] = 4;
        job_status_city[4][county[i]]++; // Workplace is based on county, not
                                         // city.
      } else {
        job_status[i] = 0;
        job_status_city[0][county[i]]++;
      }
    }
  }


  /* Print job distribution */
  int  job_dist_test[6] = { 0 };
  int *unemployed;
  unemployed = (int *)calloc(num_counties, sizeof(int));
  int *working_age;
  working_age = (int *)calloc(num_counties, sizeof(int));
  int **city_dist_test;
  city_dist_test = (int **)calloc(6, sizeof(int *));

  for (i = 0; i < 6; i++) {
    city_dist_test[i] = (int *)calloc(num_counties, sizeof(int));
  }

  for (i = 0; i < population; i++) {
    job_dist_test[job_status[i]]++;
    city_dist_test[job_status[i]][county[i]]++;

    if ((age[i] > 22) && (age[i] <= 65)) {
      if (job_status[i] == 0) {
        unemployed[county[i]]++;
      }
      working_age[county[i]]++;
    }
  }

      fprintf(stats, "Job Distribution \n");

  for (j = 0; j < num_counties; j++) {
    for (i = 0; i < 6; i++) {
      fprintf(stats,
                     "job_status %i county %i fraction_of_jobs_total %f num_jobs_in_county %i num_unemployed %i fraction_unemployed  %f \n",
              i,
              j,
              job_dist_test[i] / (double)population,
              city_dist_test[i][j],
              unemployed[j],
              (double)unemployed[j] / working_age[j]);
    }
      fprintf(stats, "\n");
  }
  fflush(stats);

    free(unemployed);
    free(working_age);

  for (i = 0; i < 6; i++) {
    free(city_dist_test[i]);
  }
    free(city_dist_test);
}

void workplace_dist(int  *workplace,
                    int  *job_status,
                    int **job_status_county,
                    int  *city,
                    int   num_cities,
                    int  *county,
                    int   num_counties,
                    int   population,
                    int  *max_num_WP,
                    int  *hosp_num,
                    FILE *stats) {
  int pp_class     = 19;  // Assumption of 15 children per class. Not used yet
  int pp_preschool = 53;  // Assumption of 200 children per school.
  int pp_school    = 220; // Assumption of 200 children per school.
  int pp_hospital  = 120; // Assumption of 120 people per hospital.
  int pp_work      = 15;  // Assumption of 15 people per close work group.
  int i;
  int j;
  int num_workplaces[6][num_cities];

  memset(num_workplaces,  0, 6 * num_cities * sizeof(int));
  int num_workplaces2[6];
  memset(num_workplaces2, 0, 6 * sizeof(int));

  /* THIS WILL NOW BE BY CITY FOR SCHOOLS AND PRESCHOOLS. */
  for (i = 0; i < num_cities; i++) {
    for (j = 0; j < 2; j++) {
      /* First only schools then workplaces */
      if (job_status_county[j][i] > 0) {
        num_workplaces[j][i] = ceil(
          job_status_county[j][i] / (double)pp_preschool);
        num_workplaces2[j]  += ceil(
          job_status_county[j][i] / (double)pp_preschool);

        if (num_workplaces2[j] > *max_num_WP) {
          *max_num_WP = num_workplaces2[j];
        }
      }
    }

    for (j = 2; j < 3; j++) {
      /* First only schools then workplaces */
      if (job_status_county[j][i] > 0) {
        num_workplaces[j][i] =
          ceil(job_status_county[j][i] / (double)pp_school);
        num_workplaces2[j] +=
          ceil(job_status_county[j][i] / (double)pp_school);

        if (num_workplaces2[j] > *max_num_WP) {
          *max_num_WP = num_workplaces2[j];
        }
      }
    }
  }

  // Broken down in case we want to do schools by municipality.
  for (i = 0; i < num_counties; i++) {
    for (j = 3; j < 6; j++) {
      if ((job_status_county[j][i] > 0) && (j == 4)) {
        num_workplaces[j][i] = ceil(job_status_county[j][i] / (double)pp_work);
        num_workplaces2[j]  += ceil(job_status_county[j][i] / (double)pp_work);

        if (num_workplaces2[j] > *max_num_WP) {
          *max_num_WP = num_workplaces2[j];
        }
      } else if ((job_status_county[j][i] > 0) && (j == 3)) {
        num_workplaces[j][i] =
          ceil(job_status_county[j][i] / (double)pp_school);
        num_workplaces2[j] +=
          ceil(job_status_county[j][i] / (double)pp_school);

        if (num_workplaces2[j] > *max_num_WP) {
          *max_num_WP = num_workplaces2[j];
        }
      } else if ((job_status_county[j][i] > 0) && (j == 5)) {
        num_workplaces[j][i] = ceil(
          job_status_county[j][i] / (double)pp_hospital);
        num_workplaces2[j]  += ceil(
          job_status_county[j][i] / (double)pp_hospital);

        if (num_workplaces2[j] > *max_num_WP) {
          *max_num_WP = num_workplaces2[j];
        }
      }
    }

    // Need to use floor+1 equation because each county should have a hospital
    // even if no one works there. */
    num_workplaces[5][i] =
      floor(job_status_county[5][i] / (double)pp_hospital) + 1;
    num_workplaces2[5] +=
      floor(job_status_county[5][i] / (double)pp_hospital) + 1;

    if (num_workplaces2[5] > *max_num_WP) {
      *max_num_WP = num_workplaces2[5];
    }
    hosp_num[i] = num_workplaces[5][i];
  }

  for (i = 0; i < population; i++) {
    /* Try to minimize necessary memory by making job_numbers independent with
      job_status. */
    int prior_workplaces = 0;

    if ((job_status[i] < 3) && (num_workplaces[job_status[i]][city[i]] > 0)) {
      for (j = 0; j < city[i]; j++) {
        prior_workplaces += num_workplaces[job_status[i]][j];
      }
      workplace[i] =
        (int)(COV_rand() *
              (num_workplaces[job_status[i]][city[i]])) + prior_workplaces;
    } else if ((job_status[i] > 3) &&
               ((num_workplaces[job_status[i]][county[i]]) > 0)) {
      for (j = 0; j < county[i]; j++) {
        prior_workplaces += num_workplaces[job_status[i]][j];
      }
      workplace[i] =
        (int)(COV_rand() *
              (num_workplaces[job_status[i]][county[i]])) + prior_workplaces;
    } else if ((num_workplaces[job_status[i]][county[i]]) > 0) {
      for (j = 0; j < county[i]; j++) {
        prior_workplaces += num_workplaces[job_status[i]][j];
      }
      workplace[i] =
        (int)(COV_rand() *
              (num_workplaces[job_status[i]][county[i]])) + prior_workplaces;
    } else {
      workplace[i] = 0;
      printf("no workplaces %i %i %i \n", i, job_status[i], county[i]);
    }
  }
}

double calc_kappa(double  t,
                  double  tau,
                  int     symptomatic,
                  double  dt,
                  double *kappa_vals,
                  int     hosp,
                  int     icu) {
  double kappa;
  double t1;
  int    t2;

  /* Determine kappa for infected person.  This is the
   * infectiousness of the person based on time since infection
   * started.  Latency period is 4.6 days.  Infection starts at
   * 5.1 days and lasts for 6 days.  Sympotmatic people are twice
   * as likely to infect others as asymptomatic.  Kappa is a log
   * normal function with mean of -0.72 and standard deviation of
   * 1.8.  From Ferguson Nature 2005
   */
  if (t - tau <= 4.6) {
    kappa = 0.;
  } else if ((t - tau > 11.1) && (hosp == 0) && (icu == 0)) {
    kappa = 0.; // # Recovered or dead
  } else {
    /* First 2 lines calculates kappa on the fly, second two get precalculated
      kappa from array. */
    if (full_kappa) {
      t1    = (log(t - tau - 4.6) + 0.72) / 1.8;
      kappa = exp(-0.5 * pow(t1, 2)) / ((t - tau - 4.6) * 1.8 * sqrt(2 * M_PI));
    } else {
      t2    = (t - tau) / dt;
      kappa = kappa_vals[t2];
    }
  }

  if (symptomatic == 0) {
    kappa = kappa * 0.5;
  }
  return kappa * R0_scale;
}

void segment_population(int    *num_sus,
                        int    *num_infectious,
                        int    *num_hosp,
                        int    *num_icu,
                        int    *num_sus_HCW,
                        int    *num_infectious_HCW,
                        int    *num_hosp_HCW,
                        int    *num_icu_HCW,
                        int    *infected,
			int    *immune,
                        int    *infectious,
                        int    *sus_list,
                        int    *hosp_list,
                        int    *hosp_pop,
                        int    *icu_pop,
                        int    *icu_list,
                        double *tau,
                        int     population,
                        double  t,
                        int    *num_sus_county,
                        int    *num_infectious_county,
                        int    *num_infected_county,
                        int    *num_hosp_county,
                        int    *num_icu_county,
                        int    *num_sus_age,
                        int    *num_infectious_age,
                        int    *num_infected_age,
                        int    *num_hosp_age,
                        int    *num_icu_age,
                        double *age,
                        int    *county,
                        int     print_loc,
                        double *lat_locale,
                        double *lon_locale,
                        FILE   *lat_lon_out,
                        int    *locale_HH,
                        int    *HH,
                        int    *job_status,
                        int    *recovered,
                        int    *dead) {
  int i;

  *num_sus            = 0;
  *num_infectious     = 0;
  *num_hosp           = 0;
  *num_icu            = 0;
  *num_sus_HCW        = 0;
  *num_infectious_HCW = 0;
  *num_hosp_HCW       = 0;
  *num_icu_HCW        = 0;

  for (i = 0; i < population; i++) {
    if ((infected[i] == 0) && (immune[i] == 0)){
      sus_list[*num_sus] = i;
      *num_sus           = *num_sus + 1;
      num_sus_county[county[i]]++;
      num_sus_age[(int)floor(age[i] / 5)]++;

      if (job_status[i] == 5) {
        *num_sus_HCW = *num_sus_HCW + 1;
      }
    } else if ((tau[i] < t - 4.6) && (recovered[i] == 0) && (dead[i] == 0)) {
      infectious[*num_infectious] = i;
      *num_infectious             = *num_infectious + 1;
      num_infectious_county[county[i]]++;
      num_infectious_age[(int)floor(age[i] / 5)]++;

      if (print_loc == 1) {
            fprintf(lat_lon_out,
                "Time %f Person %i lat %f lon %f \n",
                t,
                i,
                lat_locale[locale_HH[HH[i]]],
                lon_locale[locale_HH[HH[i]]]);
      }

      if (job_status[i] == 5) {
        *num_infectious_HCW = *num_infectious_HCW + 1;
      }
    } else {
      if (print_loc == 1) {
            fprintf(lat_lon_out,
                "Time %f Person %i lat %f lon %f \n",
                t,
                i,
                lat_locale[locale_HH[HH[i]]],
                lon_locale[locale_HH[HH[i]]]);
      }
    }


    if (icu_pop[i] == 1) {
      icu_list[*num_icu] = i;
      *num_icu           = *num_icu + 1;
      num_icu_county[county[i]]++;
      num_icu_age[(int)floor(age[i] / 5)]++;

      if (job_status[i] == 5) {
        *num_icu_HCW = *num_icu_HCW + 1;
      }
    }

    if (hosp_pop[i] > 0) {
      hosp_list[*num_hosp] = i;
      *num_hosp            = *num_hosp + 1;
      num_hosp_county[county[i]]++;
      num_hosp_age[(int)floor(age[i] / 5)]++;

      if (job_status[i] == 5) {
        *num_hosp_HCW = *num_hosp_HCW + 1;
      }
    }
  }
            fprintf(lat_lon_out, "\n\n");
}

double calc_household_infect(double kappa,
                             double omega,
                             double rho,
                             int    HH_size,
                             double alpha,
                             int    severe) {
  double betah = 0.55; // Scaled from betah=0.4 in influenza pandemic with
                       // R0=1.6, COVID-19 R0=2.2 (Ferguson 2020)

  /* potentially redefine beta as a function of overcrowding */


  return betah_scale * betah * kappa * rho *
         (1 + 0 * (omega - 1)) / (pow((double)HH_size, alpha));
}

double calc_workplace_infect(int    job_status,
                             double kappa,
                             double omega,
                             double rho,
                             int    workplace_size,
                             int    severe,
                             double Iw) {
  double betap[] = { 0.0, 1.1, 1.1, 1.1, 0.55, 0.1375 }; // Spread in all types
                                                         // of schools
                                                         // (preschool to
                                                         // college) is twice
                                                         // that of workplace
  double psi[]   = { 0.0, 0.1, 0.2, 0.25, 0.5, 0.5 };    // Accounts for
                                                         // absenteeism based on
                                                         // severe illness.
                                                         // Ferguson Nature 2006

  return betaw_scale * Iw * betap[job_status] * kappa * rho * 
         (1 + (double)severe *
          (omega * psi[job_status] - 1)) / ((double)workplace_size);
}

double calc_community_infect(double kappa, double omega, double rho, int severe, double d) {
  /* need to work on this.  Perhaps we take a random distance for each two
    people based on population density, number of people in county, county area,
    etc. */
  double zeta[] =
  { 0.1,    0.25,  0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.75,
    0.50,   0.25, 0.25, 0.25 };   //   # Travel related parameter for community
                                  // transmission. Ferguson Nature 2006
  double fd;
  double betac = 0.103;           // Scaled from betac=0.075 in influenza
                                  // pandemic with R0=1.6, COVID-19 R0=2.2
                                  // (Ferguson 2020)

  fd = 1 / (1 + pow((d / 4), 3)); // kernel density function as parameterized
                                  // for GB.
  return betac_scale * betac * kappa * fd * rho * (1 + 0 * (omega - 1));
}

void hosp_entry(double  t,
                int     num_infectious,
                int    *infectious,
                double *age,
                int    *icu_pop,
                int    *hosp_pop,
                int    *symptomatic,
                double *tau,
                int    *workplace_tmp,
                int    *workplace_num,
                int    *county,
                double  dt) {
  double hosp[] =
  { 0.001, 0.003, 0.012, 0.032, 0.049, 0.102, 0.166, 0.243, 0.273 };            //
                                                                                // #
                                                                                // From
                                                                                // Ferguson
                                                                                // 2020
  double icu[] = { 0.05, 0.05, 0.05, 0.05, 0.063, 0.122, 0.274, 0.432, 0.709 }; //
                                                                                // #
                                                                                // percent
                                                                                // of
                                                                                // hospitalized
                                                                                // that
                                                                                // need
                                                                                // icu
                                                                                // from
                                                                                // Ferguson
                                                                                // 2020.
  int    age_group;
  int    infec_person;
  int    i;

  for (i = 0; i < num_infectious; i++) {
    infec_person = infectious[i];

    if ((tau[infec_person] >= t - 10) && (tau[infec_person] < t - 10 + dt) &&
        (symptomatic[infec_person] == 1)) {
      age_group = floor(age[infec_person] / 10);

      if (COV_rand() < hosp[age_group]) {
        if (COV_rand() < icu[age_group]) {
          icu_pop[infec_person]  = 1;
          hosp_pop[infec_person] = 2;
        } else {
          hosp_pop[infec_person] = 1;
        }

        /* put them in a hospital */
        workplace_tmp[infec_person] =
          (int)(COV_rand() * workplace_num[county[infec_person]]);
      }
    }
  }
}

void hosp_release(double  t,
                  int     num_hosp,
                  int    *hosp_list,
                  double *tau,
                  int    *recovered,
                  int    *hosp_pop,
                  int    *num_recovered,
                  int    *recovered_hosp,
                  int    *recovered_icu,
                  double  dt,
                  int    *num_recovered_county,
                  int    *num_recovered_age,
                  double *age,
                  int    *county,
                  int    *num_recovered_hosp_county,
                  int    *num_recovered_hosp_age,
                  int    *num_recovered_icu_county,
                  int    *num_recovered_icu_age,
                  int    *num_recovered_HCW,
                  int    *recovered_hosp_HCW,
                  int    *recovered_icu_HCW,
                  int    *job_status) {
  int i;
  int infec_person;

  for (i = 0; i < num_hosp; i++) {
    infec_person = hosp_list[i];

    /* Regular hospital patients get released after 8 days (18 days from
      infection).  ICU patients get released after 26 days but do not go into
      regular hospital until 20 days.  */
    if ((tau[infec_person] >= t - 18) && (tau[infec_person] < t - 18 + dt) &&
        (hosp_pop[infec_person] == 1)) {
      hosp_pop[infec_person]  = 0;
      recovered[infec_person] = 1;
      *num_recovered          = *num_recovered + 1;
      *recovered_hosp         = *recovered_hosp + 1;
      num_recovered_county[county[infec_person]]++;
      num_recovered_age[(int)floor(age[infec_person] / 5)]++;
      num_recovered_hosp_county[county[infec_person]]++;
      num_recovered_hosp_age[(int)floor(age[infec_person] / 5)]++;

      if (job_status[infec_person] == 5) {
        *num_recovered_HCW  = *num_recovered_HCW + 1;
        *recovered_hosp_HCW = *recovered_hosp_HCW + 1;
      }
    } else if ((tau[infec_person] >= t - 25) &&
               (tau[infec_person] < t - 25 + dt) &&
               (hosp_pop[infec_person] == 2)) {
      hosp_pop[infec_person]  = 0;
      recovered[infec_person] = 1;
      *num_recovered          = *num_recovered + 1;
      *recovered_icu          = *recovered_icu + 1;
      num_recovered_county[county[infec_person]]++;
      num_recovered_age[(int)floor(age[infec_person] / 5)]++;
      num_recovered_icu_county[county[infec_person]]++;
      num_recovered_icu_age[(int)floor(age[infec_person] / 5)]++;

      if (job_status[infec_person] == 5) {
        *num_recovered_HCW = *num_recovered_HCW + 1;
        *recovered_icu_HCW = *recovered_icu_HCW + 1;
      }
    }
  }
}

int death(double  t,
          int     num_infectious,
          int    *infectious,
          double *tau,
          int    *dead,
          int    *icu_pop,
          int    *hosp_pop,
          int    *symptomatic,
          int     num_dead,
          double *age,
          double  dt,
          int    *num_dead_county,
          int    *num_dead_age,
          int    *county,
          int    *num_dead_HCW,
          int    *job_status) {
  double fatal_in_icu        = 0.5;
  double fatal_symptomatic[] =
  { 0.0000161, 0.0000695, 0.000309, 0.000844, 0.00161, 0.00595, 0.0193, 0.0428,
    0.078 };
  int age_group;
  int infec_person;
  int i;

  for (i = 0; i < num_infectious; i++) {
    infec_person = infectious[i];

    if ((tau[infec_person] >= t - 15) && (tau[infec_person] < t - 15 + dt)) {
      if (icu_pop[infec_person] == 1) {
        icu_pop[infec_person] = 0;

        if (COV_rand() < fatal_in_icu) {
          dead[infec_person] = 1;
          num_dead++;
          hosp_pop[infec_person] = 0;
          num_dead_county[county[infec_person]]++;
          num_dead_age[(int)floor(age[infec_person] / 5)]++;

          if (job_status[infec_person] == 5) {
            *num_dead_HCW = *num_dead_HCW + 1;
          }
        }
      } else if (symptomatic[infec_person] == 1) {
        age_group = floor(age[infec_person] / 10);

        if (COV_rand() < fatal_symptomatic[age_group]) {
          dead[infec_person]     = 1;
          icu_pop[infec_person]  = 0;
          hosp_pop[infec_person] = 0;
          num_dead++;
          num_dead_county[county[infec_person]]++;
          num_dead_age[(int)floor(age[infec_person] / 5)]++;

          if (job_status[infec_person] == 5) {
            *num_dead_HCW = *num_dead_HCW + 1;
          }
        }
      }
    }
  }
  return num_dead;
}

#if COV_GPU

void locale_infectious_loop(int            num_locale,
                            int            population,
                            int            num_households,
                            int            num_infectious,
                            int           *infectious,
                            double const  *infect_kappa,
                            double         Ic,
                            int           *intervene,
                            double         t,
                            double        *tau,
                            double        *tauI,
                            double        *interIc,
                            double         dt,
                            int           *hosp_pop,
                            int           *icu_pop,
                            double        *lat_locale,
                            double        *lon_locale,
                            int           *locale_HH,
                            int           *HH,
                            struct locale *locale_list,
                            double         omega,
                            double         rho,
                            int           *severe,
                            double         betac_scale,
                            double        *commun_nom1,
                            double        *fd_tot);

#endif /* if COV_GPU */

const char *desc[] =
{
  "This program calculates the spread of COVID-19 across Sweden.",
  "Originally published in DOI:10.1101/2020.04.11.20062133",
  "Usage as follows:",
  "./covid19 -pop $pop -sim_time $days -dt $time_steps -inter $intervention -tauI $tau_onset -R0 $R0_scaling -initial_infect_file $filename -symptomatic_percent $per_symptomatic",
  "Where all arguments are optional and described below:",
  "$pop is the population size to be simulated.  Populations smaller than 100,000 and larger than 10,098,554 are not supported.  The full population and default value is 10,098,554.",
  "$days is the total simulation time in days.  Default value is 200 days.",
  "$time_step is the time step in fraction of days. Default value is 1 day time steps.",
  "$intervention is the intervention number to apply.  Default is unmitigated spread (intervention=0).",
  "$tau_onset is the onset of the intervention in days. Default is 0, intervention starts at beginning of simulation.",
  "$R0_scaling is a scaling factor for the infection transmission rate.  R0_scaling=2.2 adjusts to a 3 day doubling time and R0_scaling=1.8 adjusts to a 5 day scaling time.  Other values are not fully vetted and doubling time would need to be determined by the user.",
  "$filename is a file pointing to the initial infections per day with one value per line for 15 days prior to simulation start in reverse chronological order.",
  "$per_symptomatic is the fraction of symptomatic individuals in float format.  Default is 0.67.",
};


int main(int argc, char *argv[]) {
  /* Parameters available to change with command line */
  int population           = 10098554;                               // Size of
                                                                     // total
                                                                     // population
  int tot_time             = 200;                                    //
                                                                     // Simulation
                                                                     // time.
  int initial_infections[] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; // Initial
                                                                     // infections
                                                                     // per
                                                                     // county.
  double percent_infect = 0.001;                                     // Default
                                                                     // 1% of
                                                                     // population
                                                                     // initially
                                                                     // infected.
  double dt             = 1.00;                                      // Time
                                                                     // step.
  int    interventions  = 0;                                         // Value of
                                                                     // interventions.
  double tauI_onset     = 0;                                         // time
                                                                     // after
                                                                     // start of
                                                                     // simulation
                                                                     // that
                                                                     // interventions
                                                                     // for
                                                                     // whole
                                                                     // community
                                                                     // take
                                                                     // place.
  int    print_lat_lon  = 0;                                         // Choose
                                                                     // whether
                                                                     // to print
                                                                     // latitude
                                                                     // and
                                                                     // longitude
                                                                     // data of
                                                                     // infected
                                                                     // individuals
                                                                     // for each
                                                                     // time
                                                                     // step.
  int    ret;
  int    i, j;
  char  *initial_infect_filename = NULL;
  char  *initial_immune_filename = NULL;

  /* Timing variables */
  struct timespec T1, T2, t1, t2, t3, t4;
  double tt, step_time, nsdiv = 1000 * 1000 * 1000;

  double symptomatic_per = 0.67; // percent of people who are symptomatic.

  /* Print optional help */
  if (argc == 1)
  {
    for (i = 0; i < sizeof(desc) / sizeof(char *); i++)
    {
            fprintf(stderr, "%s\n", desc[i]);
    }
            fprintf(stderr, "Running with default arguments\n");
  }

  /* Parse command-line arguments */
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-pop")) population = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-sim_time")) tot_time = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-infect")) percent_infect = atof(argv[++i]);
    else if (!strcmp(argv[i], "-dt")) dt = atof(argv[++i]);
    else if (!strcmp(argv[i], "-inter")) interventions = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-tauI")) tauI_onset = atof(argv[++i]);
    else if (!strcmp(argv[i], "-initial")) for (j = 0; j < 21;
                                                j++) initial_infections[j] =
          atof(argv[++i]);
    else if (!strcmp(argv[i], "-print_loc")) print_lat_lon = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-betac")) betac_scale = atof(argv[++i]);
    else if (!strcmp(argv[i], "-betah")) betah_scale = atof(argv[++i]);
    else if (!strcmp(argv[i], "-betaw")) betaw_scale = atof(argv[++i]);
    else if (!strcmp(argv[i], "-R0")) R0_scale = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pop_fract")) population_fraction = atof(
        argv[++i]);
    else if (!strcmp(argv[i], "-full_kappa")) full_kappa = 1;
    else if (!strcmp(argv[i], "-use_fixed_seed")) use_fixed_seed = 1;
    else if (!strcmp(argv[i], "-use_seed")) {
      seed           = atol(argv[++i]);
      use_fixed_seed = 1;
    }
    else if (!strcmp(argv[i], "-sort_HH")) sort_household_members = 1;
    else if (!strcmp(argv[i],
                     "-initial_infect_file")) initial_infect_filename =
        argv[++i];
    else if (!strcmp(argv[i], "-symptomatic_percent")) symptomatic_per = atof(
        argv[++i]);
  }

  /* check for supported population range */
  if ((population < POP_MIN) || (population > POP_MAX))
  {
    fprintf(stderr,
            "Population %d not in supported range of %d to %d. Exiting.\n",
            population,
            POP_MIN,
            POP_MAX);
    exit(0);
  }

  /* print out population statistics to file */
  FILE *stats = fopen("stats.log", "w");

  /*Get current time and stats */
  time_t date;
  time(&date);

  fprintf(stats,
          "Date: %s \n",
          ctime(&date));
  fprintf(stats,
          "Population: %i \nSimulationTime: %i \ndt: %f \nInterventions: %i \nTau_onset: %f \n\n\n ",
          population,
          tot_time,
          dt,
          interventions,
          tauI_onset);

  ret = clock_gettime(CLOCK_MONOTONIC, &T1);

  double HH_size        = 2.2;                         // Average household size
                                                       // from OCED
  int    num_households = (int)(population / HH_size); // Number of households
                                                       // in population
  int   *HH;                                           // Households
  HH = (int *)calloc(population, sizeof(int));
  int **per_HH_members;
  per_HH_members = (int **)calloc(num_households, sizeof(int *));

  /* Population information */
  /* Specific to Sweden.  Averaging population based on total population of
    Sweden regardless of population size. */
  double tot_pop       = 10327589.; // For use of averaging counties only.
                                    //  NOTE: land scan tot is 10,098,554,
                                    // cannot use more inhabitants than that.
  char  *county_name[] =
  { "Stockholm",       "Uppsala",              "Sodermanland",
    "Ostergotland",
    "Jonkoping",
    "Kronoberg",       "Kalmar",               "Gotland",
    "Blekinge",
    "Skane",           "Halland",
    "Vastra Gotaland", "Varmland",             "Orebro",
    "Vastmanland",
    "Dalarna",
    "Gavleborg",       "Vasternorrland",       "Jamtland",
    "Vasterbotten",
    "Norrbotten" };
  double pop_county[] =
  { 2377081,   383713,  297540, 465495, 363599, 201469, 245446,  59686, 159606,
    1377827,   333848, 1725881, 282414, 304805, 275845, 287966, 287382, 245347,
    130810,    271736, 250093 };
  int num_counties = (sizeof(pop_county) / sizeof(pop_county[0]));
  int county_int[num_counties]; // Integer values for random probability
                                // distribution.
  memset(county_int, 0, num_counties * sizeof(int));
  double *lat_city;             // latitude of city i.
  double *fd_tot;               // additive kernel density function of person
                                // with respect to all other individuals
  fd_tot   = (double *)calloc(population, sizeof(double));
  lat_city = (double *)calloc(2000, sizeof(double));
  double *long_city;            // longitude of city i.
  long_city = (double *)calloc(2000, sizeof(double));
  int *city_size;               // populatino of city i.
  city_size = (int *)calloc(2000, sizeof(int));
  int *county_size;             // population of county i.
  county_size = (int *)calloc(num_counties, sizeof(int));
  char *city_names[2000];       // City name of city i.
  int  *city;                   // city of person i by integer assignment.
  city = (int *)calloc(population, sizeof(int));
  int *city_county;             // county of person i by integer assignment.
  city_county = (int *)calloc(2000, sizeof(int));
  int num_cities = 0;           // total number of cities


  /* City information for allocating schools. */

  /* Initialize and allocate arrays */
  double *age; // Age of population
  age = (double *)calloc(population, sizeof(double));
  int *age_distrib;

  /* Distribution of ages in the generated population, people probably won't be
    older than 150 */
  age_distrib = (int *)calloc(150, sizeof(int));
  int *county;        // County of each inhabitant
  county = (int *)calloc(population, sizeof(int));
  int **county_p;     /* List of persons per county */
  county_p = (int **)calloc(num_counties, sizeof(int *));
  int *job_status;    // Type of job each person holds: 0-5 for no job,
                      // preschool, elementary school, highschool/college, job,
                      // and hospital, respectively.
  job_status = (int *)calloc(population, sizeof(int));
  int *workplace;     // Workplace of each person.
  workplace = (int *)calloc(population, sizeof(int));

  int max_num_WP = 0; // max workplaces per job_status for allocating array.

  /* Parameters for infections */
  int *infected;      // 1 if person i has been infected, 0 otherwise
  infected = (int *)calloc(population, sizeof(int));
  int *severe;        // 1 if person i has severe infection, 0 otherwise
  severe = (int *)calloc(population, sizeof(int));
  double *tau;        // day that person i became infected
  tau = (double *)calloc(population, sizeof(double));
  int *symptomatic;   // 1 if person i is symptomatic, 0 otherwise
  symptomatic = (int *)calloc(population, sizeof(int));

  /* Infection parameters.  Mostly counters for looping. */
  int  num_sus         = 0; // Number susceptible.
  int  num_sus_HCW     = 0; // Number susceptible healthcare workers.
  int  num_sus_age[18] = { 0 };
  int  num_infect      = 0; // Number infected.
  int  num_infect_HCW  = 0; // Number infected.
  int *num_sus_county;
  num_sus_county = (int *)calloc(num_counties, sizeof(int));
  int *num_infect_county;
  num_infect_county = (int *)calloc(num_counties, sizeof(int));
  int  num_infect_age[18] = { 0 };
  int  num_icu            = 0; // Number in ICU.
  int  num_icu_HCW        = 0; // Number in ICU.
  int *num_icu_county;
  num_icu_county = (int *)calloc(num_counties, sizeof(int));
  int  num_icu_age[18] = { 0 };
  int  num_hosp        = 0; // Number in hospital.
  int  num_hosp_HCW    = 0; // Number in hospital.
  int *num_hosp_county;
  num_hosp_county = (int *)calloc(num_counties, sizeof(int));
  int  num_hosp_age[18] = { 0 };
  int  num_dead         = 0; // Number of deaths.
  int  num_dead_HCW     = 0; // Number of deaths.
  int *num_dead_county;
  num_dead_county = (int *)calloc(num_counties, sizeof(int));
  int  num_dead_age[18]  = { 0 };
  int  num_recovered     = 0; // Number of people recovered.
  int  num_recovered_HCW = 0; // Number of people recovered.
  int *num_recovered_county;
  num_recovered_county = (int *)calloc(num_counties, sizeof(int));
  int  num_recovered_age[18] = { 0 };
  int *num_recovered_hosp_county;
  num_recovered_hosp_county = (int *)calloc(num_counties, sizeof(int));
  int  num_recovered_hosp_age[18] = { 0 };
  int *num_recovered_icu_county;
  num_recovered_icu_county = (int *)calloc(num_counties, sizeof(int));
  int  num_recovered_icu_age[18] = { 0 };
  int  num_infectious            = 0; // Number infectious.
  int  num_infectious_HCW        = 0; // Number infectious.
  int *num_infectious_county;
  num_infectious_county = (int *)calloc(num_counties, sizeof(int));
  int num_infectious_age[18] = { 0 };
  int recovered_hosp         = 0; // Number of people recovered.
  int recovered_icu          = 0; // Number of people recovered.
  int num_recovered_hosp_HCW = 0; // Number of people recovered.
  int num_recovered_icu_HCW  = 0; // Number of people recovered.


  /* This list hold indices for each category. */
  int *icu_list;   // Indices in ICU.
  icu_list = (int *)calloc(population, sizeof(int));
  int *hosp_list;  // Indices in hospital.
  hosp_list = (int *)calloc(population, sizeof(int));
  int *infectious; // Indices of infected.
  infectious = (int *)calloc(population, sizeof(int));
  int *sus_list;   // Indices of susceptible.
  sus_list = (int *)calloc(population, sizeof(int));

  /* Parameters pertaining to the population.  These could be made into a struct
    later to look nicer. */
  int *hosp_pop;            // 1 if person i in hospital, 0 otherwise
  hosp_pop = (int *)calloc(population, sizeof(int));
  int *icu_pop;             // 1 if person i in ICU, 0 otherwise
  icu_pop = (int *)calloc(population, sizeof(int));
  int *recovered;           // 1 if person i is recovered, 0 otherwise.
  recovered = (int *)calloc(population, sizeof(int));
  int *immune;              // 0 if susceptible, >0 if not.
  immune = (int *)calloc(population, sizeof(int));
  int *dead;                // 1 if person i is dead, 0 otherwise
  dead = (int *)calloc(population, sizeof(int));
  int *workplace_tmp;       // Hospital location when people go into hospital.
  workplace_tmp = (int *)calloc(population, sizeof(int));

  double community_nom = 0; // For adding community infection.
  int    file_count;

    #define num_I 14
  double Ic = 1.0;              // Intervention constant for community
                                // transmission.
  const double(*Iw);            // Intervention constant for workplace
                                // transmission.
  double Ih = 1;                // Intervention constant for household
                                // transmission.
  double complyI[num_I];        // percent of people who comply with
                                // intervention
  double tauI[num_I];           // time after infection that intervention takes
                                // place
  double interIc[num_I];        // Intervention constants for Ic
  double interIh[num_I];        // Intervention constants for Ih
  int    personinter[num_I];    // Tells us whether the intervention needs to be
                                // calculated on a person to person basis or for
                                // the whole community.
  int   *intervene;             // 1 if person is currently undergoing
                                // interventions, 0 otherwise.
  intervene = (int *)calloc(population, sizeof(int));
  int randomized = 0;           // for interventions that keep track of
                                // randomization
  const double *interIw[num_I]; /* InterIw array.*/

  /* Setup intervention types */
  setup_intervention_types(interIc,
                           &Iw,
                           interIh,
                           complyI,
                           tauI,
                           personinter,
                           interIw);

  /**** Set random number generator seed. ****/
  COV_init_rand();

  /* Initialize age distribution */
  age_dist(age, population, stats, age_distrib);
  city_lat_long(&num_cities,
                lat_city,
                long_city,
                city_names,
                city_county,
                county_name,
                num_counties);

  /* Parse land_scan file to get population density.  */
  double *lat_locale           = (double *)calloc(max_locale + 10,
                                                  sizeof(double));
  double *lon_locale           = (double *)calloc(max_locale + 10,
                                                  sizeof(double));
  double *pop_density_init_num = (double *)calloc(max_locale + 10,
                                                  sizeof(double));
  double tmp_lat, tmp_lon, pop_den, land_pop_total_density = 0;
  const char *popfile = "land_pop_sorted.txt";
  FILE *lat_long      = fopen(popfile, "r"); // Sorted land population in
                                             // descending order.  Important
                                             // when we don't have complete
                                             // population.

  if (lat_long == NULL)
  {
    fprintf(stderr, "Error opening file %s.  Exiting\n", popfile);
    exit(0);
  }

  while ((ret =
            fscanf(lat_long, "%lf%*c%lf%*c%lf", &tmp_lon, &tmp_lat,
                   &pop_den)) == 3) {
    if (num_locale + 1 > max_locale) {
      // max_locale += 10;
      lat_locale           = (double *)realloc(lat_locale,
                                               (num_locale +
                                                10) * sizeof(double));
      lon_locale           = (double *)realloc(lon_locale,
                                               (num_locale +
                                                10) * sizeof(double));
      pop_density_init_num =
        (double *)realloc(pop_density_init_num,
                          (num_locale + 10) * sizeof(double));
    }
    lat_locale[num_locale]           = tmp_lat;
    lon_locale[num_locale]           = tmp_lon;
    pop_density_init_num[num_locale] = pop_den;
    add_locale(tmp_lat, tmp_lon, pop_den);
    land_pop_total_density += pop_den;
    num_locale++;
  }
  fclose(lat_long);

  int   original_num_locale = num_locale;
  int **locale_to_HH;
  int  *locale_HH;
  locale_to_HH = (int **)calloc(num_locale, sizeof(int *));
  locale_HH    = (int *)calloc(num_households, sizeof(int));

  /* Initialize households */
  household_lat_long(num_households,
                     HH,
                     lat_city,
                     long_city,
                     num_cities,
                     city,
                     county,
                     city_county,
                     city_size,
                     county_size,
                     population,
                     age,
                     per_HH_members,
                     county_name,
                     pop_county,
                     tot_pop,
                     stats,
                     county_p,
                     &num_locale,
                     lat_locale,
                     lon_locale,
                     pop_density_init_num,
                     locale_to_HH,
                     locale_HH,
                     land_pop_total_density);

  free(city_county);


  /* Open files */
  char  *file_beg = "county_";
  char  *file_end = ".log";
  char   file_name[1000];
  FILE **county_files = malloc(num_counties * sizeof(FILE *));

  for (i = 0; i < num_counties; i++) {
    strcpy(file_name, "");
    strcat(file_name, file_beg);
    strcat(file_name, county_name[i]);
    strcat(file_name, file_end);
    county_files[i] = fopen(file_name, "w");
  }
  file_beg = "age_";
  file_end = ".log";
  char   str;
  FILE **age_files = malloc(sizeof(FILE *) * 18);

  for (i = 0; i < 90; i += 5) {
    str = i;
    strcpy(file_name, "");
    strcat(file_name, file_beg);
    sprintf(file_name, "%s%i", file_name, i);
    strcat(file_name, file_end);
    age_files[(i / 5)] = fopen(file_name, "w");
  }
  FILE *output_file = fopen("covid19_spread.dat", "w");
  FILE *output_HCW  = fopen("healthcare_workers.log", "w");
  FILE *lat_lon_out = fopen("lat_lon.dat", "w");

  place_initial_infections(initial_infect_filename,
                           initial_immune_filename,
		           stats,
                           initial_infections,
                           tau,
                           infected,
                           immune,
                           severe,
                           symptomatic,
                           county,
                           &num_infect,
                           num_counties,
                           symptomatic_per,
                           population,
                           dt,
                           lat_locale,
                           lon_locale,
                           num_infect_county,
                           num_infect_age,
                           age,
                           county_p,
                           county_size,
                           locale_HH,
                           HH);

  for (size_t i = 0; i < num_counties; ++i) {
    arrfree(county_p[i]);
  }
  free(county_p);
  county_p = NULL;

  int **job_status_county; // Jobs per county, or city for schools
  job_status_county = (int **)calloc(6, sizeof(int *));

  for (i = 0; i < 6; i++) job_status_county[i] = (int *)calloc(num_cities,
                                                               sizeof(int));

  /* Initialize job/school status */
  job_dist(job_status,
           job_status_county,
           age,
           county,
           city,
           population,
           num_cities,
           num_counties,
           county_size,
           stats);

  int num_HCW = 0;

  // Get number of healthcare workers //
  for (i = 0; i < population; i++) {
    if (job_status[i] == 5) {
      num_HCW++;
    }
  }

  int hosp_num[num_counties]; // Number of hospitals per county
  memset(hosp_num, 0, num_counties);

  /* Initialize workplace/school */
  workplace_dist(workplace,
                 job_status,
                 job_status_county,
                 city,
                 num_cities,
                 county,
                 num_counties,
                 population,
                 &max_num_WP,
                 hosp_num,
                 stats);

  free(city);

  /* Get size of each workplace as an array.  Cannot allocate array until
    max_num_WP is known. */
  int *workplace_size[6];

  for (size_t i = 0; i < 6; ++i) {
    workplace_size[i] = (int *)calloc(max_num_WP, sizeof(int));
  }

  for (i = 0; i < population; i++) {
    workplace_size[job_status[i]][workplace[i]]++;
  }


  for (j = 1; j < 6; j++) {
    int avg_work_size = 0, num_WP = 0;

    for (i = 0; i < max_num_WP; i++) {
      if (workplace_size[j][i] > 0) {
        avg_work_size += workplace_size[j][i];
        num_WP++;
      }
    }
    fprintf(stats,
            "Job_status %i avg_WP_size %f \n",
            j,
            (double)avg_work_size / num_WP);
  }
  fclose(stats);


  printf("Starting density kernel calculations\n");
  fflush(stdout);

  ret = clock_gettime(CLOCK_MONOTONIC, &t1);
  int *npl = (int *)malloc(num_locale * sizeof(int));
#ifdef _OPENMP
# pragma omp parallel for private(i) default(shared)
#endif /* ifdef _OPENMP */

  for (i = 0; i < num_locale; i++) {
    int npi; /* number of persons in locale i */
    npi = 0;

    for (int hh = 0; hh < arrlen(locale_to_HH[i]); hh++) {
      npi += arrlen(per_HH_members[locale_to_HH[i][hh]]);
    }
    npl[i] = npi;
  }

  /* Precalculate total density kernel function for each individual */
  for (i = 0; i < num_locale; i++) {
    double itmp_fd;
    itmp_fd = 0;
#ifdef _OPENMP
# pragma omp parallel for private(j) default(shared) reduction(+:itmp_fd)
#endif /* ifdef _OPENMP */

    for (j = i + 1; j < num_locale; j++) {
      double d;
      double tmp_fd;
#if !defined(USE_LOCALE_DISTANCE)
      d = distance(lat_locale[i],
                   lon_locale[i],
                   lat_locale[j],
                   lon_locale[j],
                   'K');
#else /* if !defined(USE_LOCALE_DISTANCE) */
      d = locale_distance(locale_list[i], locale_list[j]);
#endif /* if !defined(USE_LOCALE_DISTANCE) */
      tmp_fd     = 1 / (1 + pow((d / 4), 3)); // kernel density function as
                                              // parameterized for GB.
      itmp_fd   += tmp_fd * npl[j];
      fd_tot[j] += tmp_fd * npl[i];
    }
    fd_tot[i] += itmp_fd + npl[i] - 1;
  }
  free(npl);
  ret = clock_gettime(CLOCK_MONOTONIC, &t2);
  tt  =
    ((double)t2.tv_sec + (double)t2.tv_nsec /
     nsdiv) - ((double)t1.tv_sec + (double)t1.tv_nsec / nsdiv);
  printf("Done with density kernel calculations in %5.2f\n", tt);
  fflush(stdout);

  for (size_t i = 0; i < original_num_locale; ++i) {
    arrfree(locale_to_HH[i]);
  }
  free(locale_to_HH);


  /* Initialize constants */
  double alpha = 0.8; // From Ferguson Nature 2006
  


  /* CHANGING OMEGA */
  double omega = 1;   // From Ferguson Nature 2006


  // #Leaving out rho from Ferguson 2006.  This is a measure of how infectious
  // person is.  For now, we will assume all people are the same.
  
  /* CHANGING RHO */

  const gsl_rng_type* T;
  gsl_rng* rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  double rho_mean = 1.0;
  double std_dev = 4.0;
   

  double variance = std_dev * std_dev;
  double shape = (rho_mean * rho_mean) / variance;
  double scale_parameter = variance / rho_mean;

  double rho = gsl_ran_gamma(rng, shape, scale_parameter);
  gsl_rng_free(rng);


  /* precalculate kappa */
  int count_kappa_vals = ceil(30 / dt);
  double  kappa_t      = 0;
  double  tmp_t        = 0;
  double *kappa_vals;
  double  tau1 = 0;
  kappa_vals = (double *)calloc(count_kappa_vals, sizeof(double));

  for (i = 1; i < count_kappa_vals; i++) {
    kappa_t = i * dt - 4.6;

    if (kappa_t <= 0) {
      kappa_vals[i] = 0;
    } else {
      tmp_t         = (log(kappa_t) + 0.72) / 1.8;
      kappa_vals[i] =
        exp(-0.5 * pow(tmp_t, 2)) / ((kappa_t) * 1.8 * sqrt(2 * pi));
    }
  }

  int contact_commun         = 0;
  int contact_work           = 0;
  int contact_school         = 0;
  int contact_house          = 0;
  int num_contact_commun     = 0;
  int num_contact_work       = 0;
  int num_contact_school     = 0;
  int num_contact_house      = 0;
  int num_contact_commun_HCW = 0;
  int num_contact_work_HCW   = 0;
  int num_contact_school_HCW = 0;
  int num_contact_house_HCW  = 0;
  int num_contact_commun_county[num_counties];
  int num_contact_commun_age[18];
  int num_contact_work_county[num_counties];
  int num_contact_work_age[18];
  int num_contact_school_county[num_counties];
  int num_contact_school_age[18];
  int num_contact_house_county[num_counties];
  int num_contact_house_age[18];

  for (i = 0; i < num_counties; i++) {
    num_dead_county[i]           = 0;
    num_recovered_county[i]      = 0;
    num_recovered_hosp_county[i] = 0;
    num_recovered_icu_county[i]  = 0;
  }

  for (i = 0; i < 18; i++) {
    num_dead_age[i]           = 0;
    num_recovered_age[i]      = 0;
    num_recovered_hosp_age[i] = 0;
    num_recovered_icu_age[i]  = 0;
  }

  double *commun_nom1 = (double *)malloc(num_locale * sizeof(double));
  double *work_infect[6];

  for (size_t i = 0; i < 6; ++i) {
    work_infect[i] = (double *)malloc(max_num_WP * sizeof(double));
  }
  double *house_infect = (double *)malloc(num_households * sizeof(double));


  /* NEED TO INITIALIZE PREVIOUS INFECTIONS TO RECOVERED HOSP AND ICU */
  Ic = 1;
  double t         = 0;
  double time_step = 0;

  /* Segment population into infectious, susceptible, hospitalized, and icu */
  segment_population(&num_sus,
                     &num_infectious,
                     &num_hosp,
                     &num_icu,
                     &num_sus_HCW,
                     &num_infectious_HCW,
                     &num_hosp_HCW,
                     &num_icu_HCW,
                     infected,
                     immune,
  		     infectious,
                     sus_list,
                     hosp_list,
                     hosp_pop,
                     icu_pop,
                     icu_list,
                     tau,
                     population,
                     t,
                     num_sus_county,
                     num_infectious_county,
                     num_infect_county,
                     num_hosp_county,
                     num_icu_county,
                     num_sus_age,
                     num_infectious_age,
                     num_infect_age,
                     num_hosp_age,
                     num_icu_age,
                     age,
                     county,
                     print_lat_lon,
                     lat_locale,
                     lon_locale,
                     lat_lon_out,
                     locale_HH,
                     HH,
                     job_status,
                     recovered,
                     dead);

  // After 5 days of symptoms, people are randomly put into hospital based on
  // age distributions. //
  hosp_entry(t,
             num_infectious,
             infectious,
             age,
             icu_pop,
             hosp_pop,
             symptomatic,
             tau,
             workplace_tmp,
             hosp_num,
             county,
             dt);

  // Release from hospital //
  hosp_release(t,
               num_hosp,
               hosp_list,
               tau,
               recovered,
               hosp_pop,
               &num_recovered,
               &recovered_hosp,
               &recovered_icu,
               dt,
               num_recovered_county,
               num_recovered_age,
               age,
               county,
               num_recovered_hosp_county,
               num_recovered_hosp_age,
               num_recovered_icu_county,
               num_recovered_icu_age,
               &num_recovered,
               &recovered_hosp,
               &recovered_icu,
               job_status);


  // 15 days after onset of symptoms, people randomly die based on age based
  // percentages.  ICU patients are removed from ICU into regular hospital. //

  num_dead = death(t,
                   num_infectious,
                   infectious,
                   tau,
                   dead,
                   icu_pop,
                   hosp_pop,
                   symptomatic,
                   num_dead,
                   age,
                   dt,
                   num_dead_county,
                   num_dead_age,
                   county,
                   &num_dead_HCW,
                   job_status);
  num_dead = death(t,
                   num_icu,
                   icu_list,
                   tau,
                   dead,
                   icu_pop,
                   hosp_pop,
                   symptomatic,
                   num_dead,
                   age,
                   dt,
                   num_dead_county,
                   num_dead_age,
                   county,
                   &num_dead_HCW,
                   job_status);
  num_dead = death(t,
                   num_hosp,
                   hosp_list,
                   tau,
                   dead,
                   icu_pop,
                   hosp_pop,
                   symptomatic,
                   num_dead,
                   age,
                   dt,
                   num_dead_county,
                   num_dead_age,
                   county,
                   &num_dead_HCW,
                   job_status);


  // Recovered after 11 days (6 days of symptoms) if not in hospital/ICU. //
  for (i = 0; i < num_infectious; i++) {
    int infec_person;
    infec_person = infectious[i];

    if ((tau[infec_person] < (t - 11) + dt) && (hosp_pop[infec_person] == 0) &&
        (icu_pop[infec_person] == 0)) {
      recovered[infec_person] = 1;
      num_recovered++;

      if (job_status[infec_person] == 5) {
        num_recovered_HCW++;
      }
      num_recovered_county[county[infec_person]]++;
      num_recovered_age[(int)floor(age[infec_person] / 5)]++;
    }
  }

  /* Introduce overall community interventions. */
  if (interventions == interNone) {
    Ic = interIc[interventions];
    Ih = interIh[interventions];
    Iw = interIw[interventions];
  } else if (interventions > interNone)  {
    for (i = 0; i < population; i++) {
      /* high school and university closures */
      if ((age[i] >= 15) && (age[i] < 22)) {
        intervene[i] = 1;
      } else if ((age[i] >= 70) && (COV_rand() < complyI[8])) {
        intervene[i] = 8;
      }
    }
  }

  ret = clock_gettime(CLOCK_MONOTONIC, &T2);
  tt  =
    ((double)T2.tv_sec + (double)T2.tv_nsec /
     nsdiv) - ((double)T1.tv_sec + (double)T1.tv_nsec / nsdiv);
  printf("Total initialization time %5.2f\n", tt);

  /* Start simulation */
  printf("Starting simulation\n");
  fflush(stdout);

  /* Precalculated kappa for infected, recalculated each timestep */
  double *infect_kappa    = NULL;
  int     sz_infect_kappa = 0;

  ret = clock_gettime(CLOCK_MONOTONIC, &T1);

  for (time_step = 0; time_step < (tot_time / dt); time_step++) {
    ret = clock_gettime(CLOCK_MONOTONIC, &t1);
    printf("Timestep %5.2f ", t);
    fflush(stdout);

    /* count origin of contacts */
    num_contact_commun     = 0;
    num_contact_work       = 0;
    num_contact_school     = 0;
    num_contact_house      = 0;
    num_contact_commun_HCW = 0;
    num_contact_work_HCW   = 0;
    num_contact_school_HCW = 0;
    num_contact_house_HCW  = 0;

    for (i = 0; i < num_counties; i++) {
      num_contact_commun_county[i] = 0;
      num_contact_work_county[i]   = 0;
      num_contact_school_county[i] = 0;
      num_contact_house_county[i]  = 0;
      num_sus_county[i]            = 0;
      num_hosp_county[i]           = 0;
      num_icu_county[i]            = 0;
      num_infectious_county[i]     = 0;
    }

    for (i = 0; i < 18; i++) {
      num_contact_house_age[i]  = 0;
      num_contact_school_age[i] = 0;
      num_contact_work_age[i]   = 0;
      num_contact_commun_age[i] = 0;
      num_sus_age[i]            = 0;
      num_hosp_age[i]           = 0;
      num_icu_age[i]            = 0;
      num_infectious_age[i]     = 0;
    }

    if ((interventions == interSchoolClosure) && (t >= tauI_onset) && (t <= tauI_onset + dt)) {
      for (i = 0; i < population; i++) {
        if ((age[i] > 1) && (age[i] < 15)) {
          intervene[i] = 2;
        }
      }
    } else if ((interventions == interNonEssentialClosure) && (t >= tauI_onset) &&
               (t <= tauI_onset + dt)) {
      Ic = interIc[7];
      Ih = interIh[7];
      Iw = interIw[7];
    } else if ((interventions == interClosureAndDistancing) && (t >= tauI_onset) &&
               (t <= tauI_onset + dt)) {
      /* if not complying, have same interactions as type 7 */
      Ic = interIc[7];
      Ih = interIh[7];
      Iw = interIw[7];

      for (i = 0; i < population; i++) {
        if (COV_rand() < 0.9) {
          intervene[i] = 6;
        }
      }
    } else if ((interventions == interWFH50PlusDistancing) && (t >= tauI_onset) &&
               (t <= tauI_onset + dt)) {
      Iw = interIw[9];

      /* Intervention is social distancing, 50% work from home, and schools
        open.  */
      for (i = 0; i < population; i++) {
        if (COV_rand() < complyI[9]) {
          intervene[i] = 9;
        }
      }
    } else if ((interventions == interGoogleMobilityPlusDistancing) && (t >= tauI_onset) &&
               (t <= tauI_onset + dt)) {
      /* Intervention is social distancing with numbers from Google mobility,
        24% work from home, 41% reduction in community interaction, and schools
        open.  */
      Ic = interIc[10];
      Ih = interIh[10];
      Iw = interIw[10];

      for (i = 0; i < population; i++) {
        if ((age[i] >= 15) && (age[i] < 22)) {
          intervene[i] = 0;
        }
      }
    } else if ((interventions == interVolWFH) && (t >= tauI_onset) &&
               (t <= tauI_onset + dt) && (randomized != interVolWFH)) {
      /* Intervention 6 is scalable A/C switching:   */
      /* Don't modify baseline Ic, Ih, IW.  Those are then "type A" */
      /* Only randomize once. */
      randomized = interVolWFH;
      int ctr = 0;

      for (i = 0; i < population; i++) {
        if (COV_rand() < population_fraction) {
          intervene[i] = 12;
          ctr++;
        }
      }
      printf("Modified %d people\n", ctr);
    } else if ((interventions == interVolSelfIsolation) && (t >= tauI_onset) &&
               (t <= tauI_onset + dt) && (randomized != interVolSelfIsolation)) {
      /* Intervention 7 is scalable A/D switching:   */
      /* Don't modify baseline Ic, Ih, IW.  Those are then "type A" */
      randomized = interVolSelfIsolation;

      for (i = 0; i < population; i++) {
        if (COV_rand() < population_fraction) {
          intervene[i] = 13;
        }
      }
    } else if ((interventions == interVolWFHDistancing) && (t >= tauI_onset) &&
               (t <= tauI_onset + dt) && (randomized != interVolWFHDistancing)) {
      /* Intervention 8 is scalable B/C switching:   */
      /* Don't modify baseline Ic, Ih, IW.  Those are then "type A" */
      randomized = interVolWFHDistancing;
      Ic         = interIc[11];
      Ih         = interIh[11];
      Iw         = interIw[11];

      for (i = 0; i < population; i++) {
        if (COV_rand() < population_fraction) {
          intervene[i] = 12;
        }
      }
    } else if ((interventions == interVolSelfIsolDistancing) && (t >= tauI_onset) &&
               (t <= tauI_onset + dt) && (randomized != interVolSelfIsolDistancing)) {
      /* Intervention 9 is scalable B/D switching:   */
      /* Don't modify baseline Ic, Ih, IW.  Those are then "type A" */
      randomized = interVolSelfIsolDistancing;
      Ic         = interIc[11];
      Ih         = interIh[11];
      Iw         = interIw[11];

      for (i = 0; i < population; i++) {
        if (COV_rand() < population_fraction) {
          intervene[i] = 13;
        }
      }
    }

    /* Segment population into infectious, susceptible, hospitalized, and icu */
    segment_population(&num_sus,
                       &num_infectious,
                       &num_hosp,
                       &num_icu,
                       &num_sus_HCW,
                       &num_infectious_HCW,
                       &num_hosp_HCW,
                       &num_icu_HCW,
                       infected,
                       immune,
    		       infectious,
                       sus_list,
                       hosp_list,
                       hosp_pop,
                       icu_pop,
                       icu_list,
                       tau,
                       population,
                       t,
                       num_sus_county,
                       num_infectious_county,
                       num_infect_county,
                       num_hosp_county,
                       num_icu_county,
                       num_sus_age,
                       num_infectious_age,
                       num_infect_age,
                       num_hosp_age,
                       num_icu_age,
                       age,
                       county,
                       print_lat_lon,
                       lat_locale,
                       lon_locale,
                       lat_lon_out,
                       locale_HH,
                       HH,
                       job_status,
                       recovered,
                       dead);


      memset(commun_nom1,    0, num_locale * sizeof(double));

    for (i = 0; i < 6; ++i) {
      memset(work_infect[i], 0, max_num_WP * sizeof(double));
    }
      memset(house_infect,   0, num_households * sizeof(double));

    if (num_infectious > sz_infect_kappa) {
      sz_infect_kappa = num_infectious;
      infect_kappa    = (double *)realloc(infect_kappa,
                                          sz_infect_kappa * sizeof(double));
    }

    for (i = 0; i < num_infectious; i++) {
      int infec_person;
      infec_person    = infectious[i];
      infect_kappa[i] =
        calc_kappa(t,
                   tau[infec_person],
                   symptomatic[infec_person],
                   dt,
                   kappa_vals,
                   hosp_pop[infec_person],
                   icu_pop[infec_person]);
    }

#if COV_GPU

    locale_infectious_loop(num_locale,
                           population,
                           num_households,
                           num_infectious,
                           infectious,
                           infect_kappa,
                           Ic,
                           intervene,
                           t,
                           tau,
                           tauI,
                           interIc,
                           dt,
                           hosp_pop,
                           icu_pop,
                           lat_locale,
                           lon_locale,
                           locale_HH,
                           HH,
                           locale_list,
                           omega,
                           rho,
                           severe,
                           betac_scale,
                           commun_nom1,
                           fd_tot);

#else // COV_GPU

# ifdef _OPENMP
#  pragma omp parallel for private(j, i) default(shared)
# endif /* ifdef _OPENMP */

    for (j = 0; j < num_locale; j++) {
      double tmp_comm_inf = 0;

      for (i = 0; i < num_infectious; i++) {
        int infec_person; // Counter for infected person.
        double kappa;     // #Infectiousness
        double tIc;
        infec_person = infectious[i];
        tIc          = Ic;

        if ((intervene[infec_person] > 0) &&
            (t > tau[infec_person] + tauI[intervene[infec_person]])) {
          tIc = interIc[intervene[infec_person]];
        }
        kappa = infect_kappa[i];

        if (hosp_pop[infec_person] == 0) {
          double d; // distance between people.
          // Community transmission //
# if !defined(USE_LOCALE_DISTANCE)
          d =
            distance(lat_locale[j], lon_locale[j],
                     lat_locale[locale_HH[HH[infec_person]]],
                     lon_locale[locale_HH[HH[infec_person]]], 'K');
# else /* if !defined(USE_LOCALE_DISTANCE) */
          d =
            locale_distance(locale_list[j],
                            locale_list[locale_HH[HH[infec_person]]]);
# endif /* if !defined(USE_LOCALE_DISTANCE) */
          tmp_comm_inf += tIc *
                          calc_community_infect(kappa, omega, rho,
                                                severe[infec_person], d);
        }
      }
      commun_nom1[j] = tmp_comm_inf / fd_tot[j];
    }

#endif // COV_GPU

    for (i = 0; i < num_infectious; i++) {
      int infec_person;       // Counter for infected person.
      double kappa; // #Infectiousness
      double tIw, tIh;
      int    age_group    = 0;
      int    tmp_job_stat = 0;
      double tmp_work_inf = 0, tmp_house_inf = 0;
      infec_person = infectious[i];
      tIh          = Ih;
      tIw          = Iw[job_status[infec_person]];

      if ((intervene[infec_person] > 0) &&
          (t > tau[infec_person] + tauI[intervene[infec_person]])) {
        tIw = interIw[intervene[infec_person]][job_status[infec_person]];
        tIh = interIh[intervene[infec_person]];
      }
      kappa = calc_kappa(t,
                         tau[infec_person],
                         symptomatic[infec_person],
                         dt,
                         kappa_vals,
                         hosp_pop[infec_person],
                         icu_pop[infec_person]);

      if (hosp_pop[infec_person] == 0) {
        tmp_job_stat = job_status[infec_person];

        // Workplace/School transmission: People must be in same workplace and
        // job type. //
        tmp_work_inf = calc_workplace_infect(tmp_job_stat,
                                             kappa,
                                             omega, 
                                             rho,
                                             workplace_size[tmp_job_stat][
                                               workplace[infec_person]],
                                             severe[infec_person],
                                             tIw);
      } else {
        /* In hospital, only have interaction with hospital workers and half
          interaction with family (household). */
        tmp_job_stat = 5;
        tIw          = Iw[tmp_job_stat];
        tIh          = 0.25;
        tmp_work_inf = calc_workplace_infect(tmp_job_stat,
                                             kappa,
                                             omega,
                                             rho,
                                             workplace_size[tmp_job_stat][
                                               workplace_tmp[infec_person]],
                                             severe[infec_person],
                                             tIw);
      }
      work_infect[tmp_job_stat][workplace[infec_person]] += tmp_work_inf;

      // Household transmission //
      tmp_house_inf                   = tIh * calc_household_infect(kappa,
                                                                    omega,
                                                                    rho,
                                                                    arrlen(
                                                                      per_HH_members
                                                                      [HH[
      infec_person]]),
                                                                    alpha,
                                                                    severe[
                                                                      infec_person]);
      house_infect[HH[infec_person]] += tmp_house_inf;
    }


    // #### Only Susceptible people can get the virus and infected people spread
    // it.
    //
    double inf_compil_home = 0;
    double inf_compil_work = 0;
    double inf_compil_school = 0;
    double inf_compil_commun = 0;

    for (i = 0; i < num_sus; i++) {
      int sus_person;     // Counter for susceptible person.
      int age_group;
      double infect;      // Infectiousness
      double infect_prob; // Infectious probability
      double tIw, tIh, tIc;
      sus_person = sus_list[i];
      int contact_work   = 0;
      int contact_commun = 0;
      int contact_house  = 0;
      int contact_school = 0;

    
      tIh = Ih;
      tIc = Ic;
      tIw = Iw[job_status[sus_person]];

      if ((intervene[sus_person] > 0) &&
          (t > tau[sus_person] + tauI[intervene[sus_person]])) {
        tIw = interIw[intervene[sus_person]][job_status[sus_person]];
        tIh = interIh[intervene[sus_person]];
        tIc = interIc[intervene[sus_person]];
      }
      age_group = floor(age[sus_person] / 5);
      double zeta[] =
      { 0.1,    0.25,  0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.75,
        0.50,   0.25, 0.25, 0.25 }; //   # Travel related parameter for
                                    // community transmission. Ferguson Nature
                                    // 2006
      infect  = tIh * house_infect[HH[sus_person]];
      inf_compil_home += tIh * house_infect[HH[sus_person]];
      // printf("Household infect %8.3f\n",infect); 
      infect += tIw *
                work_infect[job_status[sus_person]][workplace[sus_person]];
      inf_compil_work += tIw *
                work_infect[job_status[sus_person]][workplace[sus_person]];
      //printf("With work infect %8.3f\n", infect);
      infect += tIc * zeta[age_group] * commun_nom1[locale_HH[HH[sus_person]]];
      inf_compil_commun += tIc * zeta[age_group] * commun_nom1[locale_HH[HH[sus_person]]];
      //printf("And also community infect %8.3f\n", infect);

      // ### Probability of being infected ####
      infect_prob = (1 - exp(-infect * dt));

      // ### Monte carlo type prediction ###
      if (COV_rand() < infect_prob) {
        infected[sus_person] = 1;
        tau[sus_person]      = t;
        severe[sus_person]   = round(COV_rand());

        if (COV_rand() < symptomatic_per) {
          symptomatic[sus_person] = 1;
        }

        if (tIh * house_infect[HH[sus_person]] > 0) {
          num_contact_house++;

          if (job_status[sus_person] == 5) {
            num_contact_house_HCW++;
          }
          num_contact_house_county[county[sus_person]]++;
          num_contact_house_age[(int)floor(age[sus_person] / 5)]++;
        } else if (tIw *
                   work_infect[job_status[sus_person]][workplace[sus_person]] >
                   0) {
          if (job_status[sus_person] < 4) {
            num_contact_school++;
            num_contact_school_county[county[sus_person]]++;
            num_contact_school_age[(int)floor(age[sus_person] / 5)]++;
          } else {
            num_contact_work++;

            if (job_status[sus_person] == 5) {
              num_contact_work_HCW++;
            }
            num_contact_work_county[county[sus_person]]++;
            num_contact_work_age[(int)floor(age[sus_person] / 5)]++;
          }
        } else if (tIc * zeta[age_group] *
                   commun_nom1[locale_HH[HH[sus_person]]] > 0) {
          num_contact_commun++;

          if (job_status[sus_person] == 5) {
            num_contact_commun_HCW++;
          }
          num_contact_commun_county[county[sus_person]]++;
          num_contact_commun_age[(int)floor(age[sus_person] / 5)]++;
        }
        num_infect++;
        num_infect_county[county[sus_person]]++;
        num_infect_age[(int)floor(age[sus_person] / 5)]++;

        if (job_status[sus_person] == 5) {
          num_infect_HCW++;
        }

        /* Determine if following interventions only for interventions that
          effect individuals.*/
        if ((interventions > interNone) && (COV_rand() < complyI[3])) {
          if ((interventions == interCaseHouseIsolation) && (intervene[sus_person] == 4)) {
            intervene[sus_person] = 5;
          } else {
            intervene[sus_person] = 3;
          }

          /* Intervention 2 is household quarantine with current
            recommendations. Applicable for whole household.  */
          if ((interventions == interCaseHouseIsolation) && (t > tauI_onset)) {
            int i1, hh = HH[sus_person];

            for (i1 = 0; i1 < arrlen(per_HH_members[hh]); i1++) {
              if (COV_rand() < complyI[4]) {
                intervene[per_HH_members[hh][i1]] = 4;
              }
            }
          }
        }
      }
    }

    printf("Household infect %8.3f\n",inf_compil_home); 
    printf("Work infect %8.3f\n",inf_compil_work); 
    printf("Community infect %8.3f\n",inf_compil_commun); 
    

    // free(inf_compil_home);
    //free(inf_compil_work);
    //free(inf_compil_commun);
    // After 5 days of symptoms, people are randomly put into hospital based on
    // age distributions. //
    hosp_entry(t,
               num_infectious,
               infectious,
               age,
               icu_pop,
               hosp_pop,
               symptomatic,
               tau,
               workplace_tmp,
               hosp_num,
               county,
               dt);

    // Release from hospital //
    hosp_release(t,
                 num_hosp,
                 hosp_list,
                 tau,
                 recovered,
                 hosp_pop,
                 &num_recovered,
                 &recovered_hosp,
                 &recovered_icu,
                 dt,
                 num_recovered_county,
                 num_recovered_age,
                 age,
                 county,
                 num_recovered_hosp_county,
                 num_recovered_hosp_age,
                 num_recovered_icu_county,
                 num_recovered_icu_age,
                 &num_recovered,
                 &recovered_hosp,
                 &recovered_icu,
                 job_status);


    // 15 days after onset of symptoms, people randomly die based on age based
    // percentages.  ICU patients are removed from ICU into regular hospital. //

    num_dead = death(t,
                     num_infectious,
                     infectious,
                     tau,
                     dead,
                     icu_pop,
                     hosp_pop,
                     symptomatic,
                     num_dead,
                     age,
                     dt,
                     num_dead_county,
                     num_dead_age,
                     county,
                     &num_dead_HCW,
                     job_status);
    num_dead = death(t,
                     num_icu,
                     icu_list,
                     tau,
                     dead,
                     icu_pop,
                     hosp_pop,
                     symptomatic,
                     num_dead,
                     age,
                     dt,
                     num_dead_county,
                     num_dead_age,
                     county,
                     &num_dead_HCW,
                     job_status);
    num_dead = death(t,
                     num_hosp,
                     hosp_list,
                     tau,
                     dead,
                     icu_pop,
                     hosp_pop,
                     symptomatic,
                     num_dead,
                     age,
                     dt,
                     num_dead_county,
                     num_dead_age,
                     county,
                     &num_dead_HCW,
                     job_status);


    // Recovered after 11 days (6 days of symptoms) if not in hospital/ICU. //
    for (i = 0; i < num_infectious; i++) {
      int infec_person;
      infec_person = infectious[i];

      if ((tau[infec_person] < (t - 11) + dt) &&
          (tau[infec_person] >= t - 11) && (hosp_pop[infec_person] == 0) &&
          (icu_pop[infec_person] == 0)) {
        recovered[infec_person] = 1;
        num_recovered++;

        if (job_status[infec_person] == 5) {
          num_recovered_HCW++;
        }
        num_recovered_county[county[infec_person]]++;
        num_recovered_age[(int)floor(age[infec_person] / 5)]++;
      }
    }

    ret       = clock_gettime(CLOCK_MONOTONIC, &t2);
    step_time =
      ((double)t2.tv_sec + (double)t2.tv_nsec /
       nsdiv) - ((double)t1.tv_sec + (double)t1.tv_nsec / nsdiv);
    printf("time %5.2f\n", step_time);
    fflush(stdout);

    fprintf(output_file,
            "Walltime/timestep %6.2f Time %6.2f num_infected %i num_infectious %i num_in_hosp %i num_in_icu %i num_dead %i recovered_tot %i recovered_from_hosp %i recovered_from_icu %i contact_work %i contact_school %i contact_home %i contact_community %i \n",
            step_time,
            t,
            num_infect,
            num_infectious,
            num_hosp,
            num_icu,
            num_dead,
            num_recovered,
            recovered_hosp,
            recovered_icu,
            num_contact_work,
            num_contact_school,
            num_contact_house,
            num_contact_commun);
    fflush(output_file);

    fprintf(output_HCW,
            "Walltime/timestep %6.2f Time %6.2f num_infected %i num_infectious %i num_in_hosp %i num_in_icu %i num_dead %i recovered_tot %i recovered_from_hosp %i recovered_from_icu %i contact_work %i contact_school %i contact_home %i contact_community %i total_healthcare %i \n",
            step_time,
            t,
            num_infect_HCW,
            num_infectious_HCW,
            num_hosp_HCW,
            num_icu_HCW,
            num_dead_HCW,
            num_recovered_HCW,
            num_recovered_hosp_HCW,
            num_recovered_icu_HCW,
            num_contact_work_HCW,
            num_contact_school_HCW,
            num_contact_house_HCW,
            num_contact_commun_HCW,
            num_HCW);
    fflush(output_HCW);

    for (file_count = 0; file_count < num_counties; file_count++) {
      fprintf(county_files[file_count],
              "Walltime/timestep %6.2f Time %6.2f num_infected %i num_infectious %i num_in_hosp %i num_in_icu %i num_dead %i recovered_tot %i recovered_from_hosp %i recovered_from_icu %i contact_work %i contact_school %i contact_home %i contact_community %i total_individuals %i \n",
              step_time,
              t,
              num_infect_county[file_count],
              num_infectious_county[file_count],
              num_hosp_county[file_count],
              num_icu_county[file_count],
              num_dead_county[file_count],
              num_recovered_county[file_count],
              num_recovered_hosp_county[file_count],
              num_recovered_icu_county[file_count],
              num_contact_work_county[file_count],
              num_contact_school_county[file_count],
              num_contact_house_county[file_count],
              num_contact_commun_county[file_count],
              county_size[file_count]);
      fflush(county_files[file_count]);
    }

    for (file_count = 0; file_count < 90; file_count += 5) {
      fprintf(age_files[file_count / 5],
              "Walltime/timestep %6.2f Timestep %6.2f num_infected %i num_infectious %i num_in_hosp %i num_in_icu %i num_dead %i recovered_tot %i recovered_from_hosp %i recovered_from_icu %i contact_work %i contact_school %i contact_home %i contact_community %i total_individuals %i \n",
              step_time,
              t,
              num_infect_age[file_count / 5],
              num_infectious_age[file_count / 5],
              num_hosp_age[file_count / 5],
              num_icu_age[file_count / 5],
              num_dead_age[file_count / 5],
              num_recovered_age[file_count / 5],
              num_recovered_hosp_age[file_count / 5],
              num_recovered_icu_age[file_count / 5],
              num_contact_work_age[file_count / 5],
              num_contact_school_age[file_count / 5],
              num_contact_house_age[file_count / 5],
              num_contact_commun_age[file_count / 5],
              age_distrib[file_count / 5]);
      fflush(age_files[file_count / 5]);
    }

    t = t + dt;
  }
  ret       = clock_gettime(CLOCK_MONOTONIC, &T2);
  step_time =
    ((double)T2.tv_sec + (double)T2.tv_nsec /
     nsdiv) - ((double)T1.tv_sec + (double)T1.tv_nsec / nsdiv);
  fprintf(output_file, "Total time %8.3f\n", step_time);
  printf("Total simulation time %8.3f\n", step_time);


    free(infect_kappa);

    free(commun_nom1);
    free(house_infect);

  for (i = 0; i < 6; ++i) {
    free(work_infect[i]);
  }
    free(HH);

  for (i = 0; i < num_households; ++i) {
    arrfree(per_HH_members[i]);
  }
    free(per_HH_members);

  for (i = 0; i < 6; ++i) {
    free(workplace_size[i]);
  }
    free(lat_city);
    free(long_city);
    free(city_size);
    free(county_size);
    free(age);
    free(county);
    free(job_status);
    free(workplace);
    free(infected);
    free(severe);
    free(tau);
    free(symptomatic);

    free(icu_list);
    free(hosp_list);
    free(infectious);
    free(sus_list);

    free(hosp_pop);
    free(icu_pop);
    free(recovered);
    free(dead);
    free(workplace_tmp);

  for (i = 0; i < num_counties; i++) {
    fclose(county_files[i]);
  }

  for (i = 0; i < 85; i += 5) {
    fclose(age_files[(i / 5)]);
  }
    fclose(output_file);
    fclose(lat_lon_out);
    fclose(output_HCW);
}
