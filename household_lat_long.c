
/* Puts households in certain locations based on population density
 * distribution.  Only really correct for full population but works for
 * smaller populations.  Biases towards smaller households for smaller
 * populations.  Each household is also fed into a locality based on the
 * shortest distance to the center of that locality for purposes of
 * school and workplace choice.
 */

#include <stdlib.h>
#include <math.h>

#include "common.h"
#include "distance.h"
#include "household_lat_long.h"

#include "stb_ds.h"

int extern sort_household_members, typical_max_HH_sz;

int cmp_int(void const* a_, void const* b_) {
	int a = *(int*)a_;
	int b = *(int*)b_;
	if (a != b) {
		return (a < b) ? -1 : 1;
	}
	return 0;
}


void household_lat_long(int num_households, int * HH, double * lat_city, double * long_city, int num_cities, int * city, int * county, int * city_county, int * city_size, int * county_size, int population, double * age, int** per_HH_members, char ** county_names, double * county_pop, double tot_pop_actual, FILE* stats, int **county_p, int *num_locale, double *lat_locale, double *lon_locale, double *pop_density_init_num, int **locale_to_HH, int *locale_HH, double land_pop_total_density) {

  /* Initialize helpers for population density */
  double tot_pop_density=0;

  /* Initialize helpers for household distribution */
  int* city_HH = malloc(num_households * sizeof(int));
  int* county_HH = malloc(num_households * sizeof(int));
  int county_num;
  int city_num;

  double dist1;
  double min_dist=1000000;
  int placement;
  int i, j;

  int tmp_county_count[21]={0};
  int tmp_county_density[21]={0};

  /* Fill up households. */
  int HH_count=0 ; // Counter
  int HH_person=0 ; // Counter
  int max_HH_size=0;

  /* Initially set all households to -1 */
  for (i=0; i < population; i++) {
    HH[i]=-1;
  }

  /* Save list of households in each locale. */
  int original_num_locale = *num_locale;
  int ** county_list;
  county_list = (int**)calloc(*num_locale,sizeof(int*));
  int * locale_count;
  locale_count = (int*)calloc(*num_locale,sizeof(int));
  int * locale_HH_count;
  locale_HH_count = (int*)calloc(*num_locale,sizeof(int));


  city_num = -1;
  while (( HH_count < num_households ) && (HH_count < *num_locale)) {
    double tmp_lat, tmp_lon;
    min_dist=1000000;

    tmp_lat = lat_locale[HH_count];
    tmp_lon = lon_locale[HH_count];
    /* Scale population density */
    pop_density_init_num[HH_count] = ceil(pop_density_init_num[HH_count] * population / land_pop_total_density);
    tot_pop_density += pop_density_init_num[HH_count];
    if (tot_pop_density > population || ceil(pop_density_init_num[HH_count+1] * population / land_pop_total_density) < 0.5) {
      *num_locale = HH_count + 1;
    }

    // Determine city of each population square.  Use city data to determine which schools students attend.  Workplaces are placed by county. //
    for (j=0; j<num_cities; j++) {
      dist1=distance(tmp_lat, tmp_lon, lat_city[j], long_city[j], 'K');
      if (dist1<min_dist) {
	min_dist=dist1;
	city_num=j;
      }
    }
    if (city_num < 0) {
      fprintf(stderr, "Error in household_lat_long: city_num < 0\n");
      exit(0);
    }

    county_num=city_county[city_num];
    tmp_county_count[county_num]++;
    tmp_county_density[county_num]+=pop_density_init_num[HH_count];

    /* Set up household */
    city_HH[HH_count]=city_num;
    county_HH[HH_count]=county_num;
    arrput(county_list[HH_count], HH_count);
    locale_HH[HH_count] = HH_count;
    arrput(locale_to_HH[HH_count], HH_count);
    locale_HH_count[HH_count]+=1;

    /* Allocate an adult to the household */
    while ( age[HH_person]<20 ) {
      HH_person++;
    }
    HH[HH_person]=HH_count;
    city[HH_person]=city_HH[HH_count];
    county[HH_person]=county_HH[HH_count];
    county_size[county[HH_person]]++;
    arrput(county_p[county[HH_person]], HH_person);
    city_size[city[HH_person]]++;
    arrput(per_HH_members[HH[HH_person]], HH_person);
    locale_count[HH_count]+=1;
    if (arrlen(per_HH_members[HH[HH_person]]) > max_HH_size) {
      max_HH_size=arrlen(per_HH_members[HH[HH_person]]);
    }
    HH_person++;
    HH_count++;
  }

  /* Uncomment to test population denity WRT county.
  // Test population based on method vs actual population density on county level.  Lines up fairly well. //
  for (i=0;i<21; i++) {
  printf("counties %i %s number of locales %i calculated density %f actual density %f\n", i, county_names[i], tmp_county_count[i], tmp_county_density[i]/(double)tot_pop_density, county_pop[i]/tot_pop_actual);
  fflush(stdout);
  }
  */


  placement=0; //Keeps track of household placement after first loop of locales.
  /* First add adult head of household to the rest of the households. */
  while ( HH_count < num_households ) {
    while ( age[HH_person]<20 ) {
      HH_person++;
    }
    /* Place another household in locales with population density greater than 2.2*(households in locale). */
    if (locale_HH_count[placement]*2.2 > pop_density_init_num[placement]) {
      placement++; /* continue with next locale */
    }
    city_HH[HH_count]=city_HH[placement];
    county_HH[HH_count]=county_HH[placement];
    arrput(county_list[placement], HH_count);
    locale_HH[HH_count] = placement;
    arrput(locale_to_HH[placement], HH_count);
    locale_HH_count[placement]+=1;

    /* Set up head of household. */
    HH[HH_person]=HH_count;
    city[HH_person]=city_HH[HH_count];
    county[HH_person]=county_HH[HH_count];
    county_size[county[HH_person]]++;
    arrput(county_p[county[HH_person]], HH_person);
    city_size[city[HH_person]]++;
    locale_count[placement]+=1;
    arrput(per_HH_members[HH[HH_person]], HH_person);
    if (arrlen(per_HH_members[HH[HH_person]]) > max_HH_size) {
      max_HH_size=arrlen(per_HH_members[HH[HH_person]]);
    }
    HH_person++;
    HH_count++;
  }

  placement = 0;
  /* Distribute remaining people randomly.  This could be changed to a distribution to more realistically reflect household size in the future. */
  for ( HH_person=0; HH_person<population ; HH_person++) {
    int tmp_HH, tmp_county_HH;
    if (HH[HH_person]==-1) {

      /* Place people in random households within locale until locale is full. */
      if (locale_count[placement] + 1  > pop_density_init_num[placement]) {
	placement++; /* continue with next locale */
      }

      /* Pick a random household in the locale. */
      int only_once;
      only_once = 0;
      tmp_county_HH=(int)(COV_rand() * locale_HH_count[placement]);
      tmp_HH=county_list[placement][tmp_county_HH];
      while (arrlen(per_HH_members[tmp_HH])+1 > typical_max_HH_sz) {
	tmp_county_HH++;
	if (tmp_county_HH >= locale_HH_count[placement]) {
	  if (only_once) {
	    if (placement + 1 >= *num_locale) {
	      printf("Bailing out, no more locale to put people in\n");
	      exit(0);
	    }
	    placement++;
	    tmp_county_HH=(int)(COV_rand() * locale_HH_count[placement]);
	    only_once = 0;
	  } else {
	    /* Start over from the beginning of county_list[placement] in the search for a unfilled household */
	    /* but only do that once */
	    tmp_county_HH = 0;
	    only_once = 1;
	  }
	}
	tmp_HH=county_list[placement][tmp_county_HH];
      }
      HH[HH_person]=tmp_HH;
      fflush(stdout);
      city[HH_person]=city_HH[HH[HH_person]];
      county[HH_person]=county_HH[HH[HH_person]];
      county_size[county[HH_person]]++;
      arrput(county_p[county[HH_person]], HH_person);
      city_size[city[HH_person]]++;
      locale_count[placement]+=1;
      arrput(per_HH_members[HH[HH_person]], HH_person);
      if (arrlen(per_HH_members[HH[HH_person]]) > max_HH_size) {
	max_HH_size=arrlen(per_HH_members[HH[HH_person]]);
      }
    }
  }

  printf("max_HH_size = %d\n", max_HH_size);
  fflush(stdout);

  if (sort_household_members) {
    for (size_t hh = 0; hh < num_households; ++hh) {
      int* set = per_HH_members[hh];
      if (set) {
	qsort(set, arrlen(set), sizeof(int), cmp_int);
      }
    }
  }


  /* Print household distribution */
  int *HH_dist_test;
  HH_dist_test = (int *)calloc(max_HH_size+1, sizeof(int));
  int *HH_county_test;
  HH_county_test = (int *)calloc(21, sizeof(int));
  for (i=0; i<num_households; i++) {
    HH_dist_test[arrlen(per_HH_members[i])]++;
    HH_county_test[county_HH[i]]++;
  }

  fprintf(stats, "Household distributions \n");
  for (i=0; i<=max_HH_size; i++) {
    fprintf(stats, "household_size %i fraction_of_households %f num_households %i total_households %i \n", i, HH_dist_test[i]/(double)num_households, HH_dist_test[i], num_households) ;
    fflush(stats);
  }
  free(HH_dist_test);

  fprintf(stats, "\n\n County distribution \n");
  for (i=0; i<21; i++) {
    fprintf(stats, "%s county %i population %i fraction %f actual %f num_households %i \n", county_names[i],i, county_size[i], county_size[i]/(double)population, county_pop[i]/tot_pop_actual, HH_county_test[i])  ;
    fflush(stats);
  }
  fprintf(stats, "\n\n");
  fflush(stats);
  free(HH_county_test);

  for (i = 0; i < original_num_locale; i++) {
    arrfree(county_list[i]);
  }
  free(county_list);
  free(locale_count);
  free(locale_HH_count);

  free(city_HH);
  free(county_HH);
}
