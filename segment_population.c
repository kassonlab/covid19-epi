
#include <stdio.h>
#include <math.h>

#include "segment_population.h"

void segment_population(int* num_sus, int* num_infectious, int* num_hosp, int* num_icu, int* num_sus_HCW, int* num_infectious_HCW, int* num_hosp_HCW, int* num_icu_HCW, int* infected, int* infectious, int* sus_list, int* hosp_list, int * hosp_pop, int * icu_pop, int* icu_list, double* tau, int population, double t, int * num_sus_county, int * num_infectious_county, int * num_infected_county, int * num_hosp_county, int * num_icu_county, int * num_sus_age, int * num_infectious_age, int * num_infected_age, int * num_hosp_age, int * num_icu_age, double * age, int * county, int print_loc, double * lat_locale, double * lon_locale, FILE * lat_lon_out, int *locale_HH, int *HH, int * job_status, int * recovered, int * dead) {

  int i;

  *num_sus=0;
  *num_infectious=0;
  *num_hosp=0;
  *num_icu=0;
  *num_sus_HCW=0;
  *num_infectious_HCW=0;
  *num_hosp_HCW=0;
  *num_icu_HCW=0;
  for (i=0; i<population; i++) {
    if (infected[i]==0) {
      sus_list[*num_sus]=i;
      *num_sus=*num_sus+1;
      num_sus_county[county[i]]++;
      num_sus_age[(int)floor(age[i]/5)]++;
      if (job_status[i]==5) {
	*num_sus_HCW=*num_sus_HCW+1;
      }
    } else if ((tau[i]<t-4.6) && recovered[i]==0 && dead[i]==0) {
      infectious[*num_infectious]=i;
      *num_infectious=*num_infectious+1;
      num_infectious_county[county[i]]++;
      num_infectious_age[(int)floor(age[i]/5)]++;
      if (print_loc==1) {
	fprintf(lat_lon_out, "Time %f Person %i lat %f lon %f \n", t, i, lat_locale[locale_HH[HH[i]]], lon_locale[locale_HH[HH[i]]]);
      }
      if (job_status[i]==5) {
	*num_infectious_HCW=*num_infectious_HCW+1;
      }
    } else {
      if (print_loc==1) {
	fprintf(lat_lon_out, "Time %f Person %i lat %f lon %f \n", t, i, lat_locale[locale_HH[HH[i]]], lon_locale[locale_HH[HH[i]]]);
      }
    }


    if (icu_pop[i]==1) {
      icu_list[*num_icu]=i;
      *num_icu=*num_icu+1;
      num_icu_county[county[i]]++;
      num_icu_age[(int)floor(age[i]/5)]++;
      if (job_status[i]==5) {
	*num_icu_HCW=*num_icu_HCW+1;
      }
    }
    if (hosp_pop[i]>0) {
      hosp_list[*num_hosp]=i;
      *num_hosp=*num_hosp+1;
      num_hosp_county[county[i]]++;
      num_hosp_age[(int)floor(age[i]/5)]++;
      if (job_status[i]==5) {
	*num_hosp_HCW=*num_hosp_HCW+1;
      }
    }

  }
  fprintf(lat_lon_out, "\n\n");
}
