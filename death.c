
#include <math.h>
#include "common.h"
#include "death.h"

int death(double t, int num_infectious, int * infectious, double * tau, int * dead, int * icu_pop, int * hosp_pop, int * symptomatic, int num_dead, double * age, double dt, int * num_dead_county, int * num_dead_age, int * county, int * num_dead_HCW, int * job_status) {

  double fatal_in_icu=0.5;
  double fatal_symptomatic[]={0.0000161, 0.0000695, 0.000309, 0.000844, 0.00161, 0.00595, 0.0193, 0.0428, 0.078};
  int age_group;
  int infec_person;
  int i;

  for (i=0; i<num_infectious; i++) {
    infec_person=infectious[i];
    if (tau[infec_person]>=t-15 && tau[infec_person]<t-15+dt) {
      if (icu_pop[infec_person]==1) {
	icu_pop[infec_person]=0;
	if (COV_rand() < fatal_in_icu) {
	  dead[infec_person]=1;
	  num_dead++;
	  hosp_pop[infec_person]=0;
	  num_dead_county[county[infec_person]]++;
	  num_dead_age[(int)floor(age[infec_person]/5)]++;
	  if (job_status[infec_person]==5) {
	    *num_dead_HCW=*num_dead_HCW+1;
	  }
	}
      } else if (symptomatic[infec_person]==1) {
	age_group=floor(age[infec_person]/10);
	if (COV_rand() < fatal_symptomatic[age_group]) {
	  dead[infec_person]=1;
	  icu_pop[infec_person]=0;
	  hosp_pop[infec_person]=0;
	  num_dead++;
	  num_dead_county[county[infec_person]]++;
	  num_dead_age[(int)floor(age[infec_person]/5)]++;
	  if (job_status[infec_person]==5) {
	    *num_dead_HCW=*num_dead_HCW+1;
	  }
	}
      }
    }
  }
  return(num_dead);
}
