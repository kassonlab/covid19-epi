#include <math.h>
#include "common.h"
#include "hospital.h"

void hosp_entry(double t, int num_infectious, int * infectious, double * age, int * icu_pop, int * hosp_pop, int * symptomatic, double * tau, int * workplace_tmp, int * workplace_num, int * county, double dt) {

  double hosp[]={0.001, 0.003, 0.012, 0.032, 0.049, 0.102, 0.166, 0.243, 0.273}; //# From Ferguson 2020
  double icu[]={0.05, 0.05, 0.05, 0.05, 0.063, 0.122, 0.274, 0.432, 0.709} ; //# percent of hospitalized that need icu from Ferguson 2020.
  int age_group;
  int infec_person;
  int i;

  for (i=0; i<num_infectious; i++) {
    infec_person=infectious[i];
    if (tau[infec_person]>=t-10 && tau[infec_person]<t-10+dt && symptomatic[infec_person]==1) {
      age_group=floor(age[infec_person]/10);
      if (COV_rand() < hosp[age_group]) {
	if (COV_rand() < icu[age_group]) {
	  icu_pop[infec_person]=1;
	  hosp_pop[infec_person]=2;
	} else {
	  hosp_pop[infec_person]=1;
	}
	/* put them in a hospital */
	workplace_tmp[infec_person] = (int)(COV_rand() * workplace_num[county[infec_person]]);
      }
    }
  }
}


void hosp_release(double t, int num_hosp, int * hosp_list, double * tau, int * recovered, int * hosp_pop, int * num_recovered, int * recovered_hosp, int * recovered_icu, double dt, int * num_recovered_county, int * num_recovered_age, double * age, int * county, int * num_recovered_hosp_county, int * num_recovered_hosp_age, int * num_recovered_icu_county, int * num_recovered_icu_age, int * num_recovered_HCW, int * recovered_hosp_HCW, int * recovered_icu_HCW, int * job_status ) {

  int i;
  int infec_person;

  for (i=0; i<num_hosp; i++) {
    infec_person=hosp_list[i];
    /* Regular hospital patients get released after 8 days (18 days from infection).  ICU patients get released after 26 days but do not go into regular hospital until 20 days.  */
    if ((tau[infec_person]>=t-18) && (tau[infec_person]<t-18+dt) && (hosp_pop[infec_person]==1)) {
      hosp_pop[infec_person]=0;
      recovered[infec_person]=1;
      *num_recovered=*num_recovered+1;
      *recovered_hosp=*recovered_hosp+1;
      num_recovered_county[county[infec_person]]++;
      num_recovered_age[(int)floor(age[infec_person]/5)]++;
      num_recovered_hosp_county[county[infec_person]]++;
      num_recovered_hosp_age[(int)floor(age[infec_person]/5)]++;
      if (job_status[infec_person]==5) {
	*num_recovered_HCW=*num_recovered_HCW+1;
	*recovered_hosp_HCW=*recovered_hosp_HCW+1;
      }
    } else if ((tau[infec_person]>=t-25) && (tau[infec_person]<t-25+dt) && (hosp_pop[infec_person]==2)) {
      hosp_pop[infec_person]=0;
      recovered[infec_person]=1;
      *num_recovered=*num_recovered+1;
      *recovered_icu=*recovered_icu+1;
      num_recovered_county[county[infec_person]]++;
      num_recovered_age[(int)floor(age[infec_person]/5)]++;
      num_recovered_icu_county[county[infec_person]]++;
      num_recovered_icu_age[(int)floor(age[infec_person]/5)]++;
      if (job_status[infec_person]==5) {
	*num_recovered_HCW=*num_recovered_HCW+1;
	*recovered_icu_HCW=*recovered_icu_HCW+1;
      }
    }
  }
}
