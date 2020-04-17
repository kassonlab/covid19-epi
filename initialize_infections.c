
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "initialize_infections.h"


void initialize_infections(int * initial_infections, double * tau, int * infected, int * severe, int * symptomatic, int * county, int * num_infect, int num_counties, double symptomatic_per, int population, double dt, double t, double * lat_locale, double * lon_locale, int * num_infect_county, int * num_infect_age, double * age, int **county_p, int *county_size, int *locale_HH, int *HH, double *tmp_lat, double *tmp_lon) {

  int person_infected=0, county_person_inf = 0;
  int tmp_infect=0;
  double min_diff=1000;
  double diff_lat_lon=10;
  int i, j;

  for (i=0; i < num_counties; i++) {
    tmp_infect=0;
    while ((tmp_infect<initial_infections[i])) {
      int tmp_j=0;
      // Test up to population to see if we can find someone who fits a perviously determined cluster.  If not, leave this loop and pick a random person.
      county_person_inf = (int)(COV_rand() * county_size[i]);
      person_infected=county_p[i][county_person_inf];
      while (((county[person_infected]!=i) || (infected[person_infected]!=0) || min_diff>0.5) && tmp_j<county_size[i]) {
	/* pick first available person in the county to infect if the randomly choosen aren't free to pick */
	county_person_inf++;
	if (county_person_inf >= county_size[i]) {
	  county_person_inf = 0;
	}
	person_infected=county_p[i][county_person_inf];
	if (county[person_infected] != i) {
	  fprintf(stderr, "Error: random person is not in expected county %d, is in %d\n", i, county[person_infected]);
	}
	min_diff=1000;
	if (t<-13 || num_infect_county[i]==0 ) {
	  min_diff=0;
	} else {
	  for (j=0; j<*num_infect; j++) {
	    diff_lat_lon=(fabsf(tmp_lat[j]-lat_locale[locale_HH[HH[person_infected]]])+fabsf(tmp_lon[j]-lon_locale[locale_HH[HH[person_infected]]]));
	    if (diff_lat_lon<min_diff) {
	      min_diff=diff_lat_lon;
	    }
	  }
	}
	tmp_j++;
      }
      if (min_diff>1) {
	while ((county[person_infected]!=i) || (infected[person_infected]!=0)) {
	  county_person_inf++;
	  if (county_person_inf >= county_size[i]) {
	    county_person_inf = 0;
	  }
	  person_infected=county_p[i][county_person_inf];
	  if (county[person_infected] != i) {
	    fprintf(stderr, "Error: random person is not in expected county %d, is in %d\n", i, county[person_infected]);
	  }
	}
      }

      infected[person_infected]=1;
      severe[person_infected]=(int)(COV_rand() * 2);

      if (COV_rand() < symptomatic_per) {
	symptomatic[person_infected]=1;
      }
      tau[person_infected]=t;
      tmp_lat[*num_infect]=lat_locale[locale_HH[HH[person_infected]]];
      tmp_lon[*num_infect]=lon_locale[locale_HH[HH[person_infected]]];
      *num_infect=*num_infect+1;
      num_infect_county[i]++;
      num_infect_age[(int)(floor(age[person_infected]/5))]++;
      tmp_infect++;
    }

  }
}
