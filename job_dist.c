#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "job_dist.h"

void job_dist(int * job_status, int ** job_status_city, double * age, int * county, int * city, int population, int num_cities, int num_counties, int * county_size, FILE * stats) {

  int i, j;


  /* All currently randomly placed based on county.  Would like to do per town but each town needs inhabitants. See commented section below for per city distribution of schools */
  for (i=0; i < population; i++) {
    if (age[i]<1 || age[i]>75) {
      job_status[i]=0;
      job_status_city[0][county[i]]++;
    } else if (age[i]>=1 && age[i]<3) {
      if (COV_rand() < 0.7800) {
	job_status[i]=1;
	job_status_city[1][city[i]]++;
      } else {
	job_status[i]=0;
	job_status_city[0][city[i]]++;
      }
    } else if (age[i]>=3 && age[i]<6) {
      if (COV_rand() < 0.9500) {
	job_status[i]=1;
	job_status_city[1][city[i]]++;
      } else {
	job_status[i]=0;
	job_status_city[0][city[i]]++;
      }
    } else if (age[i]>=6 && age[i]<15) {
      job_status[i]=2;
      job_status_city[2][city[i]]++;
    } else if (age[i]>=15 && age[i]<22) {
      job_status[i]=3;
      job_status_city[3][county[i]]++;
    } else if (age[i]>=22 && age[i]<=65) {
      if (COV_rand() < 0.773) {
	// 17.25% of workforce is in healthcare from OECD 2017 statstics.  Assume 1/4 of these are in hospitals.
	if (COV_rand() < 0.04325) {
	  job_status[i]=5;
	  job_status_city[5][county[i]]++; // Workplace is based on county, not city.
	} else {
	  job_status[i]=4;
	  job_status_city[4][county[i]]++; // Workplace is based on county, not city.
	}
      } else {
	job_status[i]=0;
	job_status_city[0][county[i]]++;
      }
    } else if (age[i]>=65 && age[i]<=75) {
      if (COV_rand() < 0.172) {
	job_status[i]=4;
	job_status_city[4][county[i]]++; // Workplace is based on county, not city.
      } else {
	job_status[i]=0;
	job_status_city[0][county[i]]++;
      }
    }
  }


  /* Print job distribution */
  int job_dist_test[6]={0};
  int *unemployed;
  unemployed = (int *)calloc(num_counties, sizeof(int));
  int *working_age;
  working_age = (int *)calloc(num_counties, sizeof(int));
  int **city_dist_test;
  city_dist_test = (int **)calloc(6, sizeof(int *));
  for (i = 0; i < 6; i++) {
    city_dist_test[i] = (int *)calloc(num_counties, sizeof(int));
  }
  for (i=0; i<population; i++) {
    job_dist_test[job_status[i]]++;
    city_dist_test[job_status[i]][county[i]]++;
    if (age[i]>22 && age[i]<=65) {
      if (job_status[i]==0) {
	unemployed[county[i]]++;
      }
      working_age[county[i]]++;
    }
  }

  fprintf(stats, "Job Distribution \n");
  for (j=0; j<num_counties; j++) {
    for (i=0; i<6; i++) {
      fprintf(stats, "job_status %i county %i fraction_of_jobs_total %f num_jobs_in_county %i num_unemployed %i fraction_unemployed  %f \n", i, j, job_dist_test[i]/(double)population, city_dist_test[i][j], unemployed[j], (double)unemployed[j]/working_age[j]) ;

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
