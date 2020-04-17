#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "workplace_dist.h"

void workplace_dist(int * workplace, int * job_status, int ** job_status_county, int * city, int num_cities, int * county, int num_counties, int population, int * max_num_WP , int * hosp_num, FILE * stats) {

  int pp_class = 19; //Assumption of 15 children per class. Not used yet
  int pp_preschool = 53; //Assumption of 200 children per school.
  int pp_school = 220; //Assumption of 200 children per school.
  int pp_hospital = 120; //Assumption of 120 people per hospital.
  int pp_work = 15; //Assumption of 15 people per close work group.
  int i;
  int j;
  int num_workplaces[6][num_cities];
  memset(num_workplaces, 0, 6*num_cities*sizeof(int));
  int num_workplaces2[6];
  memset(num_workplaces2, 0, 6*sizeof(int));

  /* THIS WILL NOW BE BY CITY FOR SCHOOLS AND PRESCHOOLS. */
  for (i=0; i < num_cities; i++) {
    for (j=0; j<2; j++) {
      /* First only schools then workplaces */
      if (job_status_county[j][i]>0) {
	num_workplaces[j][i]=ceil(job_status_county[j][i]/(double)pp_preschool);
	num_workplaces2[j]+=ceil(job_status_county[j][i]/(double)pp_preschool);
	if (num_workplaces2[j]>*max_num_WP) {
	  *max_num_WP=num_workplaces2[j];
	}
      }
    }
    for (j=2; j<3; j++) {
      /* First only schools then workplaces */
      if (job_status_county[j][i]>0) {
	num_workplaces[j][i]=ceil(job_status_county[j][i]/(double)pp_school);
	num_workplaces2[j]+=ceil(job_status_county[j][i]/(double)pp_school);
	if (num_workplaces2[j]>*max_num_WP) {
	  *max_num_WP=num_workplaces2[j];
	}
      }
    }
  }

  // Broken down in case we want to do schools by municipality.
  for (i=0; i < num_counties; i++) {
    for (j=3; j<6; j++) {
      if (job_status_county[j][i]>0 && j==4) {
	num_workplaces[j][i]=ceil(job_status_county[j][i]/(double)pp_work);
	num_workplaces2[j]+=ceil(job_status_county[j][i]/(double)pp_work);
	if (num_workplaces2[j]>*max_num_WP) {
	  *max_num_WP=num_workplaces2[j];
	}
      } else if (job_status_county[j][i]>0 && j==3) {
	num_workplaces[j][i]=ceil(job_status_county[j][i]/(double)pp_school);
	num_workplaces2[j]+=ceil(job_status_county[j][i]/(double)pp_school);
	if (num_workplaces2[j]>*max_num_WP) {
	  *max_num_WP=num_workplaces2[j];
	}
      } else if (job_status_county[j][i]>0 && j==5) {
	num_workplaces[j][i]=ceil(job_status_county[j][i]/(double)pp_hospital);
	num_workplaces2[j]+=ceil(job_status_county[j][i]/(double)pp_hospital);
	if (num_workplaces2[j]>*max_num_WP) {
	  *max_num_WP=num_workplaces2[j];
	}
      }
    }
    // Need to use floor+1 equation because each county should have a hospital even if no one works there. */
    num_workplaces[5][i]=floor(job_status_county[5][i]/(double)pp_hospital)+1;
    num_workplaces2[5]+=floor(job_status_county[5][i]/(double)pp_hospital)+1;
    if (num_workplaces2[5]>*max_num_WP) {
      *max_num_WP=num_workplaces2[5];
    }
    hosp_num[i]=num_workplaces[5][i];
  }

  for (i=0; i < population; i++) {
    /* Try to minimize necessary memory by making job_numbers independent with job_status. */
    int prior_workplaces=0;
    if (job_status[i] < 3 && num_workplaces[job_status[i]][city[i]] > 0) {
      for (j=0; j < city[i]; j++) {
	prior_workplaces+=num_workplaces[job_status[i]][j];
      }
      workplace[i]=(int)(COV_rand() * (num_workplaces[job_status[i]][city[i]]))+prior_workplaces;
    } else if (job_status[i]>3 && (num_workplaces[job_status[i]][county[i]])>0) {
      for (j=0; j<county[i]; j++) {
	prior_workplaces+=num_workplaces[job_status[i]][j];
      }
      workplace[i]=(int)(COV_rand() * (num_workplaces[job_status[i]][county[i]]))+prior_workplaces;
    } else if ((num_workplaces[job_status[i]][county[i]])>0) {
      for (j=0; j<county[i]; j++) {
	prior_workplaces+=num_workplaces[job_status[i]][j];
      }
      workplace[i]=(int)(COV_rand() * (num_workplaces[job_status[i]][county[i]]))+prior_workplaces;
    } else {
      workplace[i]=0;
      printf("no workplaces %i %i %i \n", i, job_status[i], county[i]);
    }

  }
}
