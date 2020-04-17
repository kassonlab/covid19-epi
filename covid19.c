#include <stdio.h>
#include <stddef.h>
#include <stdint.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

//#ifdef _MPI
#include <mpi.h>
//#endif

#include "common.h"

#include "age_dist.h"
#include "calc_infect.h"
#include "city_lat_long.h"
#include "death.h"
#include "hospital.h"
#include "household_lat_long.h"
#include "initialize_infections.h"
#include "job_dist.h"
#include "calc_kappa.h"
#include "segment_population.h"
#include "workplace_dist.h"

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

int typical_max_HH_sz = 7;

int sort_household_members = 0;

int full_kappa = 0;

double betac_scale = 8.4, betah_scale = 2.0, betaw_scale = 1.0, R0_scale=2.2;

#define pi 3.14159265358979323846


#if COV_GPU

void locale_infectious_loop(int num_locale, int population, int num_households, int num_infectious, int* infectious, double const* infect_kappa, double Ic, int* intervene, double t, double* tau, double* tauI, double* interIc, double dt, int* hosp_pop, int* icu_pop, double* lat_locale, double* lon_locale, int* locale_HH, int* HH, struct locale* locale_list, double omega, int* severe, double betac_scale, double* commun_nom1, double* fd_tot);

#endif

/***** COVID-19 infectious spread model *****
****** (C) 2020 Jasmine Gardner, PhD    *****
This program calculates the spread of COVID-19 across Sweden.

Use:
compile as : gcc -o covid covid19.c -lm
use as: ./covid -pop 100000 -sim_time 100 % for a population of 100,000 and a simulation time of 100 days.
******/


int main (int argc, char *argv[]) {

  /* Parameters available to change with command line */
  int population = 10000;  // Size of total population
  int tot_time=200; // Simulation time.
  int initial_infections[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; //Initial infections per county.
  double percent_infect=0.001 ; // Default 1% of population initially infected.
  double dt=1.00; // Time step.
  int interventions=0; // Value of interventions.
  double tauI_onset=1; //time after start of simulation that interventions for whole community take place.
  int print_lat_lon=0; // Choose whether to print latitude and longitude data of infected individuals for each time step.
  int ret;
  int i, j;
  char *initial_infect_filename = NULL;

  /* Timing variables */
  struct timespec T1, T2, t1, t2, t3, t4;
  double tt, step_time, nsdiv = 1000*1000*1000;

  /* Parse command-line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-pop")) population=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-sim_time")) tot_time=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-infect")) percent_infect=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dt")) dt=atof(argv[++i]);
    else if (!strcmp(argv[i],"-inter")) interventions=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-tauI")) tauI_onset=atof(argv[++i]);
    else if (!strcmp(argv[i],"-initial")) for (j=0; j<21; j++) initial_infections[j]=atof(argv[++i]);
    else if (!strcmp(argv[i],"-print_loc")) print_lat_lon=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-betac")) betac_scale=atof(argv[++i]);
    else if (!strcmp(argv[i],"-betah")) betah_scale=atof(argv[++i]);
    else if (!strcmp(argv[i],"-betaw")) betaw_scale=atof(argv[++i]);
    else if (!strcmp(argv[i],"-R0")) R0_scale=atof(argv[++i]);
    else if (!strcmp(argv[i],"-full_kappa")) full_kappa = 1;
    else if (!strcmp(argv[i],"-use_fixed_seed")) use_fixed_seed = 1;
    else if (!strcmp(argv[i],"-use_seed")) seed = atol(argv[++i]);
    else if (!strcmp(argv[i],"-sort_HH")) sort_household_members = 1;
    else if (!strcmp(argv[i],"-initial_infect_file")) initial_infect_filename = argv[++i];
  }

  /* print out population statistics to file */
  FILE *stats = fopen("stats.log", "w");
  /*Get current time and stats */
  time_t date;
  time(&date);

  fprintf(stats, "Date: %s \n", ctime(&date));
  fprintf(stats, "Population: %i \nSimulationTime: %i \ndt: %f \nInterventions: %i \nTau_onset: %f \n\n\n ", population, tot_time, dt, interventions, tauI_onset);

  ret = clock_gettime(CLOCK_MONOTONIC, &T1);

  double HH_size = 2.2 ; // Average household size from OCED
  int num_households = (int)(population/HH_size); // Number of households in population
  int * HH;  // Households
  HH = (int*)calloc(population,sizeof(int));
  int **per_HH_members;
  per_HH_members = (int**)calloc(num_households, sizeof(int*));

  /* Population information */
  /* Specific to Sweden.  Averaging population based on total population of Sweden regardless of population size. */
  double tot_pop=10327589.; //For use of averaging counties only.  NOTE: land scan tot is 10,098,554, cannot use more inhabitants than that.
  char * county_name[]={"Stockholm", "Uppsala", "Sodermanland", "Ostergotland", "Jonkoping", "Kronoberg", "Kalmar", "Gotland", "Blekinge", "Skane", "Halland", "Vastra Gotaland", "Varmland", "Orebro", "Vastmanland", "Dalarna", "Gavleborg", "Vasternorrland", "Jamtland", "Vasterbotten", "Norrbotten"};
  double pop_county[]={2377081, 383713, 297540, 465495, 363599, 201469, 245446, 59686, 159606, 1377827, 333848, 1725881, 282414, 304805, 275845, 287966, 287382, 245347, 130810, 271736, 250093};
  int num_counties=(sizeof(pop_county)/sizeof(pop_county[0]));
  int county_int[num_counties]; // Integer values for random probability distribution.
  memset(county_int, 0, num_counties*sizeof(int));
  double * lat_city; //latitude of city i.
  double * fd_tot; //additive kernel density function of person with respect to all other individuals
  fd_tot = (double*)calloc(population,sizeof(double));
  lat_city = (double*)calloc(2000,sizeof(double));
  double * long_city; // longitude of city i.
  long_city = (double*)calloc(2000,sizeof(double));
  int * city_size; // populatino of city i.
  city_size = (int*)calloc(2000,sizeof(int));
  int * county_size; // population of county i.
  county_size = (int*)calloc(num_counties,sizeof(int));
  char * city_names[2000]; //City name of city i.
  int * city; // city of person i by integer assignment.
  city = (int*)calloc(population,sizeof(int));
  int * city_county; // county of person i by integer assignment.
  city_county = (int*)calloc(2000,sizeof(int));
  int num_cities=0; // total number of cities


  /* City information for allocating schools. */

  /* Initialize and allocate arrays */
  double * age;  // Age of population
  age = (double*)calloc(population,sizeof(double));
  int * age_distrib;
  /* Distribution of ages in the generated population, people probably won't be older ehan 150 */
  age_distrib = (int*)calloc(150,sizeof(int));
  int * county;  // County of each inhabitant
  county = (int*)calloc(population,sizeof(int));
  int **county_p; /* List of persons per county */
  county_p = (int **)calloc(num_counties, sizeof(int *));
  int * job_status; // Type of job each person holds: 0-5 for no job, preschool, elementary school, highschool/college, job, and hospital, respectively.
  job_status = (int*)calloc(population,sizeof(int));
  int * workplace; // Workplace of each person.
  workplace = (int*)calloc(population,sizeof(int));

  int max_num_WP=0; // max workplaces per job_status for allocating array.

  /* Parameters for infections */
  double symptomatic_per=0.67; // percent of people who are symptomatic.
  int * infected; // 1 if person i has been infected, 0 otherwise
  infected = (int*)calloc(population,sizeof(int));
  int * severe; // 1 if person i has severe infection, 0 otherwise
  severe = (int*)calloc(population,sizeof(int));
  double * tau; // day that person i became infected
  tau = (double*)calloc(population,sizeof(double));
  int * symptomatic; // 1 if person i is symptomatic, 0 otherwise
  symptomatic = (int*)calloc(population,sizeof(int));

  /* Infection parameters.  Mostly counters for looping. */
  int num_sus=0; //Number susceptible.
  int num_sus_HCW=0; //Number susceptible healthcare workers.
  int num_sus_age[18]={0};
  int num_infect=0; //Number infected.
  int num_infect_HCW=0; //Number infected.
  int *num_sus_county;
  num_sus_county = (int *)calloc(num_counties, sizeof(int));
  int * num_infect_county;
  num_infect_county = (int*)calloc(num_counties, sizeof(int));
  int num_infect_age[18]={0};
  int num_icu=0; //Number in ICU.
  int num_icu_HCW=0; //Number in ICU.
  int *num_icu_county;
  num_icu_county = (int *)calloc(num_counties, sizeof(int));
  int num_icu_age[18]={0};
  int num_hosp=0; //Number in hospital.
  int num_hosp_HCW=0; //Number in hospital.
  int *num_hosp_county;
  num_hosp_county = (int *)calloc(num_counties, sizeof(int));
  int num_hosp_age[18]={0};
  int num_dead=0; //Number of deaths.
  int num_dead_HCW=0; //Number of deaths.
  int *num_dead_county;
  num_dead_county = (int *)calloc(num_counties, sizeof(int));
  int num_dead_age[18]={0};
  int num_recovered=0; //Number of people recovered.
  int num_recovered_HCW=0; //Number of people recovered.
  int *num_recovered_county;
  num_recovered_county = (int *)calloc(num_counties, sizeof(int));
  int num_recovered_age[18]={0};
  int *num_recovered_hosp_county;
  num_recovered_hosp_county = (int *)calloc(num_counties, sizeof(int));
  int num_recovered_hosp_age[18]={0};
  int *num_recovered_icu_county;
  num_recovered_icu_county = (int *)calloc(num_counties, sizeof(int));
  int num_recovered_icu_age[18]={0};
  int num_infectious=0; //Number infectious.
  int num_infectious_HCW=0; //Number infectious.
  int *num_infectious_county;
  num_infectious_county = (int *)calloc(num_counties, sizeof(int));
  int num_infectious_age[18]={0};
  int recovered_hosp=0; //Number of people recovered.
  int recovered_icu=0; //Number of people recovered.
  int num_recovered_hosp_HCW=0; //Number of people recovered.
  int num_recovered_icu_HCW=0; //Number of people recovered.


  /* This list hold indices for each category. */
  int * icu_list; //Indices in ICU.
  icu_list = (int*)calloc(population,sizeof(int));
  int * hosp_list; //Indices in hospital.
  hosp_list = (int*)calloc(population,sizeof(int));
  int * infectious; //Indices of infected.
  infectious = (int*)calloc(population,sizeof(int));
  int * sus_list; //Indices of susceptible.
  sus_list = (int*)calloc(population,sizeof(int));

  /* Parameters pertaining to the population.  These could be made into a struct later to look nicer. */
  int * hosp_pop; // 1 if person i in hospital, 0 otherwise
  hosp_pop = (int*)calloc(population,sizeof(int));
  int * icu_pop; // 1 if person i in ICU, 0 otherwise
  icu_pop = (int*)calloc(population,sizeof(int));
  int * recovered; // 1 if person i is recovered, 0 otherwise.
  recovered = (int*)calloc(population,sizeof(int));
  int * dead; // 1 if person i is dead, 0 otherwise
  dead = (int*)calloc(population,sizeof(int));
  int * workplace_tmp; // Hospital location when people go into hospital.
  workplace_tmp = (int*)calloc(population,sizeof(int));

  double community_nom=0; // For adding community infection.
  int file_count;

  int num_I=9;
  double Ic=1.0; //Intervention constant for community transmission.
  double (*Iw); //Intervention constant for workplace transmission.
  double Ih=1; //Intervention constant for household transmission.
  double complyI[num_I]; // percent of people who comply with intervention
  double tauI[num_I]; // time after infection that intervention takes place
  double interIc[num_I]; //Intervention constants for Ic
  double interIh[num_I]; //Intervention constants for Ih
  int personinter[num_I]; // Tells us whether the intervention needs to be calculated on a person to person basis or for the whole community.
  double Ihosp=0.25;  // Accounts for increased cleanliness and infection control at hospital.
  int * intervene; // 1 if person is currently undergoing interventions, 0 otherwise.
  intervene = (int*)calloc(population,sizeof(int));


  /* current recommendations are intervention 1, 3, and 8. */
  /**** Introduce Interventions: must include documentation for values *****/
  /* No interventions. */
  interIc[0]=1.0;
  double interIw0[6]={0, 1, 1, 1, 1, Ihosp};
  Iw=interIw0;
  interIh[0]=1.00;
  complyI[0]=1.0;
  tauI[0]=0;
  personinter[0]=1;

  /* Intervention 1: school closures, only highschools and colleges.  No school transmission for job_status 3. Assuming no increase in community transmission as students would be working online at the times they would be in college class.*/
  interIc[1]=1.25;
  double interIw1[6]={0, 1, 1, 0, 1, Ihosp};
  interIh[1]=1.50;
  complyI[1]=1.0;
  tauI[1]=0;
  personinter[1]=0;

  /* Intervention 2: school closures of all schools. No school transmission for job_status 1, 2, and 3, reduction of 5% in workplace interactions to account for parents becoming childcare.  Children have 50% increase in household transmission and 25% increase in community transmission .*/
  interIc[2]=1.25;
  double interIw2[6]={0, 0, 0, 0, 1.0, Ihosp};
  interIh[2]=1.50;
  complyI[2]=1.0;
  tauI[2]=0;
  personinter[2]=0;

  /* Intervention 3: Case isolation within household. 1 day after symptoms start, 90% comply, householdi contacts remain the same, 25% contact with community, no contact with
     school or workplace. */
  interIc[3]=0.25;
  double interIw3[6]={0, 0, 0, 0, 0, Ihosp};
  interIh[3]=1.0;
  complyI[3]=0.9;
  tauI[3]=6.1;
  personinter[3]=0;

  /* Intervention 4: Case isolation of entire household if one member becomes sick.  Same as case isoloation of single person but now includes all in household.  90% of symptomatic comply and 70% household members comply. */
  interIc[4]=0.25;
  double interIw4[6]={0,0,0,0,0, Ihosp};
  interIh[4]=1.5;
  complyI[4]=0.7;
  tauI[4]=0.0;
  personinter[4]=0;

  /* Intervention 5: Case isolation of entire household if one member becomes sick.  This adds for the case of a quarantined household member getting ill.  tauI=0 */
  interIc[5]=0.25;
  double interIw5[6]={0,0,0,0,0, Ihosp};
  interIh[5]=1.0;
  complyI[5]=0.9;
  tauI[5]=0.0;
  personinter[5]=0;

  /* Intervention 6: social distancing with school closure.  Community contacts decrease by 75%, household comntact increase by 25%, 70% compliance.  essential buisnesses stay open, 75% reduction in workplace transmission. NOTE: similar to below except with minimized social interaction. */
  interIc[6]=0.25;
  double interIw6[6]={0,0,0,0,0.25, Ihosp};
  interIh[6]=1.50;
  complyI[6]=1.0;
  tauI[6]=0;
  personinter[6]=1;

  /* Intervention 7: school closures of all schools and non-essential businesses. No school transmission for job_status 1, 2, and 3, reduction of 75% workplace interactions.  50% increase in household transmission and 50% increase in community transmission .*/

  interIc[7]=1.50;
  double interIw7[6]={0, 0, 0, 0, 0.25, Ihosp};
  interIh[7]=1.50;
  complyI[7]=1.0;
  tauI[7]=0;
  personinter[7]=1;

  /* Intervention 8: Social distancing of people over 70 years old. Reduction of 75% workplace interactions. decrease of 75% of community contacts, household contacts increases 25%. 80% comply*/

  interIc[8]=0.25;
  double interIw8[6]={0, 0, 0, 0, 0.25, Ihosp};
  interIh[8]=1.25;
  complyI[8]=0.8;
  tauI[8]=0;
  personinter[8]=0;

  /* Make interIw array.*/
  double *interIw[9]={interIw0, interIw1, interIw2, interIw3, interIw4, interIw5, interIw6, interIw7, interIw8};

  /**** Set random number generator seed. ****/
  COV_init_rand();

  /* Initialize age distribution */
  age_dist(age, population, stats, age_distrib);
  city_lat_long(&num_cities,  lat_city,  long_city, city_names, city_county, county_name, num_counties) ;

  /* Parse land_scan file to get population density.  */
  double *lat_locale = NULL, *lon_locale = NULL, *pop_density_init_num = NULL;
  //int num_locale = 0, max_locale = 0;
  double tmp_lat, tmp_lon, pop_den, land_pop_total_density = 0;
  FILE* lat_long = fopen("land_pop_sorted.txt", "r"); // Sorted land population in descending order.  Important when we don't have complete population.
  while ((ret = fscanf(lat_long, "%lf%*c%lf%*c%lf", &tmp_lon, &tmp_lat, &pop_den)) == 3) {
    if (num_locale + 1 > max_locale) {
      //max_locale += 10;
      lat_locale = (double *)realloc(lat_locale, (max_locale+10) * sizeof(double));
      lon_locale = (double *)realloc(lon_locale, (max_locale+10) * sizeof(double));
      pop_density_init_num = (double *)realloc(pop_density_init_num, (max_locale+10) * sizeof(double));
    }
    lat_locale[num_locale] = tmp_lat;
    lon_locale[num_locale] = tmp_lon;
    pop_density_init_num[num_locale] = pop_den;
    add_locale(tmp_lat, tmp_lon, pop_den);
    land_pop_total_density += pop_den;
    num_locale++;
  }
  fclose(lat_long);

  int original_num_locale = num_locale;
  int **locale_to_HH;
  int *locale_HH;
  locale_to_HH = (int **)calloc(num_locale, sizeof(int *));
  locale_HH = (int *)calloc(num_households, sizeof(int));

  /* Initialize households */
  household_lat_long( num_households,  HH,  lat_city, long_city, num_cities, city, county, city_county, city_size, county_size, population, age, per_HH_members, county_name, pop_county, tot_pop, stats, county_p, &num_locale, lat_locale, lon_locale, pop_density_init_num, locale_to_HH, locale_HH, land_pop_total_density) ;

  free(city_county);


  /* Open files */
  char * file_beg="county_";
  char * file_end=".log";
  char file_name[1000];
  FILE** county_files = malloc(num_counties * sizeof(FILE*));
  for (i=0; i<num_counties; i++) {
    strcpy(file_name, "");
    strcat(file_name, file_beg);
    strcat(file_name, county_name[i]);
    strcat(file_name, file_end);
    county_files[i]=fopen(file_name, "w");
  }
  file_beg="age_";
  file_end=".log";
  char str;
  FILE** age_files = malloc(sizeof(FILE*) * 18);
  for (i=0; i<90; i+=5) {
    str=i;
    strcpy(file_name, "");
    strcat(file_name, file_beg);
    sprintf(file_name, "%s%i", file_name, i);
    strcat(file_name, file_end);
    age_files[(i/5)]=fopen(file_name, "w");
  }
  FILE * output_file = fopen("covid19_spread.dat", "w");
  FILE * output_HCW = fopen("healthcare_workers.log", "w");
  FILE * lat_lon_out = fopen("lat_lon.dat", "w");

  // Infections are randomly placed based on number of initial infections.  //
  // Includes infections from t=-11 to t=-1.
  // Percent per county taken from C19.se infections as of 2020/3/25.
  // Initial infections calculated from population admitted to intensive care per day from 2020/3/14 to 2020/3/24.
  double initial_per[21]={0.4234, 0.0404, 0.0336, 0.0843, 0.0257, 0.0079, 0.0071, 0.0020, 0.00475, 0.0973, 0.0261, 0.1088, 0.0178, 0.0230, 0.0115, 0.0158, 0.0127, 0.0075, 0.0233, 0.0131, 0.0139};
  /***** THIS IS THE REAL INITIALIZATION ARRAY, based on ICU numbers, day 0 is 3/26 ******/
  //	int initialize_base[15]={1667, 4231, 4181, 4407, 3051, 1808, 2599, 1469, 1695, 339, 678, 791, 678, 339, 113};
  /***** THIS IS THE REAL INITIALIZATION ARRAY, based on ICU numbers, day 0 is 3/21 ******/
  int initialize_base[15]={3955, 4068, 5198, 3955, 3616, 4633, 4633, 4859, 5085, 3051, 1921, 2712, 1469, 1695, 452};
  int *initialize = initialize_base;
  double tmp_t;

  if (initial_infect_filename != NULL) {
    initialize = (int *)calloc(15, sizeof(int));
    FILE *iif = fopen(initial_infect_filename, "r");
    for (i = 0; i < 15; i++) {
      fscanf(iif, "%d", &initialize[i]);
    }
    fclose(iif);
  }
  fprintf(stats, "Initial Infections by county \n");
  fflush(stats);
  double *tmp_lat_v = (double*)calloc(population,sizeof(double));
  double *tmp_lon_v = (double*)calloc(population,sizeof(double));
  for ( tmp_t=-14; tmp_t<=0; tmp_t++ ) {
    int l=-tmp_t;
    for ( j=0; j<21; j++ ) {
      initial_infections[j]=initial_per[j]*initialize[l];
      fprintf(stats, "time %f county %i initial_infections %i fraction_of_infections %f total_intialized %i\n", tmp_t, j, initial_infections[j], initial_per[j], initialize[l]);
    }
    fprintf(stats, "\n\n");
    fflush(stats);
    /* Randomly assign initial infections */
    initialize_infections( initial_infections,  tau,  infected,  severe,  symptomatic,  county,  &num_infect,  num_counties,  symptomatic_per,  population, dt, tmp_t, lat_locale, lon_locale, num_infect_county, num_infect_age, age, county_p, county_size, locale_HH, HH, tmp_lat_v, tmp_lon_v) ;
  }
  if (initial_infect_filename != NULL) {
    free(initialize);
  }
  free(tmp_lat_v);
  free(tmp_lon_v);
  fflush(stats);
  printf("All infections initialized\n");
  fflush(stdout);

  for (size_t i = 0; i < num_counties; ++i) {
    arrfree(county_p[i]);
  }
  free(county_p);
  county_p = NULL;

  int ** job_status_county; // Jobs per county, or city for schools
  job_status_county = (int**)calloc(6,sizeof(int*));
  for (i=0;i<6;i++) job_status_county[i] = (int*)calloc(num_cities,sizeof(int)) ;
  /* Initialize job/school status */
  job_dist(job_status, job_status_county, age, county, city, population, num_cities, num_counties, county_size, stats);

  int num_HCW=0;
  //Get number of healthcare workers //
  for (i=0; i<population; i++) {
    if (job_status[i]==5) {
      num_HCW++;
    }
  }

  int hosp_num[num_counties]; //Number of hospitals per county
  memset(hosp_num, 0, num_counties);
  /* Initialize workplace/school */
  workplace_dist(workplace, job_status, job_status_county, city, num_cities, county, num_counties, population, &max_num_WP , hosp_num, stats);

  free(city);

  /* Get size of each workplace as an array.  Cannot allocate array until max_num_WP is known. */
  int* workplace_size[6];
  for (size_t i = 0; i < 6; ++i) {
    workplace_size[i] = (int*)calloc(max_num_WP, sizeof(double));
  }
  for (i=0; i < population; i++) {
    workplace_size[job_status[i]][workplace[i]]++;
  }


  for (j=1; j<6; j++) {
    double avg_work_size=0, num_WP=0;
    for (i=0; i<max_num_WP; i++) {
      if (workplace_size[j][i]>0) {
	avg_work_size+=workplace_size[j][i];
	num_WP++;
      }
    }
    fprintf(stats, "Job_status %i avg_WP_size %f \n", j, avg_work_size/num_WP);
  }
  fclose(stats);


  printf("Starting density kernel calculations\n");
  fflush(stdout);

  ret = clock_gettime(CLOCK_MONOTONIC, &t1);
  int *npl = (int *)malloc(num_locale * sizeof(int));
#ifdef _OPENMP
#pragma omp parallel for private(i) default(shared)
#endif
  for (i = 0; i < num_locale; i++) {
    int npi; /* number of persons in locale i */
    npi = 0;
    for (int hh = 0; hh < arrlen(locale_to_HH[i]); hh++) {
      npi += arrlen(per_HH_members[locale_to_HH[i][hh]]);
    }
    npl[i] = npi;
  }

  /* Precalculate total density kernel function for each individual */
  for (i=0; i<num_locale; i++) {
    double itmp_fd;
    itmp_fd = 0;
#ifdef _OPENMP
#pragma omp parallel for private(j) default(shared) reduction(+:itmp_fd)
#endif
    for (j=i+1; j<num_locale; j++) {
      double d;
      double tmp_fd;
#if !defined(USE_LOCALE_DISTANCE)
      d=distance(lat_locale[i], lon_locale[i], lat_locale[j], lon_locale[j], 'K');
#else
      d = locale_distance(locale_list[i], locale_list[j]);
#endif
      tmp_fd = 1/(1+pow((d/4), 3)); //kernel density function as parameterized for GB.
      itmp_fd += tmp_fd * npl[j];
      fd_tot[j] += tmp_fd * npl[i];
    }
    fd_tot[i] += itmp_fd + npl[i] - 1;
  }
  free(npl);
  ret = clock_gettime(CLOCK_MONOTONIC, &t2);
  tt = ((double)t2.tv_sec + (double)t2.tv_nsec/nsdiv) - ((double)t1.tv_sec + (double)t1.tv_nsec/nsdiv);
  printf("Done with density kernel calculations in %5.2f\n", tt);
  fflush(stdout);

  for (size_t i = 0; i < original_num_locale; ++i) {
    arrfree(locale_to_HH[i]);
  }
  free(locale_to_HH);


  /* Initialize constants */
  double alpha=0.8 ; // From Ferguson Nature 2006
  double omega=2 ; // From Ferguson Nature 2006
  // #Leaving out rho from Ferguson 2006.  This is a measure of how infectious person is.  For now, we will assume all people are the same.

  /* precalculate kappa */
  int count_kappa_vals=ceil(30/dt);
  double kappa_t=0;
  double * kappa_vals;
  double tau1=0;
  kappa_vals = (double*)calloc(count_kappa_vals,sizeof(double));
  for (i=1; i<count_kappa_vals; i++) {
    kappa_t=i*dt-4.6;
    if (kappa_t <= 0) {
      kappa_vals[i] = 0;
    } else {
      tmp_t=(log(kappa_t)+0.72)/1.8;
      kappa_vals[i]=exp(-0.5*pow(tmp_t,2))/((kappa_t)*1.8*sqrt(2*pi));
    }
  }

  int contact_commun=0;
  int contact_work=0;
  int contact_school=0;
  int contact_house=0;
  int num_contact_commun=0;
  int num_contact_work=0;
  int num_contact_school=0;
  int num_contact_house=0;
  int num_contact_commun_HCW=0;
  int num_contact_work_HCW=0;
  int num_contact_school_HCW=0;
  int num_contact_house_HCW=0;
  int num_contact_commun_county[num_counties];
  int num_contact_commun_age[18];
  int num_contact_work_county[num_counties];
  int num_contact_work_age[18];
  int num_contact_school_county[num_counties];
  int num_contact_school_age[18];
  int num_contact_house_county[num_counties];
  int num_contact_house_age[18];
  for (i=0; i<num_counties; i++) {
    num_dead_county[i]=0;
    num_recovered_county[i]=0;
    num_recovered_hosp_county[i]=0;
    num_recovered_icu_county[i]=0;
  }
  for (i=0; i<18; i++) {
    num_dead_age[i]=0;
    num_recovered_age[i]=0;
    num_recovered_hosp_age[i]=0;
    num_recovered_icu_age[i]=0;
  }

  double *commun_nom1 = (double*)malloc(num_locale * sizeof(double));
  double* work_infect[6];
  for (size_t i = 0; i < 6; ++i) {
    work_infect[i] = (double*)malloc(max_num_WP * sizeof(double));
  }
  double *house_infect = (double *)malloc(num_households * sizeof(double));


  /* NEED TO INITIALIZE PREVIOUS INFECTIONS TO RECOVERED HOSP AND ICU */
  Ic=1;
  double t=0;
  double time_step=0;

  /* Segment population into infectious, susceptible, hospitalized, and icu */
  segment_population( &num_sus,  &num_infectious,  &num_hosp,  &num_icu,  &num_sus_HCW,  &num_infectious_HCW,  &num_hosp_HCW,  &num_icu_HCW, infected,  infectious,  sus_list,  hosp_list,  hosp_pop,  icu_pop,  icu_list,  tau,  population,  t,  num_sus_county,  num_infectious_county,  num_infect_county,  num_hosp_county,  num_icu_county,  num_sus_age,  num_infectious_age,  num_infect_age,  num_hosp_age,  num_icu_age,  age,  county, print_lat_lon,  lat_locale,  lon_locale,  lat_lon_out, locale_HH, HH, job_status, recovered, dead);

  // After 5 days of symptoms, people are randomly put into hospital based on age distributions. //
  hosp_entry(t, num_infectious,  infectious,  age,  icu_pop,  hosp_pop,  symptomatic, tau, workplace_tmp, hosp_num, county, dt) ;

  // Release from hospital //
  hosp_release(t, num_hosp,  hosp_list,  tau,  recovered, hosp_pop, &num_recovered, &recovered_hosp, &recovered_icu, dt, num_recovered_county, num_recovered_age, age, county, num_recovered_hosp_county, num_recovered_hosp_age, num_recovered_icu_county, num_recovered_icu_age, &num_recovered, &recovered_hosp, &recovered_icu, job_status ) ;


  // 15 days after onset of symptoms, people randomly die based on age based percentages.  ICU patients are removed from ICU into regular hospital. //

  num_dead = death(t, num_infectious,  infectious,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt, num_dead_county, num_dead_age, county, &num_dead_HCW, job_status) ;
  num_dead = death(t, num_icu,  icu_list,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt, num_dead_county, num_dead_age, county, &num_dead_HCW, job_status) ;
  num_dead = death(t, num_hosp,  hosp_list,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt, num_dead_county, num_dead_age, county, &num_dead_HCW, job_status) ;


  // Recovered after 11 days (6 days of symptoms) if not in hospital/ICU. //
  for (i=0; i<num_infectious; i++) {
    int infec_person;
    infec_person=infectious[i];
    if ((tau[infec_person]<(t-11)+dt) && (hosp_pop[infec_person]==0) && (icu_pop[infec_person]==0)) {
      recovered[infec_person]=1;
      num_recovered++;
      if (job_status[infec_person]==5) {
	num_recovered_HCW++;
      }
      num_recovered_county[county[infec_person]]++;
      num_recovered_age[(int)floor(age[infec_person]/5)]++;

    }
  }

  /* Introduce overall community interventions. */
  if ( interventions == 0 ) {
    Ic=interIc[interventions];
    Ih=interIh[interventions];
    Iw=interIw[interventions];
  } else if ( interventions > 0 )  {
    for (i=0; i<population; i++) {
      /* high school and university closures */
      if (age[i]>=15 && age[i]<22) {
	intervene[i]=1;
      } else if ( age[i]>=70 && COV_rand()<complyI[8] ) {
	intervene[i]=8;
      }
    }
  }

  ret = clock_gettime(CLOCK_MONOTONIC, &T2);
  tt = ((double)T2.tv_sec + (double)T2.tv_nsec/nsdiv) - ((double)T1.tv_sec + (double)T1.tv_nsec/nsdiv);
  printf("Total initialization time %5.2f\n", tt);

  /* Start simulation */
  printf("Starting simulation\n");
  fflush(stdout);

  /* Precalculated kappa for infected, recalculated each timestep */
  double *infect_kappa = NULL;
  int sz_infect_kappa = 0;

  ret = clock_gettime(CLOCK_MONOTONIC, &T1);
  for (time_step=0; time_step<(tot_time/dt); time_step++) {
    ret = clock_gettime(CLOCK_MONOTONIC, &t1);
    printf("Timestep %5.2f ", t);
    fflush(stdout);

    /* count origin of contacts */
    num_contact_commun=0;
    num_contact_work=0;
    num_contact_school=0;
    num_contact_house=0;
    num_contact_commun_HCW=0;
    num_contact_work_HCW=0;
    num_contact_school_HCW=0;
    num_contact_house_HCW=0;
    for (i=0; i<num_counties; i++) {
      num_contact_commun_county[i]=0;
      num_contact_work_county[i]=0;
      num_contact_school_county[i]=0;
      num_contact_house_county[i]=0;
      num_sus_county[i]=0;
      num_hosp_county[i]=0;
      num_icu_county[i]=0;
      num_infectious_county[i]=0;
    }
    for (i=0; i<18; i++) {
      num_contact_house_age[i]=0;
      num_contact_school_age[i]=0;
      num_contact_work_age[i]=0;
      num_contact_commun_age[i]=0;
      num_sus_age[i]=0;
      num_hosp_age[i]=0;
      num_icu_age[i]=0;
      num_infectious_age[i]=0;
    }

    if (interventions == 3 && t >= tauI_onset && t <= tauI_onset+dt) {
      for ( i=0; i<population; i++ ) {
	if (age[i]>1 && age[i]<15) {
	  intervene[i]=2;
	}
      }
    } else if (interventions == 4 && t >= tauI_onset && t <= tauI_onset+dt) {
      Ic=interIc[7];
      Ih=interIh[7];
      Iw=interIw[7];
    } else if (interventions == 5 && t >= tauI_onset && t <= tauI_onset+dt) {
      /* if not complying, have same interactions as type 7 */
      Ic=interIc[7];
      Ih=interIh[7];
      Iw=interIw[7];
      for ( i=0; i<population; i++ ) {
	if ( COV_rand()<0.9 ) {
	  intervene[i]=6;
	}
      }
    }

    /* Segment population into infectious, susceptible, hospitalized, and icu */
    segment_population( &num_sus,  &num_infectious,  &num_hosp,  &num_icu,  &num_sus_HCW,  &num_infectious_HCW,  &num_hosp_HCW,  &num_icu_HCW, infected,  infectious,  sus_list,  hosp_list,  hosp_pop,  icu_pop,  icu_list,  tau,  population,  t,  num_sus_county,  num_infectious_county,  num_infect_county,  num_hosp_county,  num_icu_county,  num_sus_age,  num_infectious_age,  num_infect_age,  num_hosp_age,  num_icu_age,  age,  county, print_lat_lon,  lat_locale,  lon_locale,  lat_lon_out, locale_HH, HH, job_status, recovered, dead);


    memset(commun_nom1, 0, num_locale*sizeof(double));
    for (i = 0; i < 6; ++i) {
      memset(work_infect[i], 0, max_num_WP*sizeof(double));
    }
    memset(house_infect, 0, num_households*sizeof(double));

    if (num_infectious > sz_infect_kappa) {
      sz_infect_kappa = num_infectious;
      infect_kappa = (double *)realloc(infect_kappa, sz_infect_kappa * sizeof(double));
    }
    for(i = 0; i < num_infectious; i++) {
      int infec_person;
      infec_person=infectious[i];
      infect_kappa[i] = calc_kappa( t,  tau[infec_person], symptomatic[infec_person], dt, kappa_vals, hosp_pop[infec_person], icu_pop[infec_person]);
    }

#if COV_GPU

    locale_infectious_loop(num_locale, population, num_households, num_infectious, infectious, infect_kappa, Ic, intervene, t, tau, tauI, interIc, dt, hosp_pop, icu_pop, lat_locale, lon_locale, locale_HH, HH, locale_list, omega, severe, betac_scale, commun_nom1, fd_tot);

#else // COV_GPU

#ifdef _OPENMP
#pragma omp parallel for private(j, i) default(shared)
#endif
    for (j=0; j<num_locale; j++) {
      double tmp_comm_inf=0;
      for (i=0; i<num_infectious; i++) {
	int infec_person; //Counter for infected person.
	double kappa; // #Infectiousness
	double tIc;
	infec_person=infectious[i];
	tIc = Ic;
	if ( intervene[infec_person] > 0 && t>tau[infec_person]+tauI[intervene[infec_person]]) {
	  tIc=interIc[intervene[infec_person]];
	}
	kappa = infect_kappa[i];

	if (hosp_pop[infec_person]==0) {
	  double d; //distance between people.
	  // Community transmission //
#if !defined(USE_LOCALE_DISTANCE)
	  d=distance(lat_locale[j], lon_locale[j], lat_locale[locale_HH[HH[infec_person]]], lon_locale[locale_HH[HH[infec_person]]], 'K');
#else
	  d=locale_distance(locale_list[j], locale_list[locale_HH[HH[infec_person]]]);
#endif
	  tmp_comm_inf+=tIc*calc_community_infect( kappa, omega, severe[infec_person], d);
	}
      }
      commun_nom1[j]=tmp_comm_inf/fd_tot[j];
    }

#endif // COV_GPU

    for (i=0; i<num_infectious; i++) {
      int infec_person; //Counter for infected person.
      double kappa; // #Infectiousness
      double tIw, tIh;
      int age_group=0;
      int tmp_job_stat=0;
      double tmp_work_inf=0, tmp_house_inf=0;
      infec_person=infectious[i];
      tIh = Ih;
      tIw = Iw[job_status[infec_person]];
      if ( intervene[infec_person] > 0 && t>tau[infec_person]+tauI[intervene[infec_person]]) {
	tIw=interIw[intervene[infec_person]][job_status[infec_person]];
	tIh=interIh[intervene[infec_person]];
      }
      kappa = calc_kappa( t,  tau[infec_person], symptomatic[infec_person], dt, kappa_vals, hosp_pop[infec_person], icu_pop[infec_person]);

      if (hosp_pop[infec_person]==0) {
	tmp_job_stat=job_status[infec_person];
	// Workplace/School transmission: People must be in same workplace and job type. //
	tmp_work_inf=calc_workplace_infect(tmp_job_stat, kappa, omega, workplace_size[tmp_job_stat][workplace[infec_person]], severe[infec_person], tIw) ;
      } else {
	/* In hospital, only have interaction with hospital workers and half interaction with family (household). */
	tmp_job_stat=5;
	tIw = Iw[tmp_job_stat];
	tIh=0.25;
	tmp_work_inf=calc_workplace_infect(tmp_job_stat, kappa, omega, workplace_size[tmp_job_stat][workplace_tmp[infec_person]], severe[infec_person], tIw) ;
      }
      work_infect[tmp_job_stat][workplace[infec_person]]+=tmp_work_inf;
      // Household transmission //
      tmp_house_inf=tIh*calc_household_infect(kappa, omega, arrlen(per_HH_members[HH[infec_person]]), alpha, severe[infec_person]);
      house_infect[HH[infec_person]]+=tmp_house_inf;
    }


    //#### Only Susceptible people can get the virus and infected people spread it.
    for (i=0; i<num_sus; i++) {
      int sus_person; //Counter for susceptible person.
      int age_group;
      double infect; //Infectiousness
      double infect_prob; // Infectious probability
      double tIw, tIh, tIc;
      sus_person=sus_list[i];
      int contact_work=0;
      int contact_commun=0;
      int contact_house=0;
      int contact_school=0;

      tIh = Ih;
      tIc=Ic;
      tIw = Iw[job_status[sus_person]];
      if ( intervene[sus_person] > 0 && t>tau[sus_person]+tauI[intervene[sus_person]]) {
	tIw=interIw[intervene[sus_person]][job_status[sus_person]];
	tIh=interIh[intervene[sus_person]];
	tIc=interIc[intervene[sus_person]];
      }
      age_group=floor(age[sus_person]/5);
      double zeta[]={0.1, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.75, 0.50, 0.25, 0.25, 0.25} ; //   # Travel related parameter for community transmission. Ferguson Nature 2006
      infect = tIh*house_infect[HH[sus_person]];
      infect += tIw*work_infect[job_status[sus_person]][workplace[sus_person]];
      infect += tIc*zeta[age_group]*commun_nom1[locale_HH[HH[sus_person]]];

      //### Probability of being infected ####
      infect_prob=(1-exp(-infect*dt));
      //### Monte carlo type prediction ###
      if (COV_rand() < infect_prob) {
	infected[sus_person]=1;
	tau[sus_person]=t;
	severe[sus_person]=round(COV_rand());

	if (COV_rand() < symptomatic_per) {
	  symptomatic[sus_person]=1;
	}

	if (tIh*house_infect[HH[sus_person]]>0) {
	  num_contact_house++;
	  if (job_status[sus_person]==5) {
	    num_contact_house_HCW++;
	  }
	  num_contact_house_county[county[sus_person]]++;
	  num_contact_house_age[(int)floor(age[sus_person]/5)]++;
	} else if (tIw*work_infect[job_status[sus_person]][workplace[sus_person]] > 0) {
	  if (job_status[sus_person]<4) {
	    num_contact_school++;
	    num_contact_school_county[county[sus_person]]++;
	    num_contact_school_age[(int)floor(age[sus_person]/5)]++;
	  } else {
	    num_contact_work++;
	    if (job_status[sus_person]==5) {
	      num_contact_work_HCW++;
	    }
	    num_contact_work_county[county[sus_person]]++;
	    num_contact_work_age[(int)floor(age[sus_person]/5)]++;
	  }
	} else if (tIc * zeta[age_group]*commun_nom1[locale_HH[HH[sus_person]]] > 0) {
	  num_contact_commun++;
	  if (job_status[sus_person]==5) {
	    num_contact_commun_HCW++;
	  }
	  num_contact_commun_county[county[sus_person]]++;
	  num_contact_commun_age[(int)floor(age[sus_person]/5)]++;
	}
	num_infect++;
	num_infect_county[county[sus_person]]++;
	num_infect_age[(int)floor(age[sus_person]/5)]++;
	if (job_status[sus_person]==5) {
	  num_infect_HCW++;
	}

	/* Determine if following interventions only for interventions that effect individuals.*/
	if ( interventions > 0 && COV_rand() < complyI[3] ) {
	  if ( interventions == 2 && intervene[sus_person] == 4 ) {
	    intervene[sus_person]=5;
	  } else {
	    intervene[sus_person]=3;
	  }
	  /* Intervention 2 is household quarantine with current recommendations. Applicable for whole household.  */
	  if ( interventions == 2 && t>tauI_onset ) {
	    int i1, hh = HH[sus_person];
	    for (i1=0; i1 < arrlen(per_HH_members[hh]); i1++) {
	      if ( COV_rand()<complyI[4] ) {
		intervene[per_HH_members[hh][i1]]=4;
	      }
	    }
	  }
	}
      }
    }

    // After 5 days of symptoms, people are randomly put into hospital based on age distributions. //
    hosp_entry(t, num_infectious,  infectious,  age,  icu_pop,  hosp_pop,  symptomatic, tau, workplace_tmp, hosp_num, county, dt) ;

    // Release from hospital //
    hosp_release(t, num_hosp,  hosp_list,  tau,  recovered, hosp_pop, &num_recovered, &recovered_hosp, &recovered_icu, dt, num_recovered_county, num_recovered_age, age, county, num_recovered_hosp_county, num_recovered_hosp_age, num_recovered_icu_county, num_recovered_icu_age, &num_recovered, &recovered_hosp, &recovered_icu, job_status ) ;


    // 15 days after onset of symptoms, people randomly die based on age based percentages.  ICU patients are removed from ICU into regular hospital. //

    num_dead = death(t, num_infectious,  infectious,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt, num_dead_county, num_dead_age, county, &num_dead_HCW, job_status) ;
    num_dead = death(t, num_icu,  icu_list,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt, num_dead_county, num_dead_age, county, &num_dead_HCW, job_status) ;
    num_dead = death(t, num_hosp,  hosp_list,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt, num_dead_county, num_dead_age, county, &num_dead_HCW, job_status) ;


    // Recovered after 11 days (6 days of symptoms) if not in hospital/ICU. //
    for (i=0; i<num_infectious; i++) {
      int infec_person;
      infec_person=infectious[i];
      if ((tau[infec_person]<(t-11)+dt) && (tau[infec_person]>=t-11) && (hosp_pop[infec_person]==0) && (icu_pop[infec_person]==0)) {
	recovered[infec_person]=1;
	num_recovered++;
	if (job_status[infec_person]==5) {
	  num_recovered_HCW++;
	}
	num_recovered_county[county[infec_person]]++;
	num_recovered_age[(int)floor(age[infec_person]/5)]++;

      }
    }

    ret = clock_gettime(CLOCK_MONOTONIC, &t2);
    step_time = ((double)t2.tv_sec + (double)t2.tv_nsec/nsdiv) - ((double)t1.tv_sec + (double)t1.tv_nsec/nsdiv);
    printf("time %5.2f\n", step_time);
    fflush(stdout);

    fprintf(output_file, "Walltime/timestep %6.2f Time %6.2f num_infected %i num_infectious %i num_in_hosp %i num_in_icu %i num_dead %i recovered_tot %i recovered_from_hosp %i recovered_from_icu %i contact_work %i contact_school %i contact_home %i contact_community %i \n", step_time, t, num_infect, num_infectious, num_hosp, num_icu, num_dead, num_recovered, recovered_hosp, recovered_icu, num_contact_work, num_contact_school, num_contact_house, num_contact_commun);
    fflush(output_file);

    fprintf(output_HCW, "Walltime/timestep %6.2f Time %6.2f num_infected %i num_infectious %i num_in_hosp %i num_in_icu %i num_dead %i recovered_tot %i recovered_from_hosp %i recovered_from_icu %i contact_work %i contact_school %i contact_home %i contact_community %i total_healthcare %i \n", step_time, t, num_infect_HCW, num_infectious_HCW, num_hosp_HCW, num_icu_HCW, num_dead_HCW, num_recovered_HCW, num_recovered_hosp_HCW, num_recovered_icu_HCW, num_contact_work_HCW, num_contact_school_HCW, num_contact_house_HCW, num_contact_commun_HCW, num_HCW);
    fflush(output_HCW);

    for (file_count=0; file_count<num_counties; file_count++) {
      fprintf(county_files[file_count], "Walltime/timestep %6.2f Time %6.2f num_infected %i num_infectious %i num_in_hosp %i num_in_icu %i num_dead %i recovered_tot %i recovered_from_hosp %i recovered_from_icu %i contact_work %i contact_school %i contact_home %i contact_community %i total_individuals %i \n", step_time, t, num_infect_county[file_count], num_infectious_county[file_count], num_hosp_county[file_count], num_icu_county[file_count], num_dead_county[file_count], num_recovered_county[file_count], num_recovered_hosp_county[file_count], num_recovered_icu_county[file_count], num_contact_work_county[file_count], num_contact_school_county[file_count], num_contact_house_county[file_count], num_contact_commun_county[file_count], county_size[file_count]);
      fflush(county_files[file_count]);
    }

    for (file_count=0; file_count<90; file_count+=5) {
      fprintf(age_files[file_count/5], "Walltime/timestep %6.2f Timestep %6.2f num_infected %i num_infectious %i num_in_hosp %i num_in_icu %i num_dead %i recovered_tot %i recovered_from_hosp %i recovered_from_icu %i contact_work %i contact_school %i contact_home %i contact_community %i total_individuals %i \n", step_time, t, num_infect_age[file_count/5], num_infectious_age[file_count/5], num_hosp_age[file_count/5], num_icu_age[file_count/5], num_dead_age[file_count/5], num_recovered_age[file_count/5], num_recovered_hosp_age[file_count/5], num_recovered_icu_age[file_count/5], num_contact_work_age[file_count/5], num_contact_school_age[file_count/5], num_contact_house_age[file_count/5], num_contact_commun_age[file_count/5], age_distrib[file_count/5]);
      fflush(age_files[file_count/5]);
    }

    t=t+dt;
  }
  ret = clock_gettime(CLOCK_MONOTONIC, &T2);
  step_time = ((double)T2.tv_sec + (double)T2.tv_nsec/nsdiv) - ((double)T1.tv_sec + (double)T1.tv_nsec/nsdiv);
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

  for (i=0; i<num_counties; i++) {
    fclose(county_files[i]);
  }

  for (i=0; i<85; i+=5) {
    fclose(age_files[(i/5)]);
  }
  fclose(output_file);
  fclose(lat_lon_out);
  fclose(output_HCW);

}
