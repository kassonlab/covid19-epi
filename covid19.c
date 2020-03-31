#include <stdio.h>
 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "common.h"

int typical_max_HH_sz = 7;

int full_fd = 0, full_kappa = 0;

float betac_scale = 8.4, betah_scale = 2.0, betaw_scale = 1.0;

//Taken from geodatasource.com //
/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::                                                                         :*/
/*::  This routine calculates the distance between two points (given the     :*/
/*::  latitude/longitude of those points). It is being used to calculate     :*/
/*::  the distance between two locations using GeoDataSource(TM) products.   :*/
/*::                                                                         :*/
/*::  Definitions:                                                           :*/
/*::    South latitudes are negative, east longitudes are positive           :*/
/*::                                                                         :*/
/*::  Passed to function:                                                    :*/
/*::    lat1, lon1 = Latitude and Longitude of point 1 (in decimal degrees)  :*/
/*::    lat2, lon2 = Latitude and Longitude of point 2 (in decimal degrees)  :*/
/*::    unit = the unit you desire for results                               :*/
/*::           where: 'M' is statute miles (default)                         :*/
/*::                  'K' is kilometers                                      :*/
/*::                  'N' is nautical miles                                  :*/
/*::  Worldwide cities and other features databases with latitude longitude  :*/
/*::  are available at https://www.geodatasource.com                         :*/
/*::                                                                         :*/
/*::  For enquiries, please contact sales@geodatasource.com                  :*/
/*::                                                                         :*/
/*::  Official Web site: https://www.geodatasource.com                       :*/
/*::                                                                         :*/
/*::           GeoDataSource.com (C) All Rights Reserved 2018                :*/
/*::                                                                         :*/
/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

#define pi 3.14159265358979323846

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  Function prototypes                                           :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double deg2rad(double);
double rad2deg(double);

double distance(double lat1, double lon1, double lat2, double lon2, char unit) {
  double theta, dist;
  if ((lat1 == lat2) && (lon1 == lon2)) {
    return 0;
  }
  else {
    theta = lon1 - lon2;
    dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
    dist = acos(dist);
    dist = rad2deg(dist);
    dist = dist * 60 * 1.1515;
    switch(unit) {
      case 'M':
        break;
      case 'K':
        dist = dist * 1.609344;
        break;
      case 'N':
        dist = dist * 0.8684;
        break;
    }
    return (dist);
  }
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts decimal degrees to radians             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double deg2rad(double deg) {
  return (deg * pi / 180);
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts radians to decimal degrees             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double rad2deg(double rad) {
  return (rad * 180 / pi);
}
/// End of code from GEODATASOURCE


void age_dist (float * age, int population, FILE* stats, int * age_distrib) {

	int i; // Counter
        int ret, age_sz;

	int age_start[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
	double age_dist[] = {0.10721, 0.11553, 0.12439, 0.13475, 0.12578, 0.12704, 0.10782, 0.09952, 0.05796};
        double *age_dst;

        age_sz = (int)(sizeof(age_dist)/sizeof(age_dist[0]));
        ret = generate_inc_distr_vec(&age_dst, age_dist, age_sz, "age");
        if (ret) {
            fprintf(stderr, "Bailing out on age initialization\n");
            exit(1);
        }

	for (i=0; i<population; i++) {
                int indx;
                indx = rand_distr_indx(age_dst, 0, age_sz-1);
		age[i] = age_start[indx] + COV_rand()*10;
		age_distrib[(int)floor(age[i]/5)]++;
	}
        free(age_dst);


// Uncomment to test age distribution.  Tested JMG 2020-03-20.
	int age_dist_test[9]={0};
	for (i=0; i<population; i++) {
		age_dist_test[(int)floor(age[i]/10)]++;
	}	
	
        fprintf(stats, "Testing age distribution\n");
	for (i=0; i<9; i++) {
		fprintf(stats, "age %i percent %f num %f \n", i, (float)age_dist_test[i]/population, age_dist[i]) ; 
	}
	fprintf(stats, "\n\n");
        fflush(stats);
       // exit(0);
//
}



/* Puts households in certain locations based on population density distribution.  Only really correct for full population but works for smaller populations.  Biases towards smaller households for smaller populations.  Each household is also fed into a locality based on the shortest distance to the center of that locality for purposes of school and workplace choice. */
void household_lat_long(int num_households, int * HH, float * lat_city, float * long_city, int num_cities, int * city, int * county, int * city_county, int * city_size, int * county_size, int population, float * age, int * per_HH_size, char ** county_names, float * county_pop, float tot_pop_actual, FILE* stats, int **county_p, int *num_locale, float *lat_locale, float *lon_locale, float *pop_density_init_num, int **locale_to_HH, int *locale_to_HH_n, int *locale_HH) {

	/* Initialize helpers for population density */
	int num_km=*num_locale;
	// FOR TEST TEXT 
	//int num_km=63;
	float tot_pop_density=0;

	/* Initialize helpers for household distribution */
	float lat_HH[num_households];	
	float lon_HH[num_households];	
	int city_HH[num_households];	
	int county_HH[num_households];	
	int county_num;	
	int city_num;	
	
	float dist1;
	float min_dist=1000000;
	int placement;
        int i, j;

	int tmp_county_count[21]={0};	
	int tmp_county_density[21]={0};	

	/* find how many locales will have households */
	if (num_households<num_km) {
		*num_locale=num_households;
	} else {
		*num_locale=num_km;
	}

	/* Fill up households. */
	int HH_count=0 ; // Counter
	int HH_person=0 ; // Counter
	int max_HH_size=0;

	/* Initially set all households to -1 */
	for (i=0; i < population; i++) {
		HH[i]=-1;
	}

	/* Save list of households in each locale. */
	int ** county_list;
        county_list = (int**)calloc(*num_locale,sizeof(int*));
	for (i=0;i<*num_locale;i++) county_list[i] = (int*)calloc(num_households,sizeof(int));
	int * locale_count;
	locale_count = (int*)calloc(*num_locale,sizeof(int));
	int * locale_HH_count;
	locale_HH_count = (int*)calloc(*num_locale,sizeof(int));
	

        city_num = -1;
	while (( HH_count < num_households ) && (HH_count < num_km)) {
                float tmp_lat, tmp_lon;
		min_dist=1000000;

                tmp_lat = lat_locale[HH_count];
                tmp_lon = lon_locale[HH_count];
		tot_pop_density+=pop_density_init_num[HH_count];

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
		lat_HH[HH_count]=tmp_lat;
		lon_HH[HH_count]=tmp_lon;
		city_HH[HH_count]=city_num;
		county_HH[HH_count]=county_num;
		county_list[HH_count][locale_count[HH_count]]=HH_count;
                locale_HH[HH_count] = HH_count;
                locale_to_HH[HH_count][locale_to_HH_n[HH_count]++] = HH_count;
		locale_HH_count[HH_count]+=1;

                /* Allocate an adult to the household */
		while ( age[HH_person]<20 ) {
			HH_person++;
		}
		HH[HH_person]=HH_count;
		city[HH_person]=city_HH[HH_count];
		county[HH_person]=county_HH[HH_count];
                county_p[county[HH_person]][county_size[county[HH_person]]++] = HH_person;
		city_size[city[HH_person]]++;
		per_HH_size[HH[HH_person]]++;
		locale_count[HH_count]+=1;
		if (per_HH_size[HH[HH_person]]>max_HH_size) {
			max_HH_size=per_HH_size[HH[HH_person]];
		}
		HH_person++;
		HH_count++;
	}

/* Uncomment to test population denity WRT county. 
	// Test population based on method vs actual population density on county level.  Lines up fairly well. //
	for (i=0;i<21; i++) {
		printf("counties %i %s number of locales %i calculated density %f actual density %f\n", i, county_names[i], tmp_county_count[i], tmp_county_density[i]/(float)tot_pop_density, county_pop[i]/tot_pop_actual);
		fflush(stdout);	
	}
*/

	int list=1; //This keeps track of how many times we can been through the locale list.	
	placement=0; //Keeps track of household placement after first loop of locales.
	/* First add adult head of household to the rest of the households. */
	while ( HH_count < num_households ) {
		while ( age[HH_person]<20 ) {
			HH_person++;
		}
		/* Place another household in locales with population density greater than 2*(households in locale). */
		if (pop_density_init_num[placement]<=2*list) {
				list++;	
				placement=0; // Start over at top of list.  Remember list is sorted by biggest to smallest locales.
		}	
		lat_HH[HH_count]=lat_HH[placement];	
		lon_HH[HH_count]=lon_HH[placement];
		city_HH[HH_count]=city_HH[placement];
		county_HH[HH_count]=county_HH[placement];
		county_list[placement][locale_count[placement]]=HH_count;
                locale_HH[HH_count] = placement;
                locale_to_HH[placement][locale_to_HH_n[placement]++] = HH_count;
		locale_HH_count[placement]+=1;

		/* Set up head of household. */
		HH[HH_person]=HH_count;
		city[HH_person]=city_HH[HH_count];	
		county[HH_person]=county_HH[HH_count];	
                county_p[county[HH_person]][county_size[county[HH_person]]++] = HH_person;
		city_size[city[HH_person]]++;
		locale_count[placement]+=1;
		per_HH_size[HH[HH_person]]++;
		if (per_HH_size[HH[HH_person]]>max_HH_size) {
			max_HH_size=per_HH_size[HH[HH_person]];
		}
		placement+=1;
		HH_person++;
		HH_count++;
	}

        placement = 0;
	/* Distribute remaining people randomly.  This could be changed to a distribution to more realistically reflect household size in the future. */
	for ( HH_person=0; HH_person<population ; HH_person++) {
                int tmp_HH;
		if (HH[HH_person]==-1) {
			
			/* Place people in random households within locale until locale is full. */
			if (placement>=*num_locale) {
				placement=0;
			} else if (locale_count[placement]>=pop_density_init_num[placement]) {
					placement=0; // Start over at top of list.  Remember list is sorted by biggest to smallest locales.
			}	

			/* Pick a random household in the locale. */
			tmp_HH=county_list[placement][(int)(COV_rand() * locale_HH_count[placement])];
                        while (per_HH_size[tmp_HH]+1 > typical_max_HH_sz) {
                            tmp_HH=county_list[placement][(int)(COV_rand() * locale_HH_count[placement])];
                        }
			HH[HH_person]=tmp_HH;
		//	printf("HH %i %i %i %f \n", HH[HH_person], placement, locale_HH_count[placement], HH_person);
			city[HH_person]=city_HH[HH[HH_person]];	
			county[HH_person]=county_HH[HH[HH_person]];	
                        county_p[county[HH_person]][county_size[county[HH_person]]++] = HH_person;
			city_size[city[HH_person]]++;
			locale_count[placement]+=1;
                        per_HH_size[HH[HH_person]]++;
                        if (per_HH_size[HH[HH_person]]>max_HH_size) {
                                max_HH_size=per_HH_size[HH[HH_person]];
                        }
			placement+=1;
		}
	}

        printf("max_HH_size = %d\n", max_HH_size);
        fflush(stdout);


// Uncomment to test household distribution.  Tested JMG 2020-03-22.
	int *HH_dist_test;
        HH_dist_test = (int *)calloc(max_HH_size+1, sizeof(int));
	for (i=0; i<num_households; i++) {
		HH_dist_test[per_HH_size[i]]++;
	}	
	
	fprintf(stats, "Household distributions \n");
	for (i=0; i<=max_HH_size; i++) {
		fprintf(stats, "household_size %i percent_households %f num_households %i total_households %i \n", i, HH_dist_test[i]/(float)num_households, HH_dist_test[i], num_households) ; 
		fflush(stats);
	}
        free(HH_dist_test);
	
	fprintf(stats, "\n\n County distribution \n");
	for (i=0; i<21; i++) {
		fprintf(stats, "%s county %i population %i percent %f actual %f \n", county_names[i],i, county_size[i], county_size[i]/(float)population, county_pop[i]/tot_pop_actual)  ;
		fflush(stats);
	} 
	fprintf(stats, "\n\n");
        fflush(stats);

        for (i = 0; i < *num_locale; i++) {
            free(county_list[i]);
        }
	free(county_list);
	free(locale_count);
	free(locale_HH_count);
}

void city_lat_long(int *num_cities, float * lat_city, float * long_city, char ** cities, int * county, char ** county_names, int num_county) {

	/* Get longitude and latitude of cities in Sweden from CSV */
	FILE* fp = fopen("cities_all.csv", "r");  // Not sure about the validity of this file.  Could use a better source.

	/* Counters for parsing file. */
	char tmp[199]; // municipality name
	char tmp1[199]; // locality name
	char tmp2[199]; // county name
	float tmp_lat; // tmp latitude
	float tmp_lon; // tmp longitude
        int i, j;

	for (i=0; i < 2000; i++) {
		int got = fscanf(fp, "%[^,\n]%*c%[^,\n]%*c%[^,\n]%*c%f%*c%f\n", tmp1, tmp, tmp2, &tmp_lat, &tmp_lon);
  		if (got != 5) break; // wrong number of tokens - maybe end of file
			
		/* Save city name, longitude, latitude, and increase number of cities. Looking at unique municipalities.*/
		cities[*num_cities]=tmp1;
		lat_city[*num_cities]=tmp_lat;
		long_city[*num_cities]=tmp_lon;
		*num_cities=*num_cities+1;
	
		/* Save county that the city is in. */
		for (j=0; j<num_county; j++) {
			if ((int)(strcmp(county_names[j], tmp2))==0) {
				county[i]=j;
				break;
			} else if (j==num_county-1) {
				printf("no_county_match %s %s %i %i\n", county_names[j], tmp2, j, strcmp(county_names[j], tmp2));
			}	
		}
	}

	fclose(fp);
}


void job_dist(int * job_status, int ** job_status_city, float * age, int * county, int * city, int population, int num_cities, int num_counties, int * county_size, FILE * stats) {

	int i, j; 
	

	/* All currently randomly placed based on county.  Would like to do per town but each town needs inhabitants. See commented section below for per city distribution of schools */ 
	for (i=0; i < population; i++) {
		if (age[i]<1 || age[i]>75) {
			job_status[i]=0;
			job_status_city[0][county[i]]++;
		} else if (age[i]>=1 && age[i]<3) {
			if (COV_rand() < 0.7800) {
				job_status[i]=1;
				job_status_city[1][county[i]]++;
			} else {
				job_status[i]=0;
				job_status_city[0][county[i]]++;
			}	
		} else if (age[i]>=3 && age[i]<6) {
			if (COV_rand() < 0.9500) {
				job_status[i]=1;
				job_status_city[1][county[i]]++;
			} else {
				job_status[i]=0;
				job_status_city[0][county[i]]++;
			}	
		} else if (age[i]>=6 && age[i]<15) {
			job_status[i]=2;
			job_status_city[2][county[i]]++;
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


// Uncomment to test job distribution.  Tested JMG 2020-03-20. 
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
			fprintf(stats, "job_status %i county %i percent_of_jobs_total %f num_jobs_in_county %i unemployed %i percent_unemployed  %f \n", i, j, job_dist_test[i]/(float)population, city_dist_test[i][j], unemployed[j], (float)unemployed[j]/working_age[j]) ;
			 
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

//	

}

void workplace_dist(int * workplace, int * job_status, int ** job_status_county, int * city, int num_cities, int * county, int num_counties, int population, int * max_num_WP , int * hosp_num, int* class, FILE * stats) {

	int pp_class = 19; //Assumption of 15 children per class.
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
				num_workplaces[j][i]=ceil(job_status_county[j][i]/(float)pp_preschool);
				num_workplaces2[j]+=ceil(job_status_county[j][i]/(float)pp_preschool);
				if (num_workplaces2[j]>*max_num_WP) {
					*max_num_WP=num_workplaces2[j];
				}
			}
		}
		for (j=2; j<3; j++) {
			/* First only schools then workplaces */
			if (job_status_county[j][i]>0) {
				num_workplaces[j][i]=ceil(job_status_county[j][i]/(float)pp_school);
				num_workplaces2[j]+=ceil(job_status_county[j][i]/(float)pp_school);
				if (num_workplaces2[j]>*max_num_WP) {
					*max_num_WP=num_workplaces2[j];
				}
			}
		}
	}

	// Broken down in case we want to do schools by municipality.
	for (i=0; i < num_counties; i++) {
		for (j=3; j<5; j++) {
                    if (job_status_county[j][i]>0) {
                            num_workplaces[j][i]=ceil(job_status_county[j][i]/(float)pp_work);
                            num_workplaces2[j]+=ceil(job_status_county[j][i]/(float)pp_work);
                            if (num_workplaces2[j]>*max_num_WP) {
                                    *max_num_WP=num_workplaces2[j];
                            }
                    }
                }
		// Need to use floor+1 equation because each county should have a hospital even if no one works there. */
		num_workplaces[5][i]=floor(job_status_county[5][i]/(float)pp_hospital)+1;
		num_workplaces2[5]+=floor(job_status_county[5][i]/(float)pp_hospital)+1;
		if (num_workplaces2[5]>*max_num_WP) {
			*max_num_WP=num_workplaces2[5];
		}
		hosp_num[i]=num_workplaces[5][i];
	}

	int max_num_classes=0;
	for (i=0; i < population; i++) {
		//Try to minimize necessary memory by making job_numbers independent with job_status. //
		int prior_workplaces=0;
		if (job_status[i] < 3 && num_workplaces[job_status[i]][city[i]] > 0) {
                        for (j=0; j < city[i]; j++) {
                                prior_workplaces+=num_workplaces[job_status[i]][j];
                        }
			workplace[i]=(int)(COV_rand() * (num_workplaces[job_status[i]][city[i]]))+prior_workplaces;
			class[i]=(int)(COV_rand() * (ceil(num_workplaces[job_status[i]][city[i]]/(float)pp_class)));

			if (ceil(num_workplaces[job_status[i]][city[i]]/(float)pp_class)>max_num_classes) {
				max_num_classes=ceil(num_workplaces[job_status[i]][city[i]]/(float)pp_class);
			}
	
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


/* Uncomment to workplace job distribution.  Tested JMG 2020-03-20. 
	int k;
	int max=*max_num_WP;
	int * work_dist_test;
	work_dist_test = (int*)calloc((max)*num_counties,sizeof(int));
	int * class_dist_test;
	class_dist_test = (int*)calloc((max)*num_counties*max_num_classes,sizeof(int));
	for (i=0; i<population; i++) {
		work_dist_test[workplace[i]][(int)county[i]]++;
		class_dist_test[workplace[i]][(int)county[i]][class[i]]++;
	}	

	fprintf(stats, "Workplace Distribution \n");
	for (j=0; j<num_counties; j++) {
		for (i=0; i<*max_num_WP; i++) {
			if (work_dist_test[i][j]>0) {
				fprintf(stats, "county %i workplace %i size %i \n", j, i, work_dist_test[i][j]) ;
			}
			for (k=0; k<max_num_classes; k++) {
				if (class_dist_test[i][j][k]>0) {
					fprintf*stats, "county %i school %i class %i size %i \n", j, i, k, class_dist_test[i][j][k]) ;
				}	
			}
		}
	}
        fflush(stats);

	free(work_dist_test);
	free(class_dist_test);
*/	

}

float calc_kappa(float t, float tau, int symptomatic, float dt, float * kappa_vals) {

	float kappa;
	float t1;
	int t2;
	//###Determine kappa for infected person.  This is the infectiousness of the person based on time since infection started.  Latency period is 4.6 days.  Infection starts at 5.1 days and lasts for 6 days.  Sympotmatic people are twice as likely to infect others as asymptomatic.
	// Kappa is a log normal function with mean of -0.72 and standard deviation of 1.8.  From Ferguson Nature 2005
	if (t-tau <= 4.6) {
		kappa=0.;
	} else if (t-tau>11.1) {
		kappa=0.; //# Recovered or dead
	} else {
		/* First 2 lines calculates kappa on the fly, second two get precalculated kappa from array. */
                if (full_kappa) {
                    t1=(log(t-tau-4.6)+0.72)/1.8;
                    kappa=exp(-0.5*pow(t1,2))/((t-tau-4.6)*1.8*sqrt(2*pi));
                } else {
                    t2=(t-tau)/dt;
                    kappa=kappa_vals[t2];
                }
	}
	if (symptomatic==0) {
		kappa=kappa*0.5;
	}
	return(kappa);
}

void initialize_infections(int * initial_infections, float * tau, int * infected, int * severe, int * symptomatic, int * county, int * num_infect, int num_counties, float symptomatic_per, int population, float dt, float t, float * lat_locale, float * lon_locale, int * num_infect_county, int * num_infect_age, float * age, int **county_p, int *county_size, int *locale_HH, int *HH) {

	int person_infected=0, county_person_inf = 0;
	int tmp_infect=0;
	float min_diff=1000;
	float diff_lat_lon=10;
        float *tmp_lat, *tmp_lon;
        int i, j;

        tmp_lat = (float*)calloc(population,sizeof(float));
        tmp_lon = (float*)calloc(population,sizeof(float));
	
	for (i=0; i < num_counties; i++) {
		tmp_infect=0;
		while ((tmp_infect<initial_infections[i])) {
			int tmp_j=0;
			// Test up to population to see if we can find someone who fits a perviously determined cluster.  If not, leave this loop and pick a random person.
                        county_person_inf = (int)(COV_rand() * county_size[i]);
                        person_infected=county_p[i][county_person_inf];
			while (((county[person_infected]!=i) || (infected[person_infected]!=0) || diff_lat_lon>1.0) && tmp_j<population) {
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
				if (t<-10 || num_infect_county[i]==0 ) {
					diff_lat_lon=0;
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
			if (diff_lat_lon>1) {
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
	//		tau[person_infected]=-COV_rand() * 5;
			tau[person_infected]=t;
			tmp_lat[*num_infect]=lat_locale[locale_HH[HH[person_infected]]];
			tmp_lon[*num_infect]=lon_locale[locale_HH[HH[person_infected]]];
			*num_infect=*num_infect+1;
			num_infect_county[i]++;
			num_infect_age[(int)(floor(age[i]/5))]++;
			tmp_infect++;
		}
			
	}
        free(tmp_lat);
        free(tmp_lon);
}


void segment_population(int* num_sus, int* num_infectious, int* num_hosp, int* num_icu, int* infected, int* infectious, int* sus_list, int* hosp_list, int * hosp_pop, int * icu_pop, int* icu_list, float* tau, int population, float t, int * num_sus_county, int * num_infectious_county, int * num_infected_county, int * num_hosp_county, int * num_icu_county, int * num_sus_age, int * num_infectious_age, int * num_infected_age, int * num_hosp_age, int * num_icu_age, float * age, int * county, int print_loc, float * lat_locale, float * lon_locale, FILE * lat_lon_out, int *locale_HH, int *HH) {

        int i;

	*num_sus=0;
	*num_infectious=0;
	*num_hosp=0;
	*num_icu=0;
	for (i=0; i<population; i++) {
		if (infected[i]==0) {
			sus_list[*num_sus]=i;
			*num_sus=*num_sus+1;
			num_sus_county[county[i]]++;
			num_sus_age[(int)floor(age[i]/5)]++;
		} else if ((tau[i]<t-4.6) && (tau[i]>t-11.1)) {
			infectious[*num_infectious]=i;
			*num_infectious=*num_infectious+1;
			num_infectious_county[county[i]]++;
			num_infectious_age[(int)floor(age[i]/5)]++;
			if (print_loc==1) {
				fprintf(lat_lon_out, "Time %f Person %i lat %f lon %f \n", t, i, lat_locale[locale_HH[HH[i]]], lon_locale[locale_HH[HH[i]]]);
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
		} 
		if (hosp_pop[i]>0) {
			hosp_list[*num_hosp]=i;
			*num_hosp=*num_hosp+1;
			num_hosp_county[county[i]]++;
			num_hosp_age[(int)floor(age[i]/5)]++;
		}
		
	}
	fprintf(lat_lon_out, "\n\n");	

/* Uncomment for information */
//	printf("Time %i infections i %i num_sus %i infected[i] %i num_infectious %i num_hosp %i num_icu %i \n", t, i, *num_sus, infected[i], *num_infectious, *num_hosp, *num_icu);

}

float calc_household_infect(float kappa, float omega, int HH_size, float alpha, int severe) {

	float betah=0.55; // Scaled from betah=0.4 in influenza pandemic with R0=1.6, COVID-19 R0=2.2 (Ferguson 2020)

	return(betah_scale*betah*kappa*(1+(float)severe*(omega-1))/(pow((float)HH_size,alpha))); 

}

float calc_workplace_infect(int job_status, float kappa, float omega, int workplace_size, int severe, float * Iw) {

	float betap[]={0.0, 1.1, 1.1, 1.1, 0.55, 0.1375} ; // Spread in all types of schools (preschool to college) is twice that of workplace 
	float psi[]={0.0, 0.1, 0.2, 0.25, 0.5, 0.5} ; // Accounts for absenteeism based on severe illness. Ferguson Nature 2006

	return(betaw_scale*Iw[job_status]*betap[job_status]*kappa*(1+(float)severe*(omega*psi[job_status]-1))/((float)workplace_size));
}

float calc_community_infect(int age_group, float kappa, float omega, int severe, float d, double * fd_vals) {

	/* need to work on this.  Perhaps we take a random distance for each two people based on population density, number of people in county, county area, etc. */
	float zeta[]={0.1, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.75, 0.50, 0.25, 0.25, 0.25} ; //   # Travel related parameter for community transmission. Ferguson Nature 2006
	float fd, fd1;
	float betac=0.103 ; // Scaled from betac=0.075 in influenza pandemic with R0=1.6, COVID-19 R0=2.2 (Ferguson 2020)

        if (full_fd) {
            fd=1/(1+pow((d/4), 3)); //kernel density function as parameterized for GB.
        } else {
            fd=fd_vals[(int)(d*10)];
        }
	return(betac_scale*zeta[age_group]*betac*kappa*fd*(1+severe*(omega-1)));
}


void hosp_entry(float t, int num_infectious, int * infectious, float * age, int * icu_pop, int * hosp_pop, int * symptomatic, float * tau, int * workplace_tmp, int * workplace_num, int * county, float dt) {

	float hosp[]={0.001, 0.003, 0.012, 0.032, 0.049, 0.102, 0.166, 0.243, 0.273}; //# From Ferguson 2020
	float icu[]={0.05, 0.05, 0.05, 0.05, 0.063, 0.122, 0.274, 0.432, 0.709} ; //# percent of hospitalized that need icu from Ferguson 2020.
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


void hosp_release(float t, int num_hosp, int * hosp_list, float * tau, int * recovered, int * hosp_pop, int * num_recovered, int * recovered_hosp, int * recovered_icu, float dt, int * num_recovered_county, int * num_recovered_age, float * age, int * county, int * num_recovered_hosp_county, int * num_recovered_hosp_age, int * num_recovered_icu_county, int * num_recovered_icu_age) {

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
		} else if ((tau[infec_person]>=t-25) && (tau[infec_person]<t-25+dt) && (hosp_pop[infec_person]==2)) {
			hosp_pop[infec_person]=0;
			recovered[infec_person]=1;
			*num_recovered=*num_recovered+1;
			*recovered_icu=*recovered_icu+1;
			num_recovered_county[county[infec_person]]++;
			num_recovered_age[(int)floor(age[infec_person]/5)]++;
			num_recovered_icu_county[county[infec_person]]++;
			num_recovered_icu_age[(int)floor(age[infec_person]/5)]++;
		}
	}
}

int death(float t, int num_infectious, int * infectious, float * tau, int * dead, int * icu_pop, int * hosp_pop, int * symptomatic, int num_dead, float * age, float dt, int * num_dead_county, int * num_dead_age, int * county) {

	float fatal_in_icu=0.5;
	float fatal_symptomatic[]={0.00002, 0.00006, 0.0003, 0.0008, 0.0015, 0.006, 0.022, 0.051, 0.093};
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
				}
			}
		}
	}
	return(num_dead);
}


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
	float percent_infect=0.001 ; // Default 1% of population initially infected.
	float dt=1.00; // Time step.
	int interventions=0; // Value of interventions.
	float tauI_onset=1; //time after start of simulation that interventions for whole community take place.
	int print_lat_lon=0; // Choose whether to print latitude and longitude data of infected individuals for each time step. 
        int ret;
        int i, j;

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
    		else if (!strcmp(argv[i],"-full_fd")) full_fd = 1;
    		else if (!strcmp(argv[i],"-full_kappa")) full_kappa = 1;
    		else if (!strcmp(argv[i],"-use_fixed_seed")) use_fixed_seed = 1;
  	}

	/* print out population statistics to file */
	FILE *stats = fopen("stats.log", "w");
	/*Get current time and stats */
	time_t date;
	time(&date);

	fprintf(stats, "Date: %s \n", ctime(&date));
	fprintf(stats, "Population: %i \nSimulationTime: %i \ndt: %f \nInterventions: %i \nTau_onset: %f \n\n\n ", population, tot_time, dt, interventions, tauI_onset);	
 
        ret = clock_gettime(CLOCK_MONOTONIC, &T1);

	float HH_size = 2.2 ; // Average household size from OCED 
	int num_households = (int)(population/HH_size); // Number of households in population
	int * HH;  // Households 
	HH = (int*)calloc(population,sizeof(int));
	int * per_HH_size; // Size of each household.  Need for infectiousness calculations.
	per_HH_size = (int*)calloc(num_households,sizeof(int));

	/* Population information */
	/* Specific to Sweden.  Averaging population based on total population of Sweden regardless of population size. */
	float tot_pop=10327589.; //For use of averaging counties only.  NOTE: land scan tot is 10,098,554, cannot use more inhabitants than that. 
	char * county_name[]={"Stockholm", "Uppsala", "Sodermanland", "Ostergotland", "Jonkoping", "Kronoberg", "Kalmar", "Gotland", "Blekinge", "Skane", "Halland", "Vastra Gotaland", "Varmland", "Orebro", "Vastmanland", "Dalarna", "Gavleborg", "Vasternorrland", "Jamtland", "Vasterbotten", "Norrbotten"}; 
	float pop_county[]={2377081, 383713, 297540, 465495, 363599, 201469, 245446, 59686, 159606, 1377827, 333848, 1725881, 282414, 304805, 275845, 287966, 287382, 245347, 130810, 271736, 250093};
	int num_counties=(sizeof(pop_county)/sizeof(pop_county[0]));
	int county_int[num_counties]; // Integer values for random probability distribution. 
	memset(county_int, 0, num_counties*sizeof(int));
	float * lat_city; //latitude of city i.
	double * fd_tot; //additive kernel density function of person with respect to all other individuals
	fd_tot = (double*)calloc(population,sizeof(double));
	lat_city = (float*)calloc(2000,sizeof(float));
	float * long_city; // longitude of city i.
	long_city = (float*)calloc(2000,sizeof(float));
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
	float * age;  // Age of population
	age = (float*)calloc(population,sizeof(float));
	int * age_distrib;
	/* Distribution of ages in the generated population, people probably won't be older ehan 150 */
	age_distrib = (int*)calloc(150,sizeof(int));
	int * county;  // County of each inhabitant
	county = (int*)calloc(population,sizeof(int));
        int **county_p; /* List of persons per county */
        county_p = (int **)malloc(num_counties * sizeof(int *));
        for(i=0; i < num_counties; i++) {
            county_p[i] = (int *) malloc(population * sizeof(int)); /* AS: can be reduced to max number of persons in the most populated county */
        }
	int * job_status; // Type of job each person holds: 0-4 for no job, preschool, elementary school, highschool/college, and job, respectively. 	
	job_status = (int*)calloc(population,sizeof(int));
	int * workplace; // Workplace of each person.
	workplace = (int*)calloc(population,sizeof(int));
	int * class; // Classroom for students 	
	class = (int*)calloc(population,sizeof(int));

	int max_num_WP=0; // max workplaces per job_status for allocating array.

	/* Parameters for infections */
//	int num_infections=(int)population*percent_infect; // Default is 10% of population has illness.
	float symptomatic_per=0.67; // percent of people who are symptomatic.
	int * infected; // 1 if person i has been infected, 0 otherwise
	infected = (int*)calloc(population,sizeof(int));
	int * severe; // 1 if person i has severe infection, 0 otherwise
	severe = (int*)calloc(population,sizeof(int));
	float * tau; // day that person i became infected
	tau = (float*)calloc(population,sizeof(float));
	int * symptomatic; // 1 if person i is symptomatic, 0 otherwise
	symptomatic = (int*)calloc(population,sizeof(int));

	/* Infection parameters.  Mostly counters for looping. */
	int num_sus=0; //Number susceptible.
	int num_sus_age[18]={0};
	int num_infect=0; //Number infected.
	int *num_sus_county;
        num_sus_county = (int *)calloc(num_counties, sizeof(int));
	int * num_infect_county;
	num_infect_county = (int*)calloc(num_counties, sizeof(int));
	int num_infect_age[18]={0};
	int num_icu=0; //Number in ICU.
	int *num_icu_county;
        num_icu_county = (int *)calloc(num_counties, sizeof(int));
	int num_icu_age[18]={0};
	int num_hosp=0; //Number in hospital. 
	int *num_hosp_county;
        num_hosp_county = (int *)calloc(num_counties, sizeof(int));
	int num_hosp_age[18]={0};
	int num_dead=0; //Number of deaths. 
	int *num_dead_county;
        num_dead_county = (int *)calloc(num_counties, sizeof(int));
	int num_dead_age[18]={0};
	int num_recovered=0; //Number of people recovered. 
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
	int *num_infectious_county;
        num_infectious_county = (int *)calloc(num_counties, sizeof(int));
	int num_infectious_age[18]={0};
	int recovered_hosp=0; //Number of people recovered. 
	int recovered_icu=0; //Number of people recovered.


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

	float infect_prob=0; // Infectious probability
	float community_nom=0; // For adding community infection.
	float community_den=0; // For adding community infection.
	float infect=0; //Infectiousness
	int file_count;

	int num_I=9;
	float Ic=1.0; //Intervention constant for community transmission.
	float (*Iw); //Intervention constant for workplace transmission.
	float Ih=1; //Intervention constant for household transmission.
	float complyI[num_I]; // percent of people who comply with intervention
	float tauI[num_I]; // time after infection that intervention takes place
	float interIc[num_I]; //Intervention constants for Ic
	float interIh[num_I]; //Intervention constants for Ih
	int personinter[num_I]; // Tells us whether the intervention needs to be calculated on a person to person basis or for the whole community. 
	float Ihosp=0.25;  // Accounts for increased cleanliness and infection control at hospital.  
	int * intervene; // 1 if person is currently undergoing interventions, 0 otherwise. 
	intervene = (int*)calloc(population,sizeof(int));
	

	/* current recommendations are intervention 1, 3, and 8. */
	/**** Introduce Interventions: must include documentation for values *****/
	/* No interventions. */
	interIc[0]=1.0;
	float interIw0[6]={0, 1, 1, 1, 1, Ihosp};
	Iw=interIw0;
	interIh[0]=1.00;
	complyI[0]=1.0;
	tauI[0]=0;
	personinter[0]=1;

	/* Intervention 1: school closures, only highschools and colleges.  No school transmission for job_status 3. Assuming no increase in community transmission as students would be working online at the times they would be in college class.*/
	interIc[1]=1.25;
	float interIw1[6]={0, 1, 1, 0, 1, Ihosp};
	interIh[1]=1.50;
	complyI[1]=1.0;
	tauI[1]=0;
	personinter[1]=0;
	
	/* Intervention 2: school closures of all schools. No school transmission for job_status 1, 2, and 3, reduction of 5% in workplace interactions to account for parents becoming childcare.  Children have 50% increase in household transmission and 25% increase in community transmission .*/
	interIc[2]=1.25;
	float interIw2[6]={0, 0, 0, 0, 1.0, Ihosp};
	interIh[2]=1.50;
	complyI[2]=1.0;
	tauI[2]=0;
	personinter[2]=0;

	/* Intervention 3: Case isolation within household. 1 day after symptoms start, 90% comply, householdi contacts remain the same, 25% contact with community, no contact with 
school or workplace. */
	interIc[3]=0.25;
	float interIw3[6]={0, 0, 0, 0, 0, Ihosp};
	interIh[3]=1.0;
	complyI[3]=0.9;
	tauI[3]=6.1;
	personinter[3]=0;

	/* Intervention 4: Case isolation of entire household if one member becomes sick.  Same as case isoloation of single person but now includes all in household.  90% of symptomatic comply and 70% household members comply. */
	interIc[4]=0.25;
	float interIw4[6]={0,0,0,0,0, Ihosp};
	interIh[4]=1.5;
	complyI[4]=0.7;
	tauI[4]=6.1;
	personinter[4]=0;

	/* Intervention 5: Case isolation of entire household if one member becomes sick.  This adds for the case of a quarantined household member getting ill.  tauI=0 */
	interIc[5]=0.25;
	float interIw5[6]={0,0,0,0,0, Ihosp};
	interIh[5]=1.0;
	complyI[5]=0.9;
	tauI[5]=0.0;
	personinter[5]=0;

	/* Intervention 5: social distancing.  workplace contact reduces 25%, household contact increases 25%, community contact reduces 75%. For whole community or subset. 90% comply
	interIc[5]=0.25;
	float interIw5[6]={0, 1.00, 1.00, 1.00, 0.75, Ihosp};
	interIh[5]=1.50;
	complyI[5]=0.90;
	tauI[5]=0;
	personinter[5]=1;
*/
	/* Intervention 6: social distancing with school closure.  Community contacts decrease by 75%, household comntact increase by 25%, 70% compliance.  essential buisnesses stay open, 75% reduction in workplace transmission. NOTE: similar to below except with minimized social interaction. */
	interIc[6]=0.25;
	float interIw6[6]={0,0,0,0,0.25, Ihosp};
	interIh[6]=1.50;
	complyI[6]=0.9;
	tauI[6]=0;
	personinter[6]=1;

	/* Intervention 7: school closures of all schools and non-essential businesses. No school transmission for job_status 1, 2, and 3, reduction of 75% workplace interactions.  50% increase in household transmission and 50% increase in community transmission .*/

	interIc[7]=1.50;
	float interIw7[6]={0, 0, 0, 0, 0.25, Ihosp};
	interIh[7]=1.50;
	complyI[7]=1.0;
	tauI[7]=0;
	personinter[7]=1;

	/* Intervention 8: Social distancing of people over 70 years old. Reduction of 75% workplace interactions. decrease of 75% of community contacts, household contacts increases 25%. 80% comply*/

	interIc[8]=0.25;
	float interIw8[6]={0, 0, 0, 0, 0.25, Ihosp};
	interIh[8]=1.25;
	complyI[8]=0.8;
	tauI[8]=0;
	personinter[8]=0;

	/* Make interIw array.*/
	float *interIw[9]={interIw0, interIw1, interIw2, interIw3, interIw4, interIw5, interIw6, interIw7, interIw8};

	/**** Set random number generator seed. ****/
	COV_init_rand();

	/* Initialize age distribution */
	age_dist(age, population, stats, age_distrib);
	city_lat_long(&num_cities,  lat_city,  long_city, city_names, city_county, county_name, num_counties) ;

	/* Checking data */
	for (i=0; i<num_cities; i++) {
//		printf("cities %i %i %i lat %f lon %f \n", i, city_county[i], num_cities, lat_city[i], long_city[i]);
	}


	/* Parse land_scan file to get population density.  */
        float *lat_locale = NULL, *lon_locale = NULL, *pop_density_init_num = NULL;
        int num_locale = 0, max_locale = 0;
        float tmp_lat, tmp_lon, pop_den;
	FILE* lat_long = fopen("land_pop_sorted.txt", "r"); // Sorted land population in descending order.  Important when we don't have complete population.   
	while ((ret = fscanf(lat_long, "%f%*c%f%*c%f", &tmp_lon, &tmp_lat, &pop_den)) == 3) {
            if (num_locale + 1 > max_locale) {
                max_locale += 10;
                lat_locale = (float *)realloc(lat_locale, max_locale * sizeof(float));
                lon_locale = (float *)realloc(lon_locale, max_locale * sizeof(float));
                pop_density_init_num = (float *)realloc(pop_density_init_num, max_locale * sizeof(float));
            }
            lat_locale[num_locale] = tmp_lat;
            lon_locale[num_locale] = tmp_lon;
            pop_density_init_num[num_locale++] = pop_den;
        }
        fclose(lat_long);
        int **locale_to_HH;
        int *locale_to_HH_n;
        int *locale_HH;
        locale_to_HH = (int **) malloc(num_locale * sizeof(int *));
        for (i = 0; i < num_locale; i++) {
            locale_to_HH[i] = (int *)calloc(num_households, sizeof(int));
        }
        locale_to_HH_n = (int *)calloc(num_locale, sizeof(int));
        locale_HH = (int *)calloc(num_households, sizeof(int));

	/* Initialize households */
	household_lat_long( num_households,  HH,  lat_city, long_city, num_cities, city, county, city_county, city_size, county_size, population, age, per_HH_size, county_name, pop_county, tot_pop, stats, county_p, &num_locale, lat_locale, lon_locale, pop_density_init_num, locale_to_HH, locale_to_HH_n, locale_HH) ;

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

	FILE * lat_lon_out = fopen("lat_lon.dat", "w");
	
	// Infections are randomly placed based on number of initial infections.  //
	// Includes infections from t=-11 to t=-1.
	// Percent per county taken from C19.se infections as of 2020/3/25.
	// Initial infections calculated from population admitted to intensive care per day from 2020/3/14 to 2020/3/24.
	float initial_per[21]={0.4234, 0.0404, 0.0336, 0.0843, 0.0257, 0.0079, 0.0071, 0.0020, 0.00475, 0.0973, 0.0261, 0.1088, 0.0178, 0.0230, 0.0115, 0.0158, 0.0127, 0.0075, 0.0233, 0.0131, 0.0139}; 
	/***** THIS IS THE REAL INITIALIZATION ARRAY, based on ICU numbers, day 0 is 3/26 ******/
	float initialize[15]={1667, 4231, 4181, 4407, 3051, 1808, 2599, 1469, 1695, 339, 678, 791, 678, 339, 113};
	/**** TMP INTIALIZATION ARRAY ***/
//	float initialize[15]={100, 50, 40, 20, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
	float tmp_t;
	fprintf(stats, "Initial Infections by county \n");
        fflush(stats);
	for ( tmp_t=-14; tmp_t<=0; tmp_t++ ) {
		int l=-tmp_t;
		for ( j=0; j<21; j++ ) {
			initial_infections[j]=initial_per[j]*initialize[l];
			fprintf(stats, "time %f county %i initial_infections %i percent %f total_intialized %f \n", tmp_t, j, initial_infections[j], initial_per[j], initialize[l]);
		}
		fprintf(stats, "\n\n");
                fflush(stats);
                /* Randomly assign initial infections */
                initialize_infections( initial_infections,  tau,  infected,  severe,  symptomatic,  county,  &num_infect,  num_counties,  symptomatic_per,  population, dt, tmp_t, lat_locale, lon_locale, num_infect_county, num_infect_age, age, county_p, county_size, locale_HH, HH) ;
	}
        fflush(stats);
        printf("All infections initialized\n");
        fflush(stdout);

	// Uncomment to see initial distribution of infections by county.
	for (i=0; i<num_counties; i++) {
	}



	int ** job_status_county; // Jobs per county, or city for schools
	job_status_county = (int**)calloc(6,sizeof(int*));
	for (i=0;i<6;i++) job_status_county[i] = (int*)calloc(num_cities,sizeof(int)) ;
	/* Initialize job/school status */
	job_dist(job_status, job_status_county, age, county, city, population, num_cities, num_counties, county_size, stats); 

	int hosp_num[num_counties]; //Number of hospitals per county
	memset(hosp_num, 0, num_counties);
	/* Initialize workplace/school */	
	workplace_dist(workplace, job_status, job_status_county, city, num_cities, county, num_counties, population, &max_num_WP , hosp_num, class, stats); 

	free(city);

        fclose(stats);

	/* Get size of each workplace as an array.  Cannot allocate array until max_num_WP is known. */
	int workplace_size[6][max_num_WP];
	memset(workplace_size, 0, 6*max_num_WP*sizeof(int));
	for (i=0; i < population; i++) {
		workplace_size[job_status[i]][workplace[i]]++;
	}
	
	
	double fd_calc[22000];
	/*Precalculate density kernel */
	for (i=0; i<22000; i++) {
		fd_calc[i]=1/(1+pow((((double)i/10.)/4), 3)); //kernel density function as parameterized for GB.	
	}	


        printf("Starting density kernel calculations\n");
        fflush(stdout);
        ret = clock_gettime(CLOCK_MONOTONIC, &t1);
	/* Precalculate total density kernel function for each individual */
	for (i=0; i<num_locale; i++) {
                float itmp_fd;
                int npi; /* number of persons in locale i */
                itmp_fd = 0;
                npi = 0;
                for (int hh = 0; hh < locale_to_HH_n[i]; hh++) {
                    npi += per_HH_size[locale_to_HH[i][hh]];
                }
#ifdef _OPENMP
#pragma omp parallel for private(j) reduction(+:itmp_fd)
#endif
		for (j=i+1; j<num_locale; j++) {
                        double d;
                        double tmp_fd;
                        int npj; /* number of persons in locale j */
                        npj = 0;
                        for (int hh = 0; hh < locale_to_HH_n[j]; hh++) {
                            npj += per_HH_size[locale_to_HH[j][hh]];
                        }
			d=distance(lat_locale[i], lon_locale[i], lat_locale[j], lon_locale[j], 'K');

                        if (full_fd) {
                            tmp_fd = 1/(1+pow((d/4), 3)); //kernel density function as parameterized for GB.
                        } else {
                            tmp_fd = fd_calc[(int)(d*10)]; //kernel density function as parameterized for GB.
                        }
                        itmp_fd += tmp_fd * npj;
			fd_tot[j] += tmp_fd * npi;
		}
                fd_tot[i] += itmp_fd + npi - 1;
	}
        ret = clock_gettime(CLOCK_MONOTONIC, &t2);
        tt = ((double)t2.tv_sec + (double)t2.tv_nsec/nsdiv) - ((double)t1.tv_sec + (double)t1.tv_nsec/nsdiv);
        printf("Done with density kernel calculations in %5.2f\n", tt);
        fflush(stdout);

	/* Initialization complete... start simulation */



	/* Initialize constants */
	float alpha=0.8 ; // From Ferguson Nature 2006
	float omega=2 ; // From Ferguson Nature 2006
	// #Leaving out rho from Ferguson 2006.  This is a measure of how infectious person is.  For now, we will assume all people are the same.

	/* precalculate kappa */
	int count_kappa_vals=ceil(12/dt);
	float kappa_t=0;
	float * kappa_vals;
	float tau1=0;	
	kappa_vals = (float*)calloc(count_kappa_vals,sizeof(float));
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
	

	Ic=1;
	float t=0;
	float time_step=0;

        ret = clock_gettime(CLOCK_MONOTONIC, &T2);
        tt = ((double)T2.tv_sec + (double)T2.tv_nsec/nsdiv) - ((double)T1.tv_sec + (double)T1.tv_nsec/nsdiv);
        printf("Total initialization time %5.2f\n", tt);

	/* Start simulation */			
        printf("Starting simulation\n");
        fflush(stdout);

        ret = clock_gettime(CLOCK_MONOTONIC, &T1);
	for (time_step=0; time_step<(tot_time/dt); time_step++) {
                ret = clock_gettime(CLOCK_MONOTONIC, &t1);
		t=t+dt;
                printf("Timestep %5.2f ", t);
                fflush(stdout);

		/* count origin of contacts */	
		 num_contact_commun=0;
		 num_contact_work=0;
		 num_contact_school=0;	
		 num_contact_house=0;
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
	//		for ( i=0; i<population; i++ ) {
	//			intervene[i]=1;	
	//		}
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
		segment_population( &num_sus,  &num_infectious,  &num_hosp,  &num_icu,  infected,  infectious,  sus_list,  hosp_list,  hosp_pop,  icu_pop,  icu_list,  tau,  population,  t,  num_sus_county,  num_infectious_county,  num_infect_county,  num_hosp_county,  num_icu_county,  num_sus_age,  num_infectious_age,  num_infect_age,  num_hosp_age,  num_icu_age,  age,  county, print_lat_lon,  lat_locale,  lon_locale,  lat_lon_out, locale_HH, HH); 

		//#### Only Susceptible people can get the virus and infected people spread it.
		for (i=0; i<num_sus; i++) {
                        int sus_person; //Counter for susceptible person.
                        int age_group;
			sus_person=sus_list[i];
			community_nom=0;
			community_den=0;
			infect = 0;
			contact_commun=0;
			contact_work=0;
			contact_house=0;
			contact_school=0;
                        age_group=floor(age[sus_person]/5);

#ifdef _OPENMP
#pragma omp parallel for private(j) default(shared) reduction(+:infect) reduction(+:community_nom) reduction(+:contact_commun) reduction(+:contact_work) reduction(+:contact_house) reduction(+:contact_school)
#endif
			for (j=0; j<num_infectious; j++) {


                                int infec_person; //Counter for infected person.
                                float kappa; // #Infectiousness
                                float tIc, tIh, *tIw;
				infec_person=infectious[j];
                                tIc = Ic;
                                tIh = Ih;
                                tIw = Iw;
				/* Determine if person is under individual interventions and set parameters */
				if ( intervene[infec_person] > 0 && t>tau[infec_person]+tauI[intervene[infec_person]]) {
					tIh=interIh[intervene[infec_person]];
					tIc=interIc[intervene[infec_person]];
					tIw=interIw[intervene[infec_person]];
				} 
			//	printf("inter %f %f %f \n", tIc, tIh, tIw[job_status[infec_person]]);	
				/* This will probably have to move outside to a pair list.  NOTE: The list of coworkers/classmates and community members within contact may not completely overlap. i.e. a coworker could be outside of the realm of commumnity transmission if someone lives on the edge of a county. */	
				kappa = calc_kappa( t,  tau[infec_person], symptomatic[infec_person], dt, kappa_vals);

				
				if (hosp_pop[infec_person]==0) {
                                        float d; //distance between people.
					if (HH[sus_person]==HH[infec_person]) {
						// Household transmission //
						infect+=tIh*calc_household_infect(kappa, omega, per_HH_size[HH[sus_person]], alpha, severe[infec_person]); 
						contact_house++;
					}

					// Workplace/School transmission: People must be in same workplace and job type. // 
					if ((workplace[sus_person]==workplace[infec_person]) && (job_status[sus_person]==job_status[infec_person]) && (job_status[sus_person]>0) && tIw[job_status[sus_person]]>0) {
						infect+=calc_workplace_infect(job_status[sus_person], kappa, omega, workplace_size[(job_status[sus_person])][(workplace[sus_person])], severe[infec_person], tIw) ;
						if (job_status[sus_person]<4) {
							contact_school++;
							if (class[sus_person]==class[infec_person]) {
								infect+=calc_workplace_infect(job_status[sus_person], kappa, omega, 15.00, severe[infec_person], tIw) ;
							}
						} else {
							contact_work++;
						}
					}

					// Community transmission // 
                                        d=distance(lat_locale[locale_HH[HH[sus_person]]], lon_locale[locale_HH[HH[sus_person]]], lat_locale[locale_HH[HH[infec_person]]], lon_locale[locale_HH[HH[infec_person]]], 'K');
					community_nom+=tIc*calc_community_infect( age_group, kappa, omega, severe[infec_person], d, fd_calc);
					contact_commun++;
				} else {
					/* In hospital, only have interaction with hospital workers and half interaction with family (household). */
					// Workplace/School transmission: People must be in same workplace and job type. // 
					if ((workplace[sus_person]==workplace_tmp[infec_person]) && (job_status[sus_person]==job_status[infec_person])) {
						infect+=calc_workplace_infect(job_status[sus_person], kappa, omega, workplace_size[(job_status[sus_person])][(workplace[sus_person])], severe[infec_person], tIw) ;
					}
					// Household transmission //
					if (HH[sus_person]==HH[infec_person]) {
						infect+=0.25*calc_household_infect(kappa, omega, per_HH_size[HH[sus_person]], alpha, severe[infec_person]); 
					}

				}
			}
                        infect+=community_nom/fd_tot[locale_HH[HH[sus_person]]]; // Community spread is additive nominator and denominator.  Must be outside of infectious persons loop.


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
				if (contact_work>0) {
					num_contact_work++;
					num_contact_work_county[county[sus_person]]++;
					num_contact_work_age[(int)floor(age[sus_person]/5)]++;
				} 
				if (contact_commun>0) {
					num_contact_commun++;
					num_contact_commun_county[county[sus_person]]++;
					num_contact_commun_age[(int)floor(age[sus_person]/5)]++;
				}
				if (contact_house>0) {
					num_contact_house++;
					num_contact_house_county[county[sus_person]]++;
					num_contact_house_age[(int)floor(age[sus_person]/5)]++;
				}
				if (contact_school>0) {
					num_contact_school++;
					num_contact_school_county[county[sus_person]]++;
					num_contact_school_age[(int)floor(age[sus_person]/5)]++;
				}
				num_infect++;
				num_infect_county[county[sus_person]]++;
				num_infect_age[(int)floor(age[sus_person]/5)]++;
				/* Determine if following interventions only for interventions that effect individuals.*/
				if ( interventions > 0 && COV_rand() < complyI[3] ) {
					if ( interventions == 2 && intervene[sus_person] == 4 ) {
						intervene[sus_person]=5;
					} else {
						intervene[sus_person]=3;
					}
					/* Intervention 2 is household quarantine with current recommendations. Applicable for whole household.  */
					if ( interventions == 2 && t>tauI_onset ) {
						int i1;
						for (i1=0; i1<population; i1++) {
							if ( HH[sus_person] == HH[i1] && COV_rand()<complyI[4] ) {
								intervene[i1]=4;
							}
						}
			
					}
				}
			}	

		}


		// After 5 days of symptoms, people are randomly put into hospital based on age distributions. //
		hosp_entry(t, num_infectious,  infectious,  age,  icu_pop,  hosp_pop,  symptomatic, tau, workplace_tmp, hosp_num, county, dt) ;
		
		// Release from hospital // 
		hosp_release(t, num_hosp,  hosp_list,  tau,  recovered, hosp_pop, &num_recovered, &recovered_hosp, &recovered_icu, dt, num_recovered_county, num_recovered_age, age, county, num_recovered_hosp_county, num_recovered_hosp_age, num_recovered_icu_county, num_recovered_icu_age) ;
		
		
		// 15 days after onset of symptoms, people randomly die based on age based percentages.  ICU patients are removed from ICU into regular hospital. // 
		
		num_dead = death(t, num_infectious,  infectious,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt, num_dead_county, num_dead_age, county) ;
		num_dead = death(t, num_icu,  icu_list,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt, num_dead_county, num_dead_age, county) ;
		num_dead = death(t, num_hosp,  hosp_list,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt, num_dead_county, num_dead_age, county) ;


		// Recovered after 11 days (6 days of symptoms) if not in hospital/ICU. // 
		for (i=0; i<num_infectious; i++) {
                        int infec_person;
			infec_person=infectious[i];
			if ((tau[infec_person]>t-11) && (hosp_pop[infec_person]==0) && (icu_pop[infec_person]==0)) {
				recovered[infec_person]=1;
				num_recovered++;
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
	
	for (file_count=0; file_count<num_counties; file_count++) {
		fprintf(county_files[file_count], "Walltime/timestep %6.2f Time %6.2f num_infected %i num_infectious %i num_in_hosp %i num_in_icu %i num_dead %i recovered_tot %i recovered_from_hosp %i recovered_from_icu %i contact_work %i contact_school %i contact_home %i contact_community %i total_individuals %i \n", step_time, t, num_infect_county[file_count], num_infectious_county[file_count], num_hosp_county[file_count], num_icu_county[file_count], num_dead_county[file_count], num_recovered_county[file_count], num_recovered_hosp_county[file_count], num_recovered_icu_county[file_count], num_contact_work_county[file_count], num_contact_school_county[file_count], num_contact_house_county[file_count], num_contact_commun_county[file_count], county_size[file_count]);
		fflush(county_files[file_count]);
	}

	for (file_count=0; file_count<90; file_count+=5) {
		fprintf(age_files[file_count/5], "Walltime/timestep %6.2f Timestep %6.2f num_infected %i num_infectious %i num_in_hosp %i num_in_icu %i num_dead %i recovered_tot %i recovered_from_hosp %i recovered_from_icu %i contact_work %i contact_school %i contact_home %i contact_community %i total_individuals %i \n", step_time, t, num_infect_age[file_count/5], num_infectious_age[file_count/5], num_hosp_age[file_count/5], num_icu_age[file_count/5], num_dead_age[file_count/5], num_recovered_age[file_count/5], num_recovered_hosp_age[file_count/5], num_recovered_icu_age[file_count/5], num_contact_work_age[file_count/5], num_contact_school_age[file_count/5], num_contact_house_age[file_count/5], num_contact_commun_age[file_count/5], age_distrib[file_count/5]);
		fflush(age_files[file_count/5]);
	}

	}
        ret = clock_gettime(CLOCK_MONOTONIC, &T2);
        step_time = ((double)T2.tv_sec + (double)T2.tv_nsec/nsdiv) - ((double)T1.tv_sec + (double)T1.tv_nsec/nsdiv);
        fprintf(output_file, "Total time %8.3f\n", step_time);
        printf("Total simulation time %8.3f\n", step_time);


	free(HH);
	free(per_HH_size);
	free(lat_city);
	free(long_city);
	free(city_size);
	free(county_size);
	free(age);
	free(county);
	free(job_status);
	free(workplace);
	free(class);
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

}
