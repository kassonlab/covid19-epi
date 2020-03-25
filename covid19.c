#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "common.h"

int i; // Counter
int j; // Counter

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


void age_dist (float * age, int population) {

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
	}
        free(age_dst);


/* Uncomment to test age distribution.  Tested JMG 2020-03-20.
        printf("Testing age distribuition\n");
	int age_dist_test[9]={0};
	for (i=0; i<population; i++) {
		age_dist_test[(int)floor(age[i]/10)]++;
	}	

	for (i=0; i<9; i++) {
		printf("%i %f %f \n", i, (float)age_dist_test[i]/population, age_dist[i]) ; 
	}
        exit(0);
*/
}


void HH_dist(int * HH, float * age, float HH_size, int * per_HH_size, int num_households, int population) {

	int i ; // Counter
	int HH_count=0 ; // Counter
	int HH_person=0 ; // Counter

	/* Initially set all households to -1 */
	for (i=0; i < population; i++) {
		HH[i]=-1;
	}

	/* First add adult head of household to each household. */
	while ( HH_count < num_households ) {
		while ( age[HH_person]<20 ) {
			HH_person++;
		}
		HH[HH_person]=HH_count;
		HH_person++;
		HH_count++;
	}


	/* Distribute remaining people randomly.  This could be changed to a distribution to more realistically reflect household size in the future. */
	for ( HH_person=0; HH_person<population ; HH_person++) {
		if (HH[HH_person]==-1) {
			HH[HH_person]=(int)(COV_rand() * num_households);
		}
	}

	/* Get size of each household as an array */
	for (HH_person=0; HH_person<population; HH_person++) {
		per_HH_size[HH[HH_person]]++;
	}

/* Uncomment to test household distribution.  Tested JMG 2020-03-20.
	int HH_dist_test[9]={0};
	for (i=0; i<num_households; i++) {
		HH_dist_test[(int)(per_HH_size[i])]++;
	}	

	for (i=0; i<9; i++) {
		printf("%i %f \n", i, HH_dist_test[i]/(float)num_households) ; 
	}
*/
}


/* Puts households in certain locations based on population density distribution.  Only really correct for full population but works for smaller populations.  Biases towards smaller households for smaller populations.  Each household is also fed into a locality based on the shortest distance to the center of that locality for purposes of school and workplace choice. */
void household_lat_long(int num_households, int * HH, float * lat, float * lon, float * lat_city, float * long_city, int num_cities, int * city, int * county, int * city_county, int * city_size, int * county_size, int population, float * age, int * per_HH_size, int * city_int, char ** county_names, float * county_pop, float tot_pop_actual) {

	/* Get longitude and latitude with population density from CSV */
	FILE* lat_long = fopen("land_pop_sorted.txt", "r"); // Sorted land population in descending order.  Important when we don't have complete population.   

	/* Initialize helpers for population density */
	int num_km=206959;
	// FOR TEST TEXT 
	//int num_km=63;
	float pop_density_init_num[num_km];
	float tot_pop_density=0;
	float tmp_lat;
	float tmp_lon;

	/* Initialize helpers for household distribution */
	float lat_HH[num_households];	
	float lon_HH[num_households];	
	int city_HH[num_households];	
	int county_HH[num_households];	
	int county_num;	
	int city_num;	
	
	float dist1;
	int tmp_city;	
	float min_dist=1000000;
	int placement;

	int tmp_county_count[21]={0};	
	int tmp_county_density[21]={0};	

	/* find how many locales will have households */
	int num_locale=0;
	if (num_households<num_km) {
		num_locale=num_households;
	} else {
		num_locale=num_km;
	}

	/* Fill up households. */
	int HH_count=0 ; // Counter
	int HH_person=0 ; // Counter

	/* Initially set all households to -1 */
	for (i=0; i < population; i++) {
		HH[i]=-1;
	}

	/* Save list of households in each locale. */
	int ** county_list;
	county_list = (int**)calloc(num_locale,sizeof(int*));
	for (i=0;i<num_locale;i++) county_list[i] = (int*)calloc(num_households,sizeof(int));
	int * locale_count;
	locale_count = (int*)calloc(num_locale,sizeof(int));
	int * locale_HH_count;
	locale_HH_count = (int*)calloc(num_locale,sizeof(int));
	

	/* Parse land_scan file to get population density.  */
	while (( HH_count < num_households ) && (HH_count < num_km)) {
		min_dist=1000000;
		int got = fscanf(lat_long, "%f%*c%f%*c%f", &tmp_lon, &tmp_lat, &pop_density_init_num[HH_count]);
  		if (got != 3) break; // wrong number of tokens - maybe end of file

		tot_pop_density+=pop_density_init_num[HH_count];

		// Determine city of each population square.  Use city data to determine which schools students attend.  Workplaces are placed by county. //
		for (j=0; j<num_cities; j++) {
			dist1=distance(tmp_lat, tmp_lon, lat_city[j], long_city[j], 'K');	
			if (dist1<min_dist) {
				min_dist=dist1;
				tmp_city=j;
			} 
		}

		city_num=tmp_city;
		county_num=city_county[tmp_city];
		tmp_county_count[county_num]++; 
		tmp_county_density[county_num]+=pop_density_init_num[HH_count]; 

		/* Set up household */
		lat_HH[HH_count]=tmp_lat;	
		lon_HH[HH_count]=tmp_lon;
		city_HH[HH_count]=city_num;
		county_HH[HH_count]=county_num;
		county_list[HH_count][locale_count[HH_count]]=HH_count;
		locale_HH_count[HH_count]+=1;
		locale_count[HH_count]+=1;
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
	/* First add adult head of household to each household. */
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
		locale_HH_count[placement]+=1;
		locale_count[placement]+=1;
		placement+=1;

		/* Set up head of household. */
		HH[HH_person]=HH_count;
		lat[HH_person]=lat_HH[HH_count];	
		lon[HH_person]=lon_HH[HH_count];	
		city[HH_person]=city_HH[HH_count];	
		county[HH_person]=county_HH[HH_count];	
		city_size[city[HH_person]]++;
		county_size[county[HH_person]]++;
		HH_person++;
		HH_count++;
	}

	/* Distribute remaining people randomly.  This could be changed to a distribution to more realistically reflect household size in the future. */
	for ( HH_person=0; HH_person<population ; HH_person++) {
		if (HH[HH_person]==-1) {
			
			/* Place people in random households within locale until locale is full. */
			if (placement>=num_locale) {
				placement=0;
			} else if (locale_count[placement]>=pop_density_init_num[placement]) {
					placement=0; // Start over at top of list.  Remember list is sorted by biggest to smallest locales.
			}	

			/* Pick a random household in the locale. */
			HH[HH_person]=county_list[placement][(int)(COV_rand() * locale_HH_count[placement])];
			lat[HH_person]=lat_HH[HH[HH_person]];	
			lon[HH_person]=lon_HH[HH[HH_person]];	
			city[HH_person]=city_HH[HH[HH_person]];	
			county[HH_person]=county_HH[HH[HH_person]];	
			city_size[city[HH_person]]++;
			county_size[county[HH_person]]++;
			locale_count[placement]+=1;
			placement+=1;
		}
	}

	int max_HH_size=0;
	/* Get size of each household as an array */
	for (HH_person=0; HH_person<population; HH_person++) {
		per_HH_size[HH[HH_person]]++;
		if (per_HH_size[HH[HH_person]]>max_HH_size) {
			max_HH_size=per_HH_size[HH[HH_person]];
		}
	}


/* Uncomment to test household distribution.  Tested JMG 2020-03-22.
	int HH_dist_test[21]={0};
	for (i=0; i<num_households; i++) {
		HH_dist_test[per_HH_size[i]]++;
	}	

	for (i=0; i<max_HH_size; i++) {
		printf("%i %f %i \n", i, HH_dist_test[i]/(float)num_households, num_households) ; 
		fflush(stdout);
	}
	
	int counter1[num_cities];
	memset(counter1, 0, num_cities*sizeof(int));
	for (i=0; i<population; i++) {
		counter1[city[i]]++;
	}

	for (i=0; i<num_cities; i++) {
		printf("city_dist %i num %i percent %f actual %f \n", i, city_size[i], city_size[i]/(float)population, pop_density_init_num[i]/tot_pop_actual)  ;
		fflush(stdout);
	} 

	for (i=0; i<21; i++) {
		printf("county_dist %i num %i percent %f actual %f \n", i, county_size[i], county_size[i]/(float)population, county_pop[i]/tot_pop_actual)  ;
		fflush(stdout);
	} 
*/
	free(locale_HH_count);
	free(county_list);
	free(locale_count);
	fclose(lat_long);
}

void city_lat_long(int *num_cities, float * lat_city, float * long_city, char ** cities, int * city, int * county, char ** county_names, int num_county) {

	/* Get longitude and latitude of cities in Sweden from CSV */
	FILE* fp = fopen("cities_all.csv", "r");  // Not sure about the validity of this file.  Could use a better source.

	/* Counters for parsing file. */
	char tmp[199]; // municipality name
	char tmp1[199]; // locality name
	char tmp2[199]; // county name
	float tmp_lat; // tmp latitude
	float tmp_lon; // tmp longitude

	for (i=0; i < 2000; i++) {
		int got = fscanf(fp, "%[^,\n]%*c%[^,\n]%*c%[^,\n]%*c%f%*c%f\n", tmp1, tmp, tmp2, &tmp_lat, &tmp_lon);
  		if (got != 5) break; // wrong number of tokens - maybe end of file
			
		/* Save city name, longitude, latitude, and increase number of cities. Looking at unique municipalities.*/
		cities[*num_cities]=tmp1;
		city[*num_cities]=*num_cities;
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


void job_dist(int * job_status, int ** job_status_city, float * age, int * county, int * city, int population, int num_cities, int num_counties) {

	int i; 
	

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
		} else if (age[i]>=22 && age[i]<=75) {
			if (COV_rand() < 0.734) {
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
		}
	}

/* CITY VERSION NOT IN USE All currently randomly placed based on county.  Would like to do per town but each town needs inhabitants. 
	for (i=0; i < population; i++) {
		if (age[i]<1 || age[i]>75) {
			job_status[i]=0;
			job_status_city[0][city[i]]++;
		} else if (age[i]>=1 && age[i]<6) {
			if (COV_rand() < 0.9000) {
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
			job_status_city[3][city[i]]++;
		} else if (age[i]>=22 && age[i]<=75) {
			if (COV_rand() < 0.734) {
				job_status[i]=4;
				job_status_city[4][county[i]]++; // Workplace is based on county, not city.
			} else {
				job_status[i]=0;
				job_status_city[0][city[i]]++;
			}
		}
	}
	*/

/* Uncomment to test job distribution.  Tested JMG 2020-03-20. 
	int job_dist_test[5]={0};
	int city_dist_test[5][2000]={0};
	for (i=0; i<population; i++) {
		job_dist_test[job_status[i]]++;
		city_dist_test[job_status[i]][county[i]]++;
	}	

	j=0;
	for (i=0; i<5; i++) {
		for (j=0; j<num_counties; j++) {
	//		printf("jobs job_status %i county %i per_job_total %f num_jobs_in_county %i \n", i, j, job_dist_test[i]/(float)population, city_dist_test[i][j]) ; 
		}
	}

*/	

}

void workplace_dist(int * workplace, int * job_status, int ** job_status_county, int * city, int num_cities, int * county, int num_counties, int population, int * max_num_WP , int * hosp_num, int* class) {

	int pp_class = 19; //Assumption of 15 children per class.
	int pp_preschool = 53; //Assumption of 200 children per school.
	int pp_school = 220; //Assumption of 200 children per school.
	int pp_hospital = 120; //Assumption of 120 people per hospital.
	int pp_work = 15; //Assumption of 15 people per close work group.
	int i;
	int j;
	int num_workplaces[6][num_counties];
	memset(num_workplaces, 0, 6*num_counties*sizeof(int));
	int num_workplaces2[6];
	memset(num_workplaces2, 0, 6*sizeof(int));

	for (i=0; i < num_counties; i++) {
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
		for (j=2; j<4; j++) {
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

	j=4; // Broken down in case we want to do schools by municipality.
	for (i=0; i < num_counties; i++) {
		if (job_status_county[4][i]>0) {
			num_workplaces[4][i]=ceil(job_status_county[4][i]/(float)pp_work);
			num_workplaces2[4]+=ceil(job_status_county[4][i]/(float)pp_work);
//			printf("work %i %i %i %i \n", i, job_status_county[4][i], num_workplaces2[4], num_workplaces[4][i]);
			if (num_workplaces2[4]>*max_num_WP) {
				*max_num_WP=num_workplaces2[4];
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

	for (i=0; i < population; i++) {
		//Try to minimize necessary memory by making job_numbers independent with job_status. //
		int prior_workplaces=0;
		for (j=0; j<county[i]; j++) {
			prior_workplaces+=num_workplaces[job_status[i]][j];
		}
		if (((num_workplaces[job_status[i]][county[i]])>0) && (job_status[i]<4)) {
			workplace[i]=(int)(COV_rand() * (num_workplaces[job_status[i]][county[i]]))+prior_workplaces;
			class[i]=(int)(COV_rand() * (ceil(num_workplaces[job_status[i]][county[i]]/(float)pp_class)));
		} else if ((num_workplaces[job_status[i]][county[i]])>0) {
			workplace[i]=(int)(COV_rand() * (num_workplaces[job_status[i]][county[i]]))+prior_workplaces;
		} else {
			workplace[i]=0;
			printf("no workplaces %i %i %i \n", i, job_status[i], county[i]);
		} 
			
	}

}

float calc_kappa(float t, float tau, int symptomatic) {

	float kappa;
	float t1;
	//###Determine kappa for infected person.  This is the infectiousness of the person based on time since infection started.  Latency period is 4.6 days.  Infection starts at 5.1 days and lasts for 6 days.  Sympotmatic people are twice as likely to infect others as asymptomatic.
	// Kappa is a log normal function with mean of -0.72 and standard deviation of 1.8.  From Ferguson Nature 2005
	if (t-tau<4.6) {
		kappa=0.;
	} else if (t-tau>11.1) {
		kappa=0.; //# Recovered or dead
	} else {
		t1=(log(t-tau-4.6)+0.72)/1.8;
		kappa=exp(-0.5*pow(t1,2))/((t-tau-4.6)*1.8*sqrt(2*pi));
	}
	if (symptomatic==0) {
		kappa=kappa*0.5;
	}
	return(kappa);
}

int * initialize_infections(int * initial_infections, float * tau, int * infected, int * severe, int * infected_list, int * symptomatic, int * county, int * num_infect, int num_counties, float symptomatic_per, int population, float dt, float t) {

	int person_infected=0;
	int tmp_infect=0;
	int i;

		
	for (i=0; i < num_counties; i++) {
		tmp_infect=0;
		while ((tmp_infect<initial_infections[i])) {
			person_infected=(int)(COV_rand() * population);

			if ((county[person_infected]==i) && (infected[person_infected]==0)) {
	
				infected[person_infected]=1;
				severe[person_infected]=(int)(COV_rand() * 2);
	
				if (COV_rand() < symptomatic_per) {
					symptomatic[person_infected]=1;
				}
		//		tau[person_infected]=-COV_rand() * 5;
				tau[person_infected]=t;
				infected_list[*num_infect]=person_infected;
				*num_infect=*num_infect+1;
				tmp_infect++;
			}
			
		}
	}
	return(infected_list);
}


void segment_population(int* num_sus, int* num_infectious, int* num_hosp, int* num_icu, int* infected, int* infectious, int* sus_list, int* hosp_list, int * hosp_pop, int * icu_pop, int* icu_list, float* tau, int population, float t) {

	*num_sus=0;
	*num_infectious=0;
	*num_hosp=0;
	*num_icu=0;

//	int i;
//	int j;
	for (i=0; i<population; i++) {
		if (infected[i]==0) {
			sus_list[*num_sus]=i;
			*num_sus=*num_sus+1;
		} else if ((tau[i]<t-4.6) && (tau[i]>t-11.1)) {
			infectious[*num_infectious]=i;
			*num_infectious=*num_infectious+1;
		}

		if (icu_pop[i]==1) {
			icu_list[*num_icu]=i;
			*num_icu=*num_icu+1;
		} 
		if (hosp_pop[i]>0) {
			hosp_list[*num_hosp]=i;
			*num_hosp=*num_hosp+1;
		}
	}

/* Uncomment for information */
//	printf("Time %i infections i %i num_sus %i infected[i] %i num_infectious %i num_hosp %i num_icu %i \n", t, i, *num_sus, infected[i], *num_infectious, *num_hosp, *num_icu);

}

float calc_household_infect(float kappa, float omega, int HH_size, float alpha, int severe) {

	float betah=0.627; // Scaled from betah=0.4 in influenza pandemic with R0=1.6, COVID-19 R0=2.4 (Ferguson 2020)

	return(betah*kappa*(1+(float)severe*(omega-1))/(pow((float)HH_size,alpha))); 

}

float calc_workplace_infect(int job_status, float kappa, float omega, int workplace_size, int severe, float * Iw) {

	float betap[]={0.0, 1.254, 1.254, 1.254, 0.627, 0.15627} ; // Spread in all types of schools (preschool to college) is twice that of workplace 
	float psi[]={0.0, 0.1, 0.2, 0.25, 0.5, 0.5} ; // Accounts for absenteeism based on severe illness. Ferguson Nature 2006

	return(Iw[job_status]*betap[job_status]*kappa*(1+(float)severe*(omega*psi[job_status]-1))/((float)workplace_size));
}

float calc_community_infect(int age_group, float kappa, float omega, int severe, float d, float *community_den) {

	/* need to work on this.  Perhaps we take a random distance for each two people based on population density, number of people in county, county area, etc. */
	float zeta[]={0.1, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.75, 0.50, 0.25, 0.25, 0.25} ; //   # Travel related parameter for community transmission. Ferguson Nature 2006
	float fd;
	float betac=0.1 ; // Scaled from betac=0.075 in influenza pandemic with R0=1.6, COVID-19 R0=2.4 (Ferguson 2020)

	fd=1/(1+pow((d/4), 3)); //kernel density function as parameterized for GB.
	*community_den=*community_den+fd;
	return(zeta[age_group]*betac*kappa*fd*(1+severe*(omega-1)));
}

void hosp_entry(float t, int num_infectious, int * infectious, float * age, int * icu_pop, int * hosp_pop, int * symptomatic, float * tau, int * workplace_tmp, int * workplace_num, int * county, float dt) {

	float hosp[]={0.001, 0.003, 0.012, 0.032, 0.049, 0.102, 0.166, 0.243, 0.273}; //# From Ferguson 2020
	float icu[]={0.05, 0.05, 0.05, 0.05, 0.063, 0.122, 0.274, 0.432, 0.709} ; //# percent of hospitalized that need icu from Ferguson 2020.
	int age_group; 
	int i;
	int j;
	int infec_person;

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


void hosp_release(float t, int num_hosp, int * hosp_list, float * tau, int * recovered, int * hosp_pop, int * num_recovered, int * recovered_hosp, int * recovered_icu, float dt) {

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
		} else if ((tau[infec_person]>=t-25) && (tau[infec_person]<t-25+dt) && (hosp_pop[infec_person]==2)) {
			hosp_pop[infec_person]=0;
			recovered[infec_person]=1;
			*num_recovered=*num_recovered+1;
			*recovered_icu=*recovered_icu+1;
		}
	}
}

int death(float t, int num_infectious, int * infectious, float * tau, int * dead, int * icu_pop, int * hosp_pop, int * symptomatic, int num_dead, float * age, float dt) {

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
				}
			} else if (symptomatic[infec_person]==1) {	
				age_group=floor(age[infec_person]/10);
				if (COV_rand() < fatal_symptomatic[age_group]) {
					dead[infec_person]=1;
					icu_pop[infec_person]=0;
					hosp_pop[infec_person]=0;
					num_dead++;
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
	int tot_time=500; // Simulation time. 
	int initial_infections[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; //Initial infections per county.
	float percent_infect=0.001 ; // Default 1% of population initially infected.
	float dt=1.00; // Time step.
	int interventions=0; // Value of interventions.
	float tauI_onset=0; //time after start of simulation that interventions for whole community take place.
        int ret;

        /* Timing variables */
        struct timespec T1, T2, t1, t2;
        double step_time, nsdiv = 1000*1000*1000;
	
	/* Parse command-line arguments */
  	for (i=1;i<argc;i++) {
    		if (!strcmp(argv[i],"-pop")) population=atoi(argv[++i]);
    		else if (!strcmp(argv[i],"-sim_time")) tot_time=atoi(argv[++i]);
    		else if (!strcmp(argv[i],"-infect")) percent_infect=atof(argv[++i]);
    		else if (!strcmp(argv[i],"-dt")) dt=atof(argv[++i]);
    		else if (!strcmp(argv[i],"-inter")) interventions=atoi(argv[++i]);
    		else if (!strcmp(argv[i],"-tauI")) tauI_onset=atof(argv[++i]);
    		else if (!strcmp(argv[i],"-initial")) for (j=0; j<21; j++) initial_infections[j]=atof(argv[++i]);
  	}

	printf("here");
	fflush(stdout);
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
//	char * county_name[]={"Stockholm/Uppsala", "Sodermanland", "Ostergotland", "Jonkoping", "Kronoberg", "Kalmar", "Gotland", "Blekinge", "Skane", "Halland", "Vastra Gotaland", "Varmland", "Orebro", "Vastmanland", "Dalarna", "Gavleborg", "Vasternorrland", "Jamtland", "Vasterbotten", "Norrbotten"}; 
//	float pop_county[]={2760794, 297540, 465495, 363599, 201469, 245446, 59686, 159606, 1377827, 333848, 1725881, 282414, 304805, 275845, 287966, 287382, 245347, 130810, 271736, 250093};
	//Combine Stockholm/Uppsala for purposes of job allocating.
	int num_counties=(sizeof(pop_county)/sizeof(pop_county[0]));
	int county_int[num_counties]; // Integer values for random probability distribution. 
	memset(county_int, 0, num_counties*sizeof(int));
	float * lat; //latitude of person
	lat = (float*)calloc(population,sizeof(float));
	float * lon; //longitude of person
	lon = (float*)calloc(population,sizeof(float));
	float * lat_city; //latitude of city i.
	lat_city = (float*)calloc(2000,sizeof(float));
	float * long_city; // longitude of city i.
	long_city = (float*)calloc(2000,sizeof(float));
	int * city_size; // populatino of city i.
	city_size = (int*)calloc(2000,sizeof(int));
	int * county_size; // population of county i.
	county_size = (int*)calloc(num_counties,sizeof(int));
	char * city_names[2000]; //City name of city i.
	char * city_county_names[2000];
	int * city; // city of person i by integer assignment.
	city = (int*)calloc(population,sizeof(int));
	int * city_county; // county of person i by integer assignment.
	city_county = (int*)calloc(2000,sizeof(int));
	int num_cities=0; // total number of cities
	int city_int[2000]; // Integer values for random probability distribution. 
	memset(city_int, 0, 2000*sizeof(int));


	printf("here");
	fflush(stdout);
	/* City information for allocating schools. */
	
	/* Initialize and allocate arrays */
	float * age;  // Age of population
	age = (float*)calloc(population,sizeof(float));
	int * county;  // County of each inhabitant
	county = (int*)calloc(population,sizeof(int));
	int * job_status; // Type of job each person holds: 0-4 for no job, preschool, elementary school, highschool/college, and job, respectively. 	
	job_status = (int*)calloc(population,sizeof(int));
	int * workplace; // Workplace of each person.
	workplace = (int*)calloc(population,sizeof(int));
	int * class; // Classroom for students 	
	class = (int*)calloc(population,sizeof(int));

	int max_num_WP=0; // max workplaces per job_status for allocating array.

	/* Parameters for infections */
	int num_infections=(int)population*percent_infect; // Default is 10% of population has illness.
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
	int num_infect=0; //Number infected.
	int num_icu=0; //Number in ICU.
	int num_hosp=0; //Number in hospital. 
	int num_dead=0; //Number of deaths. 
	int num_recovered=0; //Number of people recovered. 
	int num_infectious=0; //Number infectious.
	int recovered_hosp=0; //Number of people recovered. 
	int recovered_icu=0; //Number of people recovered.
	/* This list hold indices for each category. */ 
	int * icu_list; //Indices in ICU.
	icu_list = (int*)calloc(population,sizeof(int));
	int * hosp_list; //Indices in hospital.
	hosp_list = (int*)calloc(population,sizeof(int));
	int * infected_list; //Indices of infected.
	infected_list = (int*)calloc(population,sizeof(int));
	int * infectious; //Indices of infected.
	infectious = (int*)calloc(population,sizeof(int));
	int * sus_list; //Indices of susceptible.
	sus_list = (int*)calloc(population,sizeof(int));

	printf("here");
	fflush(stdout);
	/* Parameters pertaining to the population.  These could be made into a struct later to look nicer. */
	int sus_person; //Counter for susceptible person.
	int infec_person; //Counter for infected person.
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

	int HH_tmp;
	int age_group; 
	float infect_prob=0; // Infectious probability
	float d; //distance between people.
	float community_nom=0; // For adding community infection.
	float community_den=0; // For adding community infection.
	float infect=0; //Infectiousness
	float ran_num;

	int num_I=7;
	float Ic=1.0; //Intervention constant for community transmission.
	float (*Iw); //Intervention constant for workplace transmission.
	float Ih=1; //Intervention constant for household transmission.
	float complyI[num_I]; // percent of people who comply with intervention
	float tauI[num_I]; // time after infection that intervention takes place
	float interIc[num_I]; //Intervention constants for Ic
	float interIh[num_I]; //Intervention constants for Ih
	int personinter=0; // Tells us whether the intervention needs to be calculated on a person to person basis or for the whole community. 
	float Ihosp=0.25;  // Accounts for increased cleanliness and infection control at hospital.  
	int * intervene; // 1 if person is currently undergoing interventions, 0 otherwise. 
	intervene = (int*)calloc(population,sizeof(int));
	

	printf("here");
	fflush(stdout);
	/**** Introduce Interventions: must include documentation for values *****/
	/* No interventions. */
	interIc[0]=1.0;
	float interIw0[6]={0, 1, 1, 1, 1, Ihosp};
	Iw=interIw0;
	interIh[0]=1.00;
	complyI[0]=1.0;
	tauI[0]=0;
	personinter=0;

	/* Intervention 1: school closures, only highschools and colleges.  No school transmission for job_status 3. Assuming no increase in community transmission as students would be working online at the times they would be in college class.*/
	interIc[1]=1.25;
	float interIw1[6]={0, 1, 1, 0, 1, Ihosp};
	interIh[1]=1.50;
	complyI[1]=1.0;
	tauI[1]=0;
	personinter=1;
	
	/* Intervention 2: school closures of all schools. No school transmission for job_status 1, 2, and 3, reduction of 5% in workplace interactions to account for parents becoming childcare.  Children have 50% increase in household transmission and 25% increase in community transmission .*/
	interIc[2]=1.25;
	float interIw2[6]={0, 0, 0, 0, 1.0, Ihosp};
	interIh[2]=1.50;
	complyI[2]=1.0;
	tauI[2]=0;
	personinter=0;

	/* Intervention 3: Case isolation within household. 1 day after symptoms start, 70% comply, householdi contacts remain the same, 25% contact with community, no contact with 
school or workplace. */
	interIc[3]=0.25;
	float interIw3[6]={0, 0, 0, 0, 0, Ihosp};
	interIh[3]=1.0;
	complyI[3]=0.7;
	tauI[3]=1.5;
	personinter=0;

	/* Intervention 4: Case isolation of entire household if one member becomes sick.  Same as case isoloation of single person but now includes all in household but only 50% comply. */
	interIc[4]=0.25;
	float interIw4[6]={0,0,0,0,0, Ihosp};
	interIh[4]=2.0;
	complyI[4]=0.5;
	tauI[4]=1.5;
	personinter=0;

	/* Intervention 5: social distancing.  workplace contact reduces 25%, household contact increases 25%, community contact reduces 75%. For whole community or subset. 70% comply*/
	interIc[5]=0.25;
	float interIw5[6]={0, 1.00, 1.00, 1.00, 0.75, Ihosp};
	interIh[5]=1.25;
	complyI[5]=0.70;
	tauI[5]=0;
	personinter=1;

	/* Intervention 6: social isolation of everyone.  Community contacts decrease by 75%, household comntact increase by 25%, 70% compliance.  essential buisnesses stay open, 75% reduction in workplace transmission. NOTE: similar to below except with minimized social interaction. */
	interIc[6]=0.25;
	float interIw6[6]={0,0,0,0,0.25, Ihosp};
	interIh[6]=1.5;
	complyI[6]=0.7;
	tauI[6]=0;
	personinter=1;

	/* Intervention 7: school closures of all schools and non-essential businesses. No school transmission for job_status 1, 2, and 3, reduction of 75% workplace interactions.  50% increase in household transmission and 25% increase in community transmission .*/

	interIc[7]=1.50;
	float interIw7[6]={0, 0, 0, 0, 0.25, Ihosp};
	interIh[7]=1.50;
	complyI[7]=1.0;
	tauI[7]=0;
	personinter=1;

	/* Intervention 8: School closures of all schools and non-essential businesses combined with social distancing. No school transmission for job_status 1, 2, and 3, reduction of 75% workplace interactions. decrease of 75% of community contacts, household contacts increases 25%. NOT IN USE*/

	interIc[8]=0.25;
	float interIw8[6]={0, 0, 0, 0, 0.25, Ihosp};
	interIh[8]=1.25;
	complyI[8]=0.7;
	tauI[8]=0;
	personinter=1;

//	printf("interIc0 %f ", interIc[0]);

	/* Make interIw array.*/
	float *interIw[9]={interIw0, interIw1, interIw2, interIw3, interIw4, interIw5, interIw6, interIw7, interIw8};

	/**** Set random number generator seed. ****/
	COV_init_rand();

	printf("here");
	fflush(stdout);
	/* Initialize age distribution */
	age_dist(age, population);

	city_lat_long(&num_cities,  lat_city,  long_city, city_names, city_int, city_county, county_name, num_counties) ;

	printf("here");
	/* Checking data */
	for (i=0; i<num_cities; i++) {
//		printf("cities %i %i %i %i lat %f lon %f \n", i, city_int[i], city_county[i], num_cities, lat_city[i], long_city[i]);
	}

	printf("here");
	/* Initialize households */
	household_lat_long( num_households,  HH,  lat,  lon, lat_city, long_city, num_cities, city, county, city_county, city_size, county_size, population, age, per_HH_size, city_int, county_name, pop_county, tot_pop) ;


	/* Evenly distribute infection by population */
	double *pop_percent; // Percent for random probability distribution.
        pop_percent = (double *)malloc(num_counties * sizeof(double));
        double *pop_prob_dist;
	float tot=0;
	for (i=0; i < num_counties; i++) {
		if (county_size[i]>0) {
			pop_percent[i]=county_size[i]/(float)population;
		} else {
			pop_percent[i]=0;
		}
		tot+=county_size[i];
	}

	printf("here");
	fflush(stdout);
	// Infections are randomly placed based on number of initial infections.  //
	// Includes infections from t=-11 to t=-1.
	// Percent per county taken from C19.se infections as of 2020/3/25.
	// Initial infections calculated from population admitted to intensive care per day from 2020/3/14 to 2020/3/24.
	float initial_per[21]={0.4234, 0.0404, 0.0336, 0.0843, 0.0257, 0.0079, 0.0071, 0.0020, 0.00475, 0.0973, 0.0261, 0.1088, 0.0178, 0.0230, 0.0115, 0.0158, 0.0127, 0.0075, 0.0233, 0.0131, 0.0139}; 
//	float initialize[11]={4744, 5000, 3462, 2051, 2949, 1667, 1923, 385, 769, 897, 769};
//	float initialize[11]={0,0,0,0,0,4,0,0,0,0,0};
	float initialize[11]={0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0};
	float tmp_t;
	for ( i=1; i<12; i++ ) {
		tmp_t = -i;
		for ( j=0; j<21; j++ ) {
			initial_infections[j]=initial_per[j]*initialize[i-1]*population/tot_pop;
			printf("initial_infect %i %i %i %f \n", j, initial_infections[j], county_size[j], pop_percent[j]);
		}
			/* Randomly assign initial infections */
			infected_list = initialize_infections( initial_infections,  tau,  infected,  severe,  infected_list,  symptomatic,  county,  &num_infect,  num_counties,  symptomatic_per,  population, dt, tmp_t) ;
	}		

/* Old initialization of infection
	int county_count=0;
	for (i=0; i<num_counties; i++) {
		pop_percent[i]=(double)county_size[i]/population;
		county_int[i]=i;
	}

        ret = generate_inc_distr_vec(&pop_prob_dist, pop_percent, num_counties, "population per county");
        if (ret) {
            fprintf(stderr, "Bailing out due to county distribution errors\n");
            exit(1);
        }
	if (initial_infections[0]==0) {
		for (i=0; i < num_infections; i++) {
			j=rand_distr_indx(pop_prob_dist, 0, num_counties-1);
			initial_infections[j]++;
		}
	}
*/
	// Uncomment to see initial distribution of infections by county.
	for (i=0; i<num_counties; i++) {
	//	printf("initial_infect %i %i %i %f \n", i, initial_infections[i], county_size[i], pop_percent[i]);
	}

	int ** job_status_county; // Jobs per county
	job_status_county = (int**)calloc(6,sizeof(int*));
	for (i=0;i<6;i++) job_status_county[i] = (int*)calloc(num_counties,sizeof(int)) ;
	/* Initialize job/school status */
	job_dist(job_status, job_status_county, age, county, city, population, num_cities, num_counties); 

	int hosp_num[num_counties]; //Number of hospitals per county
	memset(hosp_num, 0, num_counties);
	/* Initialize workplace/school */	
	workplace_dist(workplace, job_status, job_status_county, city, num_cities, county, num_counties, population, &max_num_WP , hosp_num, class); 

	/* Get size of each workplace as an array.  Cannot allocate array until max_num_WP is known. */
	int workplace_size[6][max_num_WP];
	memset(workplace_size, 0, 6*max_num_WP*sizeof(int));
	for (i=0; i < population; i++) {
		workplace_size[job_status[i]][workplace[i]]++;
	}


	/* Initialization complete... start simulation */

	/* Seed infections */


	/* Initialize constants */
	float kappa=4 ; // #Infectiousness
	float alpha=0.8 ; // From Ferguson Nature 2006
	float omega=2 ; // From Ferguson Nature 2006
	// #Leaving out rho from Ferguson 2006.  This is a measure of how infectious person is.  For now, we will assume all people are the same.
	
	int contact_commun=0;
	int contact_work=0;
	int contact_school=0;	
	int contact_house=0;
	int num_contact_commun=0;
	int num_contact_work=0;
	int num_contact_school=0;	
	int num_contact_house=0;
	float infect_commun=0;
	float infect_house=0;
	float infect_work=0;

	Ic=1;
	float t=0;
	float time_step=0;
	/* Start simulation */			
        ret = clock_gettime(CLOCK_MONOTONIC, &T1);
	for (time_step=0; time_step<(tot_time/dt); time_step++) {
                ret = clock_gettime(CLOCK_MONOTONIC, &t1);
		t=t+dt;

		/* count origin of contacts */	
		 num_contact_commun=0;
		 num_contact_work=0;
		 num_contact_school=0;	
		 num_contact_house=0;

		/* Introduce overall community interventions. */
		if (interventions == 0)  {
			Ic=interIc[0];
			Ih=interIh[0];
			Iw=interIw[0];
	//		printf("here2 %f %f %f %f %f %f %f %f %f %i \n", Iw[0], Iw[1], Iw[2], Iw[3], Iw[4], Iw[5], tauI_onset, Ic, Ih, interventions);
		} else if (personinter>0 && t >= tauI_onset && t <= tauI_onset+dt) {
			Ic=interIc[interventions];
			Ih=interIh[interventions];
			Iw=interIw[interventions];
			for (i=0; i<population; i++) {
				intervene[i]=round(COV_rand()*complyI[interventions]);
			}
	//		printf("here1 %f %f %f %f %f %f %f %f %f %i \n", Iw[0], Iw[1], Iw[2], Iw[3], Iw[4], Iw[5], tauI_onset, Ic, Ih, interventions);
		} else if (interventions == 2 && t >= tauI_onset && t <= tauI_onset+dt) {
			for ( i=0; i<population; i++ ) {
				if (age[i]>1 && age[i]<22) {
					intervene[i]=1;	
				} 
			}
		} else if (interventions == 1 && t >= tauI_onset && t <= tauI_onset+dt) {
			for ( i=0; i<population; i++ ) {
				if (age[i]>=15 && age[i]<22) {
					intervene[i]=1;	
				} 
			}
		}	

		/* Segment population into infectious, susceptible, hospitalized, and icu */
		segment_population( &num_sus,  &num_infectious,  &num_hosp,  &num_icu,  infected,  infectious,  sus_list,  hosp_list,  hosp_pop,  icu_pop,  icu_list,  tau, population, t) ;

		//#### Only Susceptible people can get the virus and infected people spread it.
		for (i=0; i<num_sus; i++) {
			sus_person=sus_list[i];
			community_nom=0;
			community_den=0;
			infect = 0;
			contact_commun=0;
			contact_work=0;
			contact_house=0;
			contact_school=0;

			for (j=0; j<num_infectious; j++) {


				infec_person=infectious[j];
				/* Determine if person is under individual interventions and set parameters */
				if (intervene[infec_person] == 1 && t>tau[infec_person]+tauI[interventions]) {
					Ih=interIh[interventions];
					Ic=interIc[interventions];
					Iw=interIw[interventions];
				} 
				

				/* This will probably have to move outside to a pair list.  NOTE: The list of coworkers/classmates and community members within contact may not completely overlap. i.e. a coworker could be outside of the realm of commumnity transmission if someone lives on the edge of a county. */	
				d=distance(lat[sus_person], lon[sus_person], lat[infec_person], lon[infec_person], 'K');
				kappa = calc_kappa( t,  tau[infec_person], symptomatic[infec_person]);

				
				if (hosp_pop[infec_person]==0) {
					if (HH[sus_person]==HH[infec_person]) {
						// Household transmission //
						infect+=Ih*calc_household_infect(kappa, omega, per_HH_size[HH[sus_person]], alpha, severe[infec_person]); 
						contact_house++;
					//	printf("house %f \n", infect);
					}

					// Workplace/School transmission: People must be in same workplace and job type. // 
					if ((workplace[sus_person]==workplace[infec_person]) && (job_status[sus_person]==job_status[infec_person]) && (job_status[sus_person]>0) && Iw[job_status[sus_person]]>0) {
						infect+=calc_workplace_infect(job_status[sus_person], kappa, omega, workplace_size[(job_status[sus_person])][(workplace[sus_person])], severe[infec_person], Iw) ;
						if (job_status[sus_person]<4) {
							contact_school++;
							if (class[sus_person]==class[infec_person]) {
								infect+=calc_workplace_infect(job_status[sus_person], kappa, omega, 15.00, severe[infec_person], Iw) ;
							}
						} else {
							contact_work++;
						}
			//		printf("work  %f %f %i %f \n", infect, Iw[job_status[sus_person]], job_status[sus_person], calc_workplace_infect(job_status[sus_person], kappa, omega, workplace_size[(job_status[sus_person])][(workplace[sus_person])], severe[infec_person], Iw));
					}

					// Community transmission // 
					if (d<100) {
						age_group=floor(age[sus_person]/5);
						community_nom+=Ic*calc_community_infect( age_group, kappa, omega, severe[infec_person], d, &community_den);
						contact_commun++;

			//			printf("commm %f \n", infect);
					}
				} else {
					/* In hospital, only have interaction with hospital workers and half interaction with family (household). */
					// Workplace/School transmission: People must be in same workplace and job type. // 
					if ((workplace[sus_person]==workplace_tmp[infec_person]) && (job_status[sus_person]==job_status[infec_person])) {
						infect+=calc_workplace_infect(job_status[sus_person], kappa, omega, workplace_size[(job_status[sus_person])][(workplace[sus_person])], severe[infec_person], Iw) ;
				//	printf("hosp  %f \n", infect);
					}
					// Household transmission //
					if (HH[sus_person]==HH[infec_person]) {
						infect+=0.25*calc_household_infect(kappa, omega, per_HH_size[HH[sus_person]], alpha, severe[infec_person]); 
				//	printf("family_hosp  %f\n ", infect);
					}

				}
			}
	
                        if (community_den > 0) {
                            //printf("commun  %f %f %f \n", community_nom, community_den, community_nom/community_den);
                            infect+=community_nom/community_den; // Community spread is additive nominator and denominator.  Must be outside of infectious persons loop.
                        }


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
				} 
				if (contact_commun>0) {
					num_contact_commun++;
				}
				if (contact_house>0) {
					num_contact_house++;
				}
				if (contact_school>0) {
					num_contact_school++;
				}
				num_infect++;
				/* Determine if following interventions only for interventions that effect individuals.*/
				if ( personinter == 0 ) {
					intervene[sus_person]=round(COV_rand()*complyI[interventions]);
					/* Intervention 4 is household quarantine. Applicable for whole household.  */
					if (interventions == 4 ) {
						for (i=0; i<population; i++) {
							if ( HH[sus_person] == HH[i]) {
								intervene[i]=round(COV_rand()*complyI[interventions]);
							}
						}
					}
				}
			}	

		}


		// After 5 days of symptoms, people are randomly put into hospital based on age distributions. //
		hosp_entry(t, num_infectious,  infectious,  age,  icu_pop,  hosp_pop,  symptomatic, tau, workplace_tmp, hosp_num, county, dt) ;
		
		// Release from hospital // 
		hosp_release(t, num_hosp,  hosp_list,  tau,  recovered, hosp_pop, &num_recovered, &recovered_hosp, &recovered_icu, dt) ;
		
		
		// 15 days after onset of symptoms, people randomly die based on age based percentages.  ICU patients are removed from ICU into regular hospital. // 
		
		num_dead = death(t, num_infectious,  infectious,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt) ;
		num_dead = death(t, num_icu,  icu_list,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt) ;
		num_dead = death(t, num_hosp,  hosp_list,  tau,  dead,  icu_pop,  hosp_pop,  symptomatic, num_dead, age, dt) ;


		// Recovered after 11 days (6 days of symptoms) if not in hospital/ICU. // 
		for (i=0; i<num_infectious; i++) {
			infec_person=infectious[i];
			if ((tau[infec_person]==t-11) && (hosp_pop[infec_person]==0) && (icu_pop[infec_person]==0)) {
				recovered[infec_person]=1;
				num_recovered++;
			}
		}

        ret = clock_gettime(CLOCK_MONOTONIC, &t2);
        step_time = ((double)t2.tv_sec + (double)t2.tv_nsec/nsdiv) - ((double)t1.tv_sec + (double)t1.tv_nsec/nsdiv);
	printf("Walltime/timestep %6.2f Timestep %6.2f num_infected %i num_infectious %i num_in_hosp %i num_in_icu %i num_dead %i recovered_tot %i recovered_from_hosp %i recovered_from_icu %i contact_work %i contact_school %i contact_home %i contact_community %i \n", step_time, t, num_infect, num_infectious, num_hosp, num_icu, num_dead, num_recovered, recovered_hosp, recovered_icu, num_contact_work, num_contact_school, num_contact_house, num_contact_commun);
	fflush(stdout);
	}
        ret = clock_gettime(CLOCK_MONOTONIC, &T2);
        step_time = ((double)T2.tv_sec + (double)T2.tv_nsec/nsdiv) - ((double)T1.tv_sec + (double)T1.tv_nsec/nsdiv);
        printf("Total time %8.3f\n", step_time);
	
	free(HH);
	free(per_HH_size);
	free(lat);
	free(lon);
	free(lat_city);
	free(long_city);
	free(city_size);
	free(county_size);
	free(city);
	free(city_county);
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
	free(infected_list);
	free(infectious);
	free(sus_list);

	free(hosp_pop);
	free(icu_pop);
	free(recovered);
	free(dead);
	free(workplace_tmp);

}
