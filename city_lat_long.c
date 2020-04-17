#include <stdio.h>
#include <string.h>
#include "city_lat_long.h"

void city_lat_long(int *num_cities, double * lat_city, double * long_city, char ** cities, int * county, char ** county_names, int num_county) {

  /* Get longitude and latitude of cities in Sweden from CSV */
  FILE* fp = fopen("cities_all.csv", "r");  // Not sure about the validity of this file.  Could use a better source.

  /* Counters for parsing file. */
  char tmp[199]; // municipality name
  char tmp1[199]; // locality name
  char tmp2[199]; // county name
  double tmp_lat; // tmp latitude
  double tmp_lon; // tmp longitude
  int i, j;

  for (i=0; i < 2000; i++) {
    int got = fscanf(fp, "%[^,\n]%*c%[^,\n]%*c%[^,\n]%*c%lf%*c%lf\n", tmp1, tmp, tmp2, &tmp_lat, &tmp_lon);
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
