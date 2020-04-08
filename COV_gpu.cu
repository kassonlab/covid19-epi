#define deg2rad(deg) (deg * M_PI / 180)
#define rad2deg(rad) (rad * 180 / M_PI)

static double calc_community_infect(float kappa, float omega, int severe, double d, float betac_scale) {

	/* need to work on this.  Perhaps we take a random distance for each two people based on population density, number of people in county, county area, etc. */
	float zeta[]={0.1, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.75, 0.50, 0.25, 0.25, 0.25} ; //   # Travel related parameter for community transmission. Ferguson Nature 2006
	double fd;
	float betac=0.103 ; // Scaled from betac=0.075 in influenza pandemic with R0=1.6, COVID-19 R0=2.2 (Ferguson 2020)

	fd=1/(1+pow((d/4), 3)); //kernel density function as parameterized for GB.
	return (betac_scale*betac*kappa*fd*(1+severe*(omega-1)));
}

static double distance(double lat1, double lon1, double lat2, double lon2, char unit) {
    double theta, dist;
    if ((lat1 == lat2) && (lon1 == lon2)) {
      return 0;
    }
    else {
      theta = lon1 - lon2;
  //    dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
      double ang1,ang2;
      ang1 = deg2rad(lat1);
      ang2 = deg2rad(lat2);
      dist = cos(ang1) * cos(ang2) * ( 1.0 + cos(deg2rad(theta)) ) - cos(ang1 + ang2);
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

static float calc_kappa(float t, float tau, int symptomatic, float dt, float * kappa_vals, int hosp, int icu, int full_kappa, float R0_scale) {

	float kappa;
	float t1;
	int t2;
	//###Determine kappa for infected person.  This is the infectiousness of the person based on time since infection started.  Latency period is 4.6 days.  Infection starts at 5.1 days and lasts for 6 days.  Sympotmatic people are twice as likely to infect others as asymptomatic.
	// Kappa is a log normal function with mean of -0.72 and standard deviation of 1.8.  From Ferguson Nature 2005
	if (t-tau <= 4.6) {
		kappa=0.;
	} else if (t-tau>11.1 && hosp==0 && icu==0) {
		kappa=0.; //# Recovered or dead
	} else {
		/* First 2 lines calculates kappa on the fly, second two get precalculated kappa from array. */
                if (full_kappa) {
                    t1=(log(t-tau-4.6)+0.72)/1.8;
                    kappa=exp(-0.5*pow(t1,2))/((t-tau-4.6)*1.8*sqrt(2*M_PI));
                } else {
                    t2=(t-tau)/dt;
                    kappa=kappa_vals[t2];
                }
	}
	if (symptomatic==0) {
		kappa=kappa*0.5;
	}
	return(kappa*R0_scale);
}

extern "C" void locale_infectious_step(int j, int num_infectious, int* infectious, float Ic, int* intervene, float t, float* tau, float* tauI, float* interIc, int* symptomatic, float dt, float* kappa_vals, int* hosp_pop, int* icu_pop, float* lat_locale, float* lon_locale, int* locale_HH, int* HH, struct locale* locale_list, float omega, int* severe, double* out_tmp_comm_inf, int full_kappa, float R0_scale, float betac_scale) {
	double tmp_comm_inf = 0.0;
	int i;
	for (i=0; i<num_infectious; i++) {
		int infec_person; //Counter for infected person.
		float kappa; // #Infectiousness
		float tIc;
		infec_person = infectious[i];
		tIc = Ic;
		if ( intervene[infec_person] > 0 && t>tau[infec_person]+tauI[intervene[infec_person]]) {
			tIc = interIc[intervene[infec_person]];
		}
		kappa = calc_kappa( t,  tau[infec_person], symptomatic[infec_person], dt, kappa_vals, hosp_pop[infec_person], icu_pop[infec_person], full_kappa, R0_scale);
	
		if (hosp_pop[infec_person]==0) {
			float d; //distance between people.
			// Community transmission //
#if !defined(USE_LOCALE_DISTANCE)
			d = distance(lat_locale[j], lon_locale[j], lat_locale[locale_HH[HH[infec_person]]], lon_locale[locale_HH[HH[infec_person]]], 'K');
#else
			d = locale_distance(locale_list[j], locale_list[locale_HH[HH[infec_person]]]);
#endif
			tmp_comm_inf += tIc*calc_community_infect( kappa, omega, severe[infec_person], d, betac_scale);
		}
	}

	*out_tmp_comm_inf = tmp_comm_inf;
}