#if 0
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static struct timespec ts_diff(struct timespec* ts1, struct timespec* ts2) {
	struct timespec ret = {
		.tv_sec = ts2->tv_sec - ts1->tv_sec,
		.tv_nsec = ts2->tv_nsec - ts1->tv_nsec,
	};
	if (ret.tv_nsec < 0) {
		ret.tv_sec -= 1;
		ret.tv_nsec += 1000000000L;
	}

	return ret;
}

static double const pi = 3.14159265358979323846;

double distance(double lat1, double lon1, double lat2, double lon2) {
    if (lat1 == lat2 && lon1 == lon2) {
        return 0.0;
    }
    double theta = (lon1 - lon2) * pi / 180.0;
    double lat1r = lat1 * pi / 180.0;
    double lat2r = lat2 * pi / 180.0;
    double dist = sin(lat1r) * sin(lat2r) + cos(lat1r) * cos(lat2r) * cos(theta);
    dist = acos(dist);

    // Angle converted from radians to centigrades is close to the distance in kilometers
    return dist * 100.0 * 200.0 / pi;
}

double latitude_offset(double km) {
	return acos(cos(km * pi / (20000.0))) * 180.0 / pi;
}

float lerpf(float a, float b, float f) {
	return a * (1.0 - f) + b * f;
}

float unlerpf(float a, float b, float x) {
	return (x - a) / (b - a);
}

float clampf(float x, float lo, float hi) {
	if (x < lo) {
		return lo;
	}
	if (x > hi) {
		return hi;
	}
	return x;
}

struct band_member {
	uint32_t id;
	float lat;
	float lon;
};

struct band {
	size_t capacity;
	size_t count;
	struct band_member* members;
};

int order_bm_by_lon_id_lat(void const* first, void const* second) {
	struct band_member const* a = first;
	struct band_member const* b = second;
	if (a->lon != b->lon) {
		return (a->lon < b->lon) ? -1 : 1;
	}
	if (a->id != b->id) {
		return (a->id < b->id) ? -1 : 1;
	}
	return (a->lat < b->lat) ? -1 : 1;
}

size_t index_bands(float lat, float lat_min, float lat_max, size_t band_count) {
	float band_mapping = unlerpf(lat_min, lat_max, lat) * band_count;
	size_t band_id = (size_t)clampf(band_mapping, 0.0, band_count - 1.0);
	return band_id;
}

int main(int argc, char** argv) {
	if (argc != 2) {
		fprintf(stderr, "usage: %s POPULATION.bin\n", strlen(argv[0]) ? argv[0] : "spatial_bench");
		return 1;
	}

	FILE* f = fopen(argv[1], "rb");
	if (!f) {
		perror(argv[1]);
		return 1;
	}

	struct timespec ts_pre_load;
	clock_gettime(CLOCK_MONOTONIC, &ts_pre_load);

	size_t coord_count;
	{
		uint32_t n;
		fread(&n, sizeof(uint32_t), 1, f);
		coord_count = n;
	}
	float* lat = malloc(coord_count * sizeof(float));
	fread(lat, sizeof(float), coord_count, f);
	float* lon = malloc(coord_count * sizeof(float));
	fread(lon, sizeof(float), coord_count, f);
	
	fclose(f);

	printf("total count: %zu\n", coord_count);

	size_t band_count = 200;
	float lat_min = 54.0;
	float lat_max = 70.0;

	float band_width = (lat_max - lat_min) / (float)band_count;
	double band_length = distance(0.0, 0.0, band_width, 0.0);
	printf("band-width: %f° (%f km (%f°)))\n", band_width, band_length, latitude_offset(band_length));

	// Populate bands
	struct band* bands = calloc(band_count, sizeof(struct band));
	for (size_t i = 0; i < coord_count; ++i) {
		size_t band_id = index_bands(lat[i], lat_min, lat_max, band_count);
		struct band* band = bands + band_id;
		if (band->count == band->capacity) {
			band->capacity = (size_t)ceil((band->capacity + 1) * 1.6);
			band->members = realloc(band->members, band->capacity * sizeof(struct band_member));
		}
		struct band_member* bm = band->members + band->count;
		bm->id = i;
		bm->lat = lat[i];
		bm->lon = lon[i];
		++band->count;
	}

	for (size_t i = 0; i < band_count; ++i) {
		struct band* band = bands + i;
		qsort(band->members, band->count, sizeof(struct band_member), order_bm_by_lon_id_lat);
	}

	// Sanity check
	size_t* band_member_count = calloc(band_count, sizeof(size_t));
	for (size_t i = 0; i < coord_count; ++i) {
		float band_mapping = unlerpf(lat_min, lat_max, lat[i]) * band_count;
		if (band_mapping >= 0.0 && band_mapping < band_count) {
			band_member_count[(size_t)band_mapping]++;
		}
		else {
			printf("%f (%f, %f) out of band\n", band_mapping, lat[i], lon[i]);
		}
	}

	for (size_t i = 0; i < band_count; ++i) {
		if (bands[i].count != band_member_count[i]) {
			printf("Band %zu (%f to %f): %zu != %zu\n",
				i,
				lerpf(lat_min, lat_max, i / (double)band_count),
				lerpf(lat_min, lat_max, (i+1) / (double)band_count),
				band_member_count[i],
				bands[i].count);
		}
	}

	struct timespec ts_post_load;
	clock_gettime(CLOCK_MONOTONIC, &ts_post_load);

	{
		struct timespec tsd = ts_diff(&ts_pre_load, &ts_post_load);
		printf("[init] time: %ld.%09ld\n", tsd.tv_sec, tsd.tv_nsec);
	}


	struct timespec ts_pre_banded;
	clock_gettime(CLOCK_MONOTONIC, &ts_pre_banded);
	size_t considered_banded = 0;
	size_t interactions_banded = 0;
	double total_dist_banded = 0.0;
	double lat_40km = latitude_offset(40.0);
	for (size_t i1 = 0; i1 < coord_count; ++i1) {
		size_t first_band = index_bands(lat[i1] - lat_40km, lat_min, lat_max, band_count);
		size_t last_band = index_bands(lat[i1] + lat_40km, lat_min, lat_max, band_count);
		// printf("%zu %zu\n", first_band, last_band);
		for (size_t b = first_band; b <= last_band; ++b) {
			struct band* band = bands + b;
			for (size_t i2 = 0; i2 < band->count; ++i2) {
				considered_banded += 1;
				double dist = distance(lat[i1], lon[i1], band->members[i2].lat, band->members[i2].lon);
				if (dist < 40.0) {
					interactions_banded += 1;
					total_dist_banded += dist;
				}
			}
		}
	}
	struct timespec ts_post_banded;
	clock_gettime(CLOCK_MONOTONIC, &ts_post_banded);

	{
		struct timespec tsd = ts_diff(&ts_pre_banded, &ts_post_banded);
		printf("[banded] Total distance: %f, interactions: %zu, considered: %zu, hit-rate: %f%%, time %ld.%09ld\n",
			total_dist_banded, interactions_banded, considered_banded,
			100.0 * interactions_banded / (double)considered_banded,
			tsd.tv_sec, tsd.tv_nsec);
	}

	struct timespec ts_pre_naive;
	clock_gettime(CLOCK_MONOTONIC, &ts_pre_naive);
	size_t considered = 0;
	size_t interactions = 0;
	double total_dist = 0.0;
#if 1
	for (size_t i1 = 0; i1 < coord_count; ++i1) {
		for (size_t i2 = 0; i2 < coord_count; ++i2) {
			considered += 1;
			double dist = distance(lat[i1], lon[i1], lat[i2], lon[i2]);
			if (dist < 40.0) {
				interactions += 1;
				total_dist += dist;
			}
		}
	}
#endif
	struct timespec ts_post_naive;
	clock_gettime(CLOCK_MONOTONIC, &ts_post_naive);
	
	{
		struct timespec tsd = ts_diff(&ts_pre_load, &ts_post_naive);
		printf("[naive] Total distance: %f, interactions %zu, considered: %zu, hit-rate: %f%%, time: %ld.%09ld\n",
			total_dist, interactions, considered,
			100.0 * interactions / (double)considered,
			tsd.tv_sec, tsd.tv_nsec);
	}

	for (size_t i = 0; i < band_count; ++i) {
		free(bands[i].members);
	}
	free(bands);

	free(lat);
	free(lon);
	return 0;
}
#else
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "COV_geo.h"

typedef struct timepoint {
	struct timespec ts;
} t_timepoint;

t_timepoint tp_now() {
	t_timepoint ret;
	clock_gettime(CLOCK_MONOTONIC, &ret.ts);
	return ret;
}

t_timepoint tp_diff(t_timepoint a, t_timepoint b) {
	t_timepoint ret = {
		.ts = {
			.tv_sec = b.ts.tv_sec - a.ts.tv_sec,
			.tv_nsec = b.ts.tv_nsec - a.ts.tv_nsec,
		}
	};
	if (ret.ts.tv_nsec < 0) {
		ret.ts.tv_sec -= 1;
		ret.ts.tv_nsec += 1000000000L;
	}

	return ret;
}

double tp_duration(t_timepoint tp) {
	return tp.ts.tv_sec + tp.ts.tv_nsec / 1000000000.0;
}

int main(int argc, char** argv) {
	if (argc != 2) {
		fprintf(stderr, "usage: %s POPULATION.bin\n", strlen(argv[0]) ? argv[0] : "spatial_bench");
		return 1;
	}

	FILE* f = fopen(argv[1], "rb");
	if (!f) {
		perror(argv[1]);
		return 1;
	}
	
	size_t coord_count;
	{
		uint32_t n;
		fread(&n, sizeof(uint32_t), 1, f);
		coord_count = n;
	}
	float* lat = malloc(coord_count * sizeof(float));
	fread(lat, sizeof(float), coord_count, f);
	float* lon = malloc(coord_count * sizeof(float));
	fread(lon, sizeof(float), coord_count, f);
	
	fclose(f);

	printf("total count: %zu\n", coord_count);

	size_t band_count = 200;
	size_t slice_count = 50;
	float lat_min = 54.0;
	float lat_max = 70.0;
	float lon_min = 10.3;
	float lon_max = 24.6;

	printf("geo params: %zu latitude bands, %zu longitude slices, "
		   "lat = [%0.3f° - %0.3f°], lon = [%0.3f° - %0.3f°]\n",
		band_count, slice_count, lat_min, lat_max, lon_min, lon_max);

	// Init phase
	t_timepoint tp_pre_init = tp_now();

	struct COV_geo* geo = COV_create_geo(coord_count, NULL, lat, lon,
		band_count, slice_count,
		lat_min, lon_min, lat_max, lon_max);

	t_timepoint tp_post_init = tp_now();
	printf("[init] time: %lf\n", tp_duration(tp_diff(tp_pre_init, tp_post_init)));

	float cutoff_distance = 40.0;
	{
		// Query phase (sliced)
		t_timepoint tp_pre_sliced = tp_now();

		size_t considered_banded = 0;
		size_t interactions_banded = 0;
		double total_dist_banded = 0.0;
		for (size_t i1 = 0; i1 < coord_count; ++i1) {
			struct COV_geo_query_context qc = COV_geo_prepare_query(geo, lat[i1], lon[i1], cutoff_distance);
			struct COV_geo_slice* slice;
			while ((slice = COV_geo_query_next(&qc))) {
				// printf("slice %p has %zu elements\n", (void*)slice, slice->count);
				for (size_t si = 0; si < slice->count; ++si) {
					considered_banded += 1;
					double dist = COV_geo_distance(lat[i1], lon[i1], slice->lats[si], slice->lons[si]);
					if (dist < 40.0) {
						interactions_banded += 1;
						total_dist_banded += dist;
					}
				}
			}
		}
		
		t_timepoint tp_post_sliced = tp_now();
		printf("[sliced] Total distance: %f, interactions: %zu, considered: %zu, hit-rate: %f%%, time %lf\n",
			total_dist_banded, interactions_banded, considered_banded,
			100.0 * interactions_banded / (double)considered_banded,
			tp_duration(tp_diff(tp_pre_sliced, tp_post_sliced)));
	}

	return 0;
}
#endif