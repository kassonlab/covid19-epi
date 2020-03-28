#include "COV_geo.h"

#include <math.h>
#include <stdbool.h>
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

double COV_geo_vertical_distance(double lat1, double lon, double lat2) {
	if (lat1 == lat2) {
        return 0.0;
    }
    double theta = 0.0;
    double lat1r = lat1 * pi / 180.0;
    double lat2r = lat2 * pi / 180.0;
    double dist = sin(lat1r) * sin(lat2r) + cos(lat1r) * cos(lat2r);
    dist = acos(dist);

    // Angle converted from radians to centigrades is close to the distance in kilometers
    return dist * 100.0 * 200.0 / pi;
}

double COV_geo_horizontal_distance(double lat, double lon1, double lon2) {
    if (lon1 == lon2) {
        return 0.0;
    }
    double theta = (lon1 - lon2) * pi / 180.0;
    double latr = lat * pi / 180.0;
	double sin_latr = sin(latr);
	double cos_latr = cos(latr);
    double dist = sin_latr * sin_latr + cos_latr * cos_latr * cos(theta);
    dist = acos(dist);

    // Angle converted from radians to centigrades is close to the distance in kilometers
    return dist * 100.0 * 200.0 / pi;
}

double COV_geo_distance(double lat1, double lon1, double lat2, double lon2) {
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

static double latitude_offset(double km) {
	return acos(cos(km * pi / (20000.0))) * 180.0 / pi;
}

static float lerpf(float a, float b, float f) {
	return a * (1.0 - f) + b * f;
}

static float unlerpf(float a, float b, float x) {
	return (x - a) / (b - a);
}

static float clampf(float x, float lo, float hi) {
	if (x < lo) {
		return lo;
	}
	if (x > hi) {
		return hi;
	}
	return x;
}

struct person {
	uint32_t id;
	float lat;
	float lon;
};

static int order_bm_by_lon_id_lat(void const* first, void const* second) {
	struct person const* a = first;
	struct person const* b = second;
	if (a->lon != b->lon) {
		return (a->lon < b->lon) ? -1 : 1;
	}
	if (a->id != b->id) {
		return (a->id < b->id) ? -1 : 1;
	}
	return (a->lat < b->lat) ? -1 : 1;
}

static size_t index_bands(float lat, float lat_min, float lat_max, size_t band_count) {
	float band_mapping = unlerpf(lat_min, lat_max, lat) * band_count;
	size_t band_id = (size_t)clampf(band_mapping, 0.0, band_count - 1.0);
	return band_id;
}

static size_t index_slices(float lon, float lon_min, float lon_max, size_t slice_count) {
	float slice_mapping = unlerpf(lon_min, lon_max, lon) * slice_count;
	size_t slice_id = (size_t)clampf(slice_mapping, 0.0, slice_count - 1.0);
	return slice_id;
}

static float band_lower_latitude(size_t band_id, float lat_min, float lat_max, size_t band_count) {
	return lerpf(lat_min, lat_max, band_id / (float)band_count);
}

static float slice_lower_longitude(size_t slice_id, float lon_min, float lon_max, size_t slice_count) {
	return lerpf(lon_min, lon_max, slice_id / (float)slice_count);
}

struct COV_geo* COV_create_geo(size_t coord_count, uint32_t* ids, float* lat, float* lon,
	size_t band_count, size_t slice_count,
	float lat_min, float lon_min, float lat_max, float lon_max)
{
	struct COV_geo* geo = calloc(sizeof(struct COV_geo), 1);
	geo->band_count = band_count;
	geo->lat_min = lat_min;
	geo->lat_max = lat_max;
	geo->lon_min = lon_min;
	geo->lon_max = lon_max;

	float band_width = (lat_max - lat_min) / (float)band_count;
	float slice_width = (lon_max - lon_min) / (float)slice_count;
	double band_length = COV_geo_distance(0.0, 0.0, band_width, 0.0);
	printf("band-width: %f° (%f km (%f°)))\n", band_width, band_length, latitude_offset(band_length));

	// Populate bands
	geo->bands = calloc(band_count, sizeof(struct COV_geo_band));
	for (size_t band_id = 0; band_id < band_count; ++band_id) {
		struct COV_geo_band* band = geo->bands + band_id;
		band->slice_count = slice_count;
		band->slices = calloc(slice_count, sizeof(struct COV_geo_slice));
	}
	for (size_t i = 0; i < coord_count; ++i) {
		size_t band_id = index_bands(lat[i], lat_min, lat_max, band_count);
		size_t slice_id = index_slices(lon[i], lon_min, lon_max, slice_count);
		struct COV_geo_band* band = geo->bands + band_id;
		struct COV_geo_slice* slice = band->slices + slice_id;

		if (slice->count == slice->capacity) {
			slice->capacity = (size_t)ceil((slice->capacity + 1) * 1.6);
			slice->ids = realloc(slice->ids, slice->capacity * sizeof(uint32_t));
			slice->lats = realloc(slice->lats, slice->capacity * sizeof(float));
			slice->lons = realloc(slice->lons, slice->capacity * sizeof(float));
		}
		slice->ids[slice->count] = ids ? ids[i] : i;
		slice->lats[slice->count] = lat[i];
		slice->lons[slice->count] = lon[i];
		++slice->count;
	}

#if 0
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

	// Print statistics
	for (size_t i = 0; i < band_count; ++i) {
		struct COV_geo_band* band = geo->bands + i;
		if (band->count != band_member_count[i]) {
			printf("Band %zu (%f to %f): %zu != %zu\n",
				i,
				lerpf(lat_min, lat_max, i / (double)band_count),
				lerpf(lat_min, lat_max, (i+1) / (double)band_count),
				band_member_count[i],
				band->count);
		}
	}
#endif

	return geo;
}

void COV_destroy_geo(struct COV_geo* geo) {
	for (size_t i = 0; i < geo->band_count; ++i) {
		for (size_t j = 0; j < geo->bands[i].slice_count; ++j) {
			free(geo->bands[i].slices[j].ids);
			free(geo->bands[i].slices[j].lats);
			free(geo->bands[i].slices[j].lons);
		}
		free(geo->bands[i].slices);
	}
	free(geo->bands);
	free(geo);
}

struct COV_geo_query_context COV_geo_prepare_query(struct COV_geo* geo, float lat, float lon, float cutoff_distance) {
	struct COV_geo_query_context qc;
	qc.geo = geo;
	qc.cutoff_distance = cutoff_distance;
	qc.lat_cutoff = latitude_offset(cutoff_distance);
	qc.first_band = index_bands(lat - qc.lat_cutoff, geo->lat_min, geo->lat_max, geo->band_count);
	qc.last_band = index_bands(lat + qc.lat_cutoff, geo->lat_min, geo->lat_max, geo->band_count);
	qc.mid_slice = index_slices(lon, geo->lon_min, geo->lon_max, geo->bands->slice_count);
	qc.iteration_started = false;
	return qc;
}

bool set_query_band(struct COV_geo_query_context* qc, size_t band_id) {
	qc->current_band = band_id;
	if (band_id < qc->first_band || band_id > qc->last_band) {
		return false;
	}
	struct COV_geo_band* band = qc->geo->bands + band_id;

	/* The idea for finding the westmost/eastmost slices is to find the
	latitude border of the band that has the shortest length. This allows
	us to compute the number of slices that we need to touch for a conservative
	estimate of the number of slices we need to consider for the current band.
	*/
	double band_lat_south = band_lower_latitude(band_id, qc->geo->lat_min, qc->geo->lat_max, qc->geo->band_count);
	double band_lat_north = band_lower_latitude(band_id + 1, qc->geo->lat_min, qc->geo->lat_max, qc->geo->band_count);
	double band_lat_south_mag = fabs(band_lat_south);
	double band_lat_north_mag = fabs(band_lat_north);

	// The narrowest edge of the slice is the one that has the latitude closest to a pole.
	double band_lat_narrowest_edge = (band_lat_south_mag > band_lat_north_mag)
		? band_lat_south_mag
		: band_lat_north_mag;

	double slice_lon_west = slice_lower_longitude(0, qc->geo->lon_min, qc->geo->lon_max, band->slice_count);
	double slice_lon_east = slice_lower_longitude(1, qc->geo->lon_min, qc->geo->lon_max, band->slice_count);

	double narrow_width = COV_geo_horizontal_distance(band_lat_narrowest_edge, slice_lon_west, slice_lon_east);
	size_t narrow_cells = ceilf(qc->cutoff_distance / narrow_width);
	size_t narrow_cells_east = narrow_cells + 1;
		
	size_t west_cells_clamped = narrow_cells > qc->mid_slice ? qc->mid_slice : narrow_cells;
	qc->current_slice_start = qc->mid_slice - west_cells_clamped;

	qc->current_slice_end = qc->mid_slice + narrow_cells_east;
	qc->current_slice_end = qc->current_slice_end > band->slice_count
		? band->slice_count
		: qc->current_slice_end;

	qc->current_slice = qc->current_slice_start;

	return true;
}

struct COV_geo_slice* COV_geo_query_next(struct COV_geo_query_context* qc) {
	if (!qc->iteration_started) {
		set_query_band(qc, qc->first_band);
		qc->iteration_started = true;
	}
	if (qc->current_band <= qc->last_band) {
		struct COV_geo_slice* ret = &qc->geo->bands[qc->current_band].slices[qc->current_slice];
		if (++qc->current_slice == qc->current_slice_end) {
			set_query_band(qc, qc->current_band + 1);
		}
		return ret;
	}

	return NULL;
}

static void query_banded(struct COV_geo* geo, size_t coord_count, float* lat, float* lon, float cutoff_distance) {
	size_t considered_banded = 0;
	size_t interactions_banded = 0;
	double total_dist_banded = 0.0;
	for (size_t i1 = 0; i1 < coord_count; ++i1) {
		struct COV_geo_query_context qc = COV_geo_prepare_query(geo, lat[i1], lon[i1], cutoff_distance);
		struct COV_geo_slice* slice;
		while ((slice = COV_geo_query_next(&qc))) {
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
}

static void query_naive(size_t coord_count, float* lat, float* lon) {
	size_t considered = 0;
	size_t interactions = 0;
	double total_dist = 0.0;
	for (size_t i1 = 0; i1 < coord_count; ++i1) {
		for (size_t i2 = 0; i2 < coord_count; ++i2) {
			considered += 1;
			double dist = COV_geo_distance(lat[i1], lon[i1], lat[i2], lon[i2]);
			if (dist < 40.0) {
				interactions += 1;
				total_dist += dist;
			}
		}
	}
}