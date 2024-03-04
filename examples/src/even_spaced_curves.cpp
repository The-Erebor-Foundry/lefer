#include <iostream>
#include <vector>

#include "lefer.hpp"

#define FNL_IMPL
#include "./../FastNoiseLite.h"


int main (int argc, char *argv[]) {

	int flow_field_width = 120;
	int flow_field_height = 120;
	int n_steps = 30;
	int min_steps_allowed = 5;
	double step_length = 0.01 * flow_field_width;
	double d_sep = 0.8;
	int n_curves = 1500;

	double** flow_field;
	flow_field = (double**)malloc(sizeof(double*) * flow_field_width);
	for (int i = 0; i< flow_field_width; i++) {
		flow_field[i] = (double*)malloc(sizeof(double) * flow_field_height);
	}

	// Create and configure noise state
	fnl_state noise = fnlCreateState();
	noise.seed = 50;
	noise.noise_type = FNL_NOISE_PERLIN;
	for (int y = 0; y < flow_field_height; y++) {
		for (int x = 0; x < flow_field_width; x++) {
			flow_field[x][y] = fnlGetNoise2D(&noise, x, y) * 2 * M_PI;
		}
	}


	lefer::FlowField flow_field_obj = lefer::FlowField(flow_field, flow_field_width);
	lefer::DensityGrid density_grid = lefer::DensityGrid(flow_field_width, flow_field_height, d_sep, 2000);
	
	double x_start = 45.0;
	double y_start = 24.0;
	std::vector<lefer::Curve> curves = lefer::even_spaced_curves(
		x_start,
		y_start,
		n_curves,
		n_steps,
		min_steps_allowed,
		step_length,
		d_sep,
		&flow_field_obj,
		&density_grid
	);


	// Visualizing the coordinates calculated by the algorithm
	for (lefer::Curve curve: curves) {
		int curve_id = curve._curve_id;
		for (int i = 0; i < curve._steps_taken; i++) {
			std::cout << curve_id << "; "
				<< curve._x[i] << "; "
				<< curve._y[i] << "; "
				<< curve._direction[i] << "; "
				<< std::endl;
		}
	}

	return 0;
}
