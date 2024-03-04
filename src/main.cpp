// C Math Library
#include <math.h>

// C++ STD Libraries
#include <vector>
#include <iostream>


namespace lefer {

double distance (double x1, double y1, double x2, double y2) {
	double s1 = pow(x2 - x1, 2.0);
	double s2 = pow(y2 - y1, 2.0);
	return sqrt(s1 + s2);
}


static int _grid_index_as_1d(int x, int y, int grid_width) {
	return x + grid_width * y;
}


static int _get_flow_field_col(double x) {
	return (int) x;
}

static int _get_flow_field_row(double y) {
	return (int) y;
}

static bool _off_boundaries(double x, double y, int limit) {
	return (
	x <= 0 ||
	y <= 0 ||
	x >= limit ||
	y >= limit 
	);
}

static double get_angle(double** _flow_field, double x, double y) {
	int xi = _get_flow_field_col(x);
	int yi = _get_flow_field_row(y);
	return _flow_field[xi][yi];
}


struct Point {
	double x;
	double y;
};


class Curve {
public:
	int _curve_id;
	std::vector<double> _x;
	std::vector<double> _y;
	std::vector<int> _direction;
	std::vector<int> _step_id;
	int _steps_taken;
public:
	Curve(int id, int n_steps) {
		_curve_id = id;
		_steps_taken = 0;
		_x.reserve(n_steps);
		_y.reserve(n_steps);
		_direction.reserve(n_steps);
		_step_id.reserve(n_steps);
	}

	void insert_step(double x_coord, double y_coord, int direction_id) {
		_x.emplace_back(x_coord);
		_y.emplace_back(y_coord);
		_direction.emplace_back(direction_id);
		_step_id.emplace_back(_steps_taken);
		_steps_taken++;
	}
};


struct DensityCell {
	std::vector<double> x;
	std::vector<double> y;
	int capacity;
	int space_used;
};


class DensityGrid {
private:
	std::vector<DensityCell> _grid;
	int _width;
	int _height;
	int _n_elements;
	double _d_sep;
public:
	DensityGrid(int grid_width, int grid_height, double d_sep, int cell_capacity) {
		_d_sep = d_sep;
		_width = grid_width;
		_height = grid_height;
		_n_elements = grid_width * grid_height;
		_grid.reserve(grid_width * grid_height);

		for (int i = 0; i < _n_elements; i++) {
			_grid[i].x.reserve(cell_capacity);
			_grid[i].y.reserve(cell_capacity);
			_grid[i].capacity = cell_capacity;
			_grid[i].space_used = 0;
		}
	}

	int get_density_col (double x) {
		double c = (x / _d_sep);
		return (int) c;
	}

	int get_density_row (double y) {
		double r = (y / _d_sep);
		return (int) r;
	}

	int get_density_index (double x, double y) {
		int col = get_density_col(x);
		int row = get_density_row(y);
		return col + _width * row;
	}
	
	int get_density_index (int col, int row) {
		return col + _width * row;
	}

	bool off_boundaries(double x, double y) {
		int c = get_density_col(x);
		int r = get_density_row(y);
		return (
		c <= 0 ||
		r <= 0 ||
		c >= _width ||
		r >= _height
		);
	}

	void insert_coord(double x, double y) {
		if (off_boundaries(x, y)) {
			return;
		}

		int density_index = get_density_index(x, y);
		int space_used = _grid[density_index].space_used;
		int capacity = _grid[density_index].capacity;

		if ((space_used + 1) < capacity) {
			_grid[density_index].x.emplace_back(x);
			_grid[density_index].y.emplace_back(y);
			_grid[density_index].space_used++;
		} else {
			std::cout << "[ERROR]: Attempt to add coordinate in density cell that is out of capacity!" << std::endl;
		}
	}

	void insert_curve_coords(Curve* curve) {
		int steps_taken = curve->_steps_taken;
		for (int i = 0; i < steps_taken; i++) {
			insert_coord(curve->_x.at(i), curve->_y.at(i));
		}
	}

	bool is_valid_next_step(double x, double y) {
		if (off_boundaries(x, y)) {
			return 0;
		}

		int density_col = get_density_col(x);
		int density_row = get_density_row(y);
		int start_row = (density_row - 1) > 0 ? density_row - 1 : 0;
		int end_row = (density_row + 1) < _width ? density_row + 1 : density_row; 
		int start_col = (density_col - 1) > 0 ? density_col - 1 : 0;
		int end_col = (density_col + 1) < _height ? density_col + 1 : density_col;

		// Subtracting a very small amount from D_TEST, just to account for the lost of float precision
		// that happens during the calculations below, specially in the distance calc
		double d_test = _d_sep - (0.01 * _d_sep);
		for (int c = start_col; c <= end_col; c++) {
			for (int r = start_row; r <= end_row; r++) {
				int density_index = get_density_index(c, r);
				int n_elements = _grid[density_index].space_used;
				if (n_elements == 0) {
					continue;
				}

				for (int i = 0; i < n_elements; i++) {
					double x2 = _grid[density_index].x.at(i);
					double y2 = _grid[density_index].y.at(i);
					double dist = distance(x, y, x2, y2);
					if (dist <= d_test) {
						return 0;
					}
				}
			}
		}

		return 1;
	}
};


class SeedPointsQueue {
public:
	std::vector<Point> _points;
	int _capacity;
	int _space_used;

public:
	SeedPointsQueue(int n_steps) {
		_capacity = n_steps * 2;
		_space_used = 0;
		_points.reserve(n_steps * 2);
	}

	bool is_empty() {
		return _space_used == 0;
	}

	void insert_coord(double x, double y) {
		Point p = {x, y};
		_points.emplace_back(p);
		_space_used++;
	}

	void insert_point(Point p) {
		_points.emplace_back(p);
		_space_used++;
	}
};



SeedPointsQueue collect_seedpoints (Curve* curve, double d_sep) {
	int steps_taken = curve->_steps_taken;
	SeedPointsQueue queue = SeedPointsQueue(steps_taken);
	if (steps_taken == 0) {
		return queue;
	}

	for (int i = 0; i < steps_taken - 1; i++) {
		double x = curve->_x.at(i);
		double y = curve->_y.at(i);

		int ff_column_index = (int) floor(x);
		int ff_row_index = (int) floor(y);
		double angle = atan2(curve->_y.at(i + 1) - y, curve->_x.at(i + 1) - x);

		double angle_left = angle + (M_PI / 2);
		double angle_right = angle - (M_PI / 2);

		Point left_point = {
			x + (d_sep * cos(angle_left)),
			y + (d_sep * sin(angle_left))
		};
		Point right_point = {
			x + (d_sep * cos(angle_right)),
			y + (d_sep * sin(angle_right))
		};

		queue.insert_point(left_point);	
		queue.insert_point(right_point);	
	}

	return queue;
}




Curve draw_curve(int curve_id,
		 double x_start,
		 double y_start,
		 int n_steps,
		 double step_length,
		 double d_sep,
		 double** flow_field,
		 DensityGrid* density_grid) {

	Curve curve = Curve(curve_id, n_steps);
	curve.insert_step(x_start, y_start, 0);
	double x = x_start;
	double y = y_start;
	int i = 1;
	// Draw curve from right to left
	while (i < (n_steps / 2)) {
		if (_off_boundaries(x, y, flow_field_width)) {
			break;
		}

		double angle = flow_field->get_angle(x, y);
		double x_step = step_length * cos(angle);
		double y_step = step_length * sin(angle);
		x = x - x_step;
		y = y - y_step;

		if (!density_grid->is_valid_next_step(x, y)) {
			break;
		}

		curve.insert_step(x, y, 0);
		i++;
	}

	x = x_start;
	y = y_start;
	// Draw curve from left to right
	while (i < n_steps) {
		if (_off_boundaries(flow_field, x, y)) {
			break;
		}

		double angle = flow_field->get_angle(x, y);
		double x_step = step_length * cos(angle);
		double y_step = step_length * sin(angle);
		x = x + x_step;
		y = y + y_step;

		if (!density_grid->is_valid_next_step(x, y)) {
			break;
		}

		curve.insert_step(x, y, 1);
		i++;
	}

	return curve;
}


std::vector<Curve> even_spaced_curves(double x_start,
				      double y_start,
				      int n_curves,
				      int n_steps,
				      int min_steps_allowed,
				      double step_length,
				      double d_sep,
				      double** flow_field,
				      DensityGrid* density_grid) {

	std::vector<Curve> curves;
	curves.reserve(n_curves);
	double x = x_start;
	double y = y_start;
	int curve_array_index = 0;
	int curve_id = 0;
	Curve curve = draw_curve(
		curve_id,
		x, y,
		n_steps,
		step_length,
		d_sep,
		flow_field,
		density_grid
	);

	curves.emplace_back(curve);
	density_grid->insert_curve_coords(&curve);
	curve_array_index++;


	while (curve_id < n_curves && curve_array_index < n_curves) {
		SeedPointsQueue queue = SeedPointsQueue(n_steps);
		if (curve_id >= curves.size()) {
			// There is no more curves to be analyzed in the queue
			break;
		}
		queue = collect_seedpoints(&curves.at(curve_id), d_sep);
		for (Point p: queue._points) {
			// check if it is valid given the current state
			if (density_grid->is_valid_next_step(p.x, p.y)) {
				// if it is, draw the curve from it
				Curve curve = draw_curve(
					curve_array_index,
					p.x, p.y,
					n_steps,
					step_length,
					d_sep,
					flow_field,
					density_grid
				);

				if (curve._steps_taken < min_steps_allowed) {
					continue;
				}

				curves.emplace_back(curve);
				// insert this new curve into the density grid
				density_grid->insert_curve_coords(&curve);
				curve_array_index++;
			}
		}

		curve_id++;
	}



	return curves;
}




} // namespace lefer
