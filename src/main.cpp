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


static int _array_index_as_1d(int x, int y, int array_width) {
	return x + array_width * y;
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
	DensityCell _grid[DENSITY_GRID_WIDTH][DENSITY_GRID_HEIGHT];
	int _width;
	int _height;
	double _d_sep;
public:
	DensityGrid(double d_sep, int cell_capacity) {
		_d_sep = d_sep;
		_width = DENSITY_GRID_WIDTH;
		_height = DENSITY_GRID_HEIGHT;
		for (int r = 0; r < _height; r++) {
			for (int c = 0; c < _width; c++) {
				_grid[c][r].x.reserve(cell_capacity);
				_grid[c][r].y.reserve(cell_capacity);
				_grid[c][r].capacity = cell_capacity;
				_grid[c][r].space_used = 0;
			}
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

		int density_col = get_density_col(x);
		int density_row = get_density_row(y);

		int space_used = _grid[density_col][density_row].space_used;
		int capacity = _grid[density_col][density_row].capacity;

		if ((space_used + 1) < capacity) {
			_grid[density_col][density_row].x.emplace_back(x);
			_grid[density_col][density_row].y.emplace_back(y);
			_grid[density_col][density_row].space_used++;
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
				int n_elements = _grid[c][r].space_used;
				if (n_elements == 0) {
					continue;
				}

				for (int i = 0; i < n_elements; i++) {
					double x2 = _grid[c][r].x.at(i);
					double y2 = _grid[c][r].y.at(i);
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



SeedPointsQueue collect_seedpoints (Curve* curve) {
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
			x + (D_SEP * cos(angle_left)),
			y + (D_SEP * sin(angle_left))
		};
		Point right_point = {
			x + (D_SEP * cos(angle_right)),
			y + (D_SEP * sin(angle_right))
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
		 double** _flow_field,
		 DensityGrid* density_grid) {

	Curve curve = Curve(curve_id, n_steps);
	curve.insert_step(x_start, y_start, 0);
	double x = x_start;
	double y = y_start;
	int i = 1;
	// Draw curve from right to left
	while (i < (n_steps / 2)) {
		if (flow_field->off_boundaries(x, y)) {
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
		if (flow_field->off_boundaries(x, y)) {
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
				      FlowField* flow_field,
				      DensityGrid* density_grid) {

	std::vector<Curve> curves;
	curves.reserve(n_curves);
	/*for (int curve_id = 0; curve_id < n_curves; curve_id++) {
		curves[curve_id] = Curve(curve_id, n_steps);
	}*/

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
		queue = collect_seedpoints(&curves.at(curve_id));
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




int main() {
	FlowField flow_field = FlowField();
	DensityGrid density_grid = DensityGrid(D_SEP, 5000);
	std::vector<Curve> curves = even_spaced_curves(
		45.0, 24.0,
		N_CURVES,
		N_STEPS,
		MIN_STEPS_ALLOWED,
		STEP_LENGTH,
		D_SEP,
		&flow_field,
		&density_grid
	);

	return 1;
}



} // namespace lefer
