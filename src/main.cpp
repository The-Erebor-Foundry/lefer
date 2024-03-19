// C Math Library
#include <math.h>

// C++ STD Libraries
#include <vector>


#include "lefer.hpp"

namespace lefer {


// Main APIs of the library ================================================================


/** Draw a curve in the flow field.
 *
 * This function draws a curve in the flow field by starting
 * in a specific point in the field (`x_start` and `y_start`), and
 * then, it starts to walk through the field, by followinf the direction
 * of the angles it encounters in the field.
 *
 * For more details check: https://pedro-faria.netlify.app/posts/2024/2024-02-19-flow-even/en/
 *
 *
 * @param curve_id the id of the curve you want to draw.
 * @param x_start the x coordinate of the starting point from which the function will start to draw your curve.
 * @param y_start the y coordinate of the starting point from which the function will start to draw your curve.
 * @param n_steps the number of steps used to draw your curve.
 * @param step_length the length/distance taken in each step.
* @param d_sep the "separation distance", i.e., the amount of distance that each curve must be from neighbouring curves.
* @param flow_field a `lefer::FlowField` that contains the 2D grid of angle values that define your flow field.
* @param density_grid the density grid to be used by the algorithm, i.e., a `lefer::DensityGrid` object.
*/
Curve draw_curve(int curve_id,
		 double x_start,
		 double y_start,
		 int n_steps,
		 double step_length,
		 double d_sep,
		 FlowField* flow_field,
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



/** Draws multiple evenly-spaced and non-overlapping curves in the flow field.
* 
* This function takes a starting point (`x_start` and `y_start`) in the flow field,
* and draws a initial curve in the flow field. After that, the function starts a loop process,
* to derivate `n_curves - 1` curves from this initial curve. All the curves that are drawn
* into the field are derived from this initial curve.
*
* In other words, it is not guaranteed that this function will draw exactly `n_curves` curves
* into the field, because, it might not have enough space for `n_curves` curves, considering your current settings. The function
* will attempt to draw as many curves as possible. As long as they are not overlapping
* each other, and they are not too close to neighbouring curves, the function will
* continue to draw curves into the field.
* 
* @param x_start the x coordinate of the starting point of the initial curve.
* @param y_start the y coordinate of the starting point of the initial curve.
* @param n_curves the number of curves that the function will attempt to draw from the initial curve.
* @param n_steps the number of steps that each curve drawn into the field will have.
* @param min_steps_allowed the minimum number of steps allowed for a curve. In other words, every curve that is drawn in the field must have at least `min_steps_allowed` steps.
* @param step_length the length (or distance) taken in each step (usually, you want to set this variable between 1% and 0.1% of the flow field width.
* @param d_sep the "separation distance", i.e., the amount of distance that each curve must be from neighbouring curves.
* @param flow_field a `lefer::FlowField` that contains the 2D grid of angle values that define your flow field.
* @param density_grid the density grid to be used by the algorithm, i.e., a `lefer::DensityGrid` object.
*/

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
	double x = x_start;
	double y = y_start;
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

	while (curve_id < n_curves && curves.size() < n_curves) {
		SeedPointsQueue queue = SeedPointsQueue(n_steps);
		if (curve_id >= curves.size()) {
			// There is no more curves to be analyzed in the queue
			break;
		}
		queue = collect_seedpoints(&curves.at(curve_id), d_sep);
		for (Point p: queue._points) {
			if (curves.size() >= n_curves) {
				break;
			}
			// check if it is valid given the current state
			if (density_grid->is_valid_next_step(p.x, p.y)) {
				// if it is, draw the curve from it
				Curve curve = draw_curve(
					curves.size(),
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
			}
		}

		curve_id++;
	}



	return curves;
}




/** Draws multiple non-overlapping curves in the flow field.
* 
* While `even_spaced_curves()` checks both the distance from the current curve to neighbouring curves,
* to ensure an even space between each curve, this function checks only if the current curve is
* overlapping or not other curves. In other words, you use this function if you care only to draw
* curves that do not overlap each other.
*
* This function takes a sequence of startings points (`starting_points`). For each starting point,
* this function will attempt to draw a curve from it. So, in this function, you have total
* control over which points exatly the curves starts from.
*
* It is not guaranteed that this function will draw exactly `n_curves` curves
* into the field, because, it might not have enough space for `n_curves` curves, considering your current settings. The function
* will attempt to draw as many curves as possible. As long as they are not overlapping
* each other, the function will
* continue to draw curves into the field.
* 
* @param starting_points a sequence of `lefer::Point` objects. Each `lefer::Point` object should describe a starting point for a single curve.
* @param n_steps the number of steps that each curve drawn into the field will have.
* @param min_steps_allowed the minimum number of steps allowed for a curve. In other words, every curve that is drawn in the field must have at least `min_steps_allowed` steps.
* @param step_length the length (or distance) taken in each step (usually, you want to set this variable between 1% and 0.1% of the flow field width.
* @param d_sep the "separation distance", i.e., the amount of distance that each curve must be from neighbouring curves.
* @param flow_field a `lefer::FlowField` that contains the 2D grid of angle values that define your flow field.
* @param density_grid the density grid to be used by the algorithm, i.e., a `lefer::DensityGrid` object.
*/


std::vector<Curve> non_overlapping_curves(std::vector<Point> starting_points,
					  int n_steps,
					  int min_steps_allowed,
					  double step_length,
					  double d_sep,
					  FlowField* flow_field,
					  DensityGrid* density_grid) {

	std::vector<Curve> curves;
	curves.reserve(starting_points.size());
	int curve_id = 0;
	for (Point start_point: starting_points) {
		double x_start = start_point.x;
		double y_start = start_point.y;
		// Check if this starting point is valid given the current state
		if (density_grid->is_valid_next_step(x_start, y_start)) {
			// if it is, draw the curve from it
			Curve curve = draw_curve(
				curve_id,
				x_start, y_start,
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
			curve_id++;
		}
	}


	return curves;
}






















// Utilitaries =======================================================

/** Calculate the distance between two points
* 
* @param x1 the x coordinate of point 1.
* @param y1 the y coordinate of point 1.
* @param x2 the x coordinate of point 2.
* @param y2 the y coordinate of point 2.
*/
double distance (double x1, double y1, double x2, double y2) {
	double s1 = pow(x2 - x1, 2.0);
	double s2 = pow(y2 - y1, 2.0);
	return sqrt(s1 + s2);
}

/** Transform a 2D index into a 1D index.
*
* @param x the x coordinate in a 2D grid.
* @param y the y coordinate in a 2D grid.
* @param grid_width the width of the 2D grid you are using.
*/
static int _grid_index_as_1d(int x, int y, int grid_width) {
	return x + grid_width * y;
}












// FlowField class =======================================================

/** The constructor for FlowField class.
*
* This constructor builds a wrapper object around a 2D grid of double values.
* i.e. a 2D array of double values.
* This 2D array of double values must be a heap-based (i.e. dinamically allocated) array.
*
* Very important, this grid must be a square, meaning that, the height and width
* of the field must be the same.
*
* @param flowfield the 2D array of double values that defines the flow field.
* @param field_width the width of the field.
*
*/
FlowField::FlowField(double** flow_field, int field_width) {
	_flow_field = flow_field;
	_field_width = field_width;
}


int FlowField::get_field_width() {
	return _field_width;
}


int FlowField::get_flow_field_col(double x) {
	return (int) x;
}

int FlowField::get_flow_field_row(double y) {
	return (int) y;
}

bool FlowField::off_boundaries(double x, double y) {
	return (
	x <= 0 ||
	y <= 0 ||
	x >= _field_width ||
	y >= _field_width
	);
}


double FlowField::get_angle(double x, double y) {
	int xi = get_flow_field_col(x);
	int yi = get_flow_field_row(y);
	return _flow_field[xi][yi];
}












// Curve class =======================================================


/** The constructor for a Curve object
 *
 * This constructor returns a empty Curve object, that you can use
 * to store the coordinates and information about a specific curve you
 * want to drawn into the field.
 *
 * @param id the id you want to give to this Curve object.
 * @param n_steps the number of steps you will use to draw this curve.
*/
Curve::Curve(int id, int n_steps) {
	_curve_id = id;
	_steps_taken = 0;
	_x.reserve(n_steps);
	_y.reserve(n_steps);
	_direction.reserve(n_steps);
	_step_id.reserve(n_steps);
}

void Curve::insert_step(double x_coord, double y_coord, int direction_id) {
	_x.emplace_back(x_coord);
	_y.emplace_back(y_coord);
	_direction.emplace_back(direction_id);
	_step_id.emplace_back(_steps_taken);
	_steps_taken++;
}









// DensityGrid class ============================================================================


/** The constructor for DensityGrid class
 *
 * The Jobard and Lefer algortihm works around a "supporting" 2D grid, called of "density grid".
 * Each cell (or coordinate) in this density grid is responsible for keeping tracking of
 * curves that are already drawn into that specific area of the flow field.
 *
 * In the `cell_capacity` argument, you can specify an amount of space you want to pre-allocate
 * for each cell in the density grid before the algorithm starts to run. This might improve
 * drawstically the performance of the algorithm, because you avoid the need for frequents
 * resizing and reallocation of the `std::vector` objects that represents each cell.
 *
 * In other words, if you pre-allocate enough space for each cell in the density grid,
 * then, the algorithm does not have to spend time resizing the cell every time it hits
 * the maximum capacity for that cell.
 *
 * @param flow_field_width the width of the flow field.
 * @param flow_field_height the height of the flow field.
* @param d_sep the "separation distance", i.e., the amount of distance that each curve must be from neighbouring curves.
* @param cell_capacity the capacity (or "space") you want to allocate for each cell in the density grid.
*/
DensityGrid::DensityGrid(int flow_field_width, int flow_field_height, double d_sep, int cell_capacity) {
	int grid_width = (int)(flow_field_width / d_sep);
	int grid_height = (int)(flow_field_height / d_sep);
	_d_sep = d_sep;
	_width = grid_width;
	_height = grid_height;
	_n_elements = grid_width * grid_height;
	_grid.reserve(_n_elements);

	for (int i = 0; i < _n_elements; i++) {
		DensityCell cell(cell_capacity);
		_grid.push_back(cell);
	}
}

int DensityGrid::get_density_col (double x) {
	double c = (x / _d_sep);
	return (int) c;
}

int DensityGrid::get_density_row (double y) {
	double r = (y / _d_sep);
	return (int) r;
}

int DensityGrid::get_density_index (double x, double y) {
	int col = get_density_col(x);
	int row = get_density_row(y);
	return col + _width * row;
}

int DensityGrid::get_density_index (int col, int row) {
	return col + _width * row;
}

bool DensityGrid::off_boundaries(double x, double y) {
	int c = get_density_col(x);
	int r = get_density_row(y);
	return (
	c <= 0 ||
	r <= 0 ||
	c >= _width ||
	r >= _height
	);
}

void DensityGrid::insert_coord(double x, double y) {
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
	}
}

void DensityGrid::insert_curve_coords(Curve* curve) {
	int steps_taken = curve->_steps_taken;
	for (int i = 0; i < steps_taken; i++) {
		insert_coord(curve->_x.at(i), curve->_y.at(i));
	}
}

bool DensityGrid::is_valid_next_step(double x, double y) {
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









// DensityCell class =============================================================================
DensityCell::DensityCell(int cell_capacity) {
	x = std::vector<double>(cell_capacity, 0.0);
	y = std::vector<double>(cell_capacity, 0.0);
	capacity = cell_capacity;
	space_used = 0;
}



 




// SeedPointsQueue class =========================================================================

SeedPointsQueue::SeedPointsQueue(int n_steps) {
	_capacity = n_steps * 2;
	_space_used = 0;
	_points.reserve(n_steps * 2);
}

bool SeedPointsQueue::is_empty() {
	return _space_used == 0;
}

void SeedPointsQueue::insert_coord(double x, double y) {
	Point p = {x, y};
	_points.emplace_back(p);
	_space_used++;
}

void SeedPointsQueue::insert_point(Point p) {
	_points.emplace_back(p);
	_space_used++;
}



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






} // namespace lefer
