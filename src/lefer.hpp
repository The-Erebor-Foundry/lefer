#include <vector>


namespace lefer {

double distance (double x1, double y1, double x2, double y2);
static int _grid_index_as_1d(int x, int y, int grid_width);


class FlowField {
private:
	double** _flow_field;
	int _field_width;
public:
	FlowField(double** flow_field, int field_width);
	int get_field_width();
	int get_flow_field_col(double x);
	int get_flow_field_row(double y);
	bool off_boundaries(double x, double y);
	double get_angle(double x, double y);
};

struct Point {
	double x;
	double y;
};


/** A class that represents a curve
* 
*/
class Curve {
public:
	/** The id that identifies the curve */
	int _curve_id; 
	std::vector<double> _x; /** The x coordinates of each point in the curve */
	std::vector<double> _y; /** The y coordinates of each point in the curve */
	std::vector<int> _direction; /** The direction id of each point in the curve (0 means direction from left to right, 1 means direction from right to left) */
	std::vector<int> _step_id; /** The id (or the number) of the step for each point in the curve */
	int _steps_taken; /** The number of steps taken to draw the curve. */
public:
	Curve(int id, int n_steps);
	void insert_step(double x_coord, double y_coord, int direction_id);
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
	DensityGrid(int flow_field_width, int flow_field_height, double d_sep, int cell_capacity);
	int get_density_col (double x);
	int get_density_row (double y);
	int get_density_index (double x, double y);
	int get_density_index (int col, int row);
	bool off_boundaries(double x, double y);
	void insert_coord(double x, double y);
	void insert_curve_coords(Curve* curve);
	bool is_valid_next_step(double x, double y);
};



class SeedPointsQueue {
public:
	std::vector<Point> _points;
	int _capacity;
	int _space_used;

public:
	SeedPointsQueue(int n_steps);
	bool is_empty();
	void insert_coord(double x, double y);
	void insert_point(Point p);
};



SeedPointsQueue collect_seedpoints (Curve* curve, double d_sep);



Curve draw_curve(int curve_id,
		 double x_start,
		 double y_start,
		 int n_steps,
		 double step_length,
		 double d_sep,
		 FlowField* flow_field,
		 DensityGrid* density_grid);


std::vector<Curve> even_spaced_curves(double x_start,
				      double y_start,
				      int n_curves,
				      int n_steps,
				      int min_steps_allowed,
				      double step_length,
				      double d_sep,
				      FlowField* flow_field,
				      DensityGrid* density_grid);



std::vector<Curve> non_overlapping_curves(std::vector<Point> starting_points,
				      int n_steps,
				      int min_steps_allowed,
				      double step_length,
				      double d_sep,
				      FlowField* flow_field,
				      DensityGrid* density_grid);



} // namespace lefer
