#ifndef __POINT_DEFECTS_H
#define __POINT_DEFECTS_H

#include "../StandardIncludes.h"
#include "../Settings.h"
#include "../context/Context.h"

class Simulation;

struct SimulationParameters;

struct PointDefect
{
	/// The next defect in the same [111] crystal row.
	PointDefect* next;
	/// The lattice position of the defect along the atomic row.
	double position;
	/// The type of defect.
	int type;
	/// Position offset along the dislocation.
	double x;/// x in lattice coordinate;
	double y;/// y in lattice coordinate;
	/// Position offset number.
	int p_type;
	/// ****** initial position in real space[modify for moving solute
	Point3 p0;
	/// ****** updated position]
	Point3 p1;
};

/**
 * The point defect cloud.
 */
class PointDefects : public UsingContext
{
public:

	/// Constructor.
	PointDefects(Context& context, Simulation* sim, SimulationParameters& params);

	/// Returns a reference to the simulation's global parameter set.
	const SimulationParameters& params() const { return _params; }

	SimulationParameters& params() { return _params; }

	/// Returns the global simulation object.
	Simulation& simulation() const { return *_sim; }

	/// Dumps the current point defect cloud to an output file.
	void dumpToVTK(std::ostream& stream, double simulationTime);

	/// Returns the list of point defects in the model.
	const std::map<std::pair<double,double>,PointDefect*>& defects() const { return _defects; }

	/// Returns the list of point defects in the model.
	std::map<std::pair<double,double>, PointDefect*>& defects() { return _defects; }

	/// Returns the world-space position of a point defect.
	Point3 getWorldPosition(double x, double y, double z) const;

	bool isSoluteOnTheDislocation(const Point3& p1, const Point3& p2);

	/********************************** Modification functions *********************************/

	/// Initializes the point defect cloud.
	void initialize();

	/// Adjusts the window around the dislocation within which point defects are explicitly modeled.
	void updateWindow();

	bool migration(PointDefect* p, double x, double y, double z, double x_line, double y_line, double x1, double y1, double z1,double x1_line, double y1_line, int t);

	double bindPointDefects(const Point3& p0,const Point3& p1,double& deltaZ);

	bool isSoluteOnTheScrewDislocation(const Point3& p1, const Point3& p2);


private:

	/// Reference to the simulation's global parameters.
	SimulationParameters& _params;

	/// Pointer to the global simulation object.
	Simulation* _sim;

	/// The point defects in the model.
	///std::map<std::pair<int,int>, PointDefect*> _defects;
	std::map<std::pair<double,double>, PointDefect*> _defects;

	/// The memory pool for allocating point defect instances.
	boost::object_pool<PointDefect> _defectPool;

	/// Measures time spent computing stresses.
	//boost::timer::cpu_times _stressTime;

	/// The offset of the atomic crystal with respect to the lattice on which the dislocation moves.
	double _latticeOffsetX, _latticeOffsetY;


};

#endif // __POINT_DEFECTS_H
