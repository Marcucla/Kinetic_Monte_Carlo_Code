#ifndef __POINT_RES_H
#define __POINT_RES_H

#include "../StandardIncludes.h"
#include "../Settings.h"
#include "../context/Context.h"

class Simulation;

struct SimulationParameters;

struct PointRe
{
	/// The next defect in the same [111] crystal row.
	PointRe* next;
	/// The lattice position of the defect along the atomic row.
	int position;
	/// The type of defect.
	int type;
};

/**
 * The point defect cloud.
 */
class PointRes : public UsingContext
{
public:

	/// Constructor.
	PointRes(Context& context, Simulation* sim, SimulationParameters& params);

	/// Returns a reference to the simulation's global parameter set.
	const SimulationParameters& params() const { return _params; }

	SimulationParameters& params() { return _params; }

	/// Returns the global simulation object.
	Simulation& simulation() const { return *_sim; }

	/// Dumps the current point defect cloud to an output file.
	void dumpToVTK(std::ostream& stream, double simulationTime);

	/// Returns the list of point defects in the model.
	const std::map<std::pair<int,int>,PointRe*>& res() const { return _res; }

	/// Returns the list of point defects in the model.
	std::map<std::pair<int,int>, PointRe*>& res() { return _res; }

	/// Returns the world-space position of a point defect.
	Point3 getWorldPosition(int x, int y, int z) const;

	bool isSoluteOnTheDislocation(const Point3& p1, const Point3& p2);

	/********************************** Modification functions *********************************/

	/// Initializes the point defect cloud.
	void initialize();

	/// Adjusts the window around the dislocation within which point defects are explicitly modeled.
	void updateWindow();

	bool migration(PointRe* p, int x, int y, int z, int x1, int y1, int z1);

	double bindPointRes(const Point3& p0,const Point3& p1,double& deltaZ);

	bool isSoluteOnTheScrewDislocation(const Point3& p1, const Point3& p2);


private:

	/// Reference to the simulation's global parameters.
	SimulationParameters& _params;

	/// Pointer to the global simulation object.
	Simulation* _sim;

	/// The point defects in the model.
	std::map<std::pair<int,int>, PointRe*> _res;

	/// The memory pool for allocating point defect instances.
	boost::object_pool<PointRe> _rePool;

	/// Measures time spent computing stresses.
	//boost::timer::cpu_times _stressTime;

	/// The offset of the atomic crystal with respect to the lattice on which the dislocation moves.
	double _latticeOffsetX, _latticeOffsetY;
};

#endif // __POINT_RES_H
