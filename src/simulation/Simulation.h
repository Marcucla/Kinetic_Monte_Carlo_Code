#ifndef __SIMULATION_H
#define __SIMULATION_H

#include "../StandardIncludes.h"
#include "../Settings.h"
#include "../context/Context.h"
#include "../dislocations/Dislocations.h"
#include "../pointdefects/PointDefects.h"

/*
 * Stores the precomputed information
 * on each possible kink direction.
 */


struct KinkDirection
{
	Vector3		h_spatial;			// The direction vector of the kink in spatial coordinates.
	Vector3		normal;				// The normal of the kink plane in spatial coordinates.
	NodalVector h;					// The direction vector of the kink in lattice coordinates.
	double		h_mag;				// The spatial magnitude of the kink [units: A].
	double		theta;				// Angle of kink plane with positive X axis.
	int			reverseDirection;	// The direction index that is reverse to this one.
};

/*
 * Stores all simulation parameters.
 */
struct SimulationParameters
{
	double latticeParameter;				// Lattice constant (Angstrom units).
	Matrix3 unitCell;						// Primitive unit cell of the crystal lattice (Angstrom units).
	Matrix3 inverseUnitCell;				// Inverse matrix for transforming spatial coordinates to lattice space coordinates.

	NodalVector b;							// The Burgers vector in terms of the primitive lattice cell.
	Vector3 b_spatial;						// The Burgers vector in spatial coordinates.
	int lineLength;							// Length of (straight) screw dislocation in units of the Burgers vector. This determines the simulation size.
	double blength;							// Length of Burgers vector (Angstrom units).
	double pbcLength;						// Size of periodic box = blength * lineLength.
	std::vector<KinkDirection> kinkDirections;// List of possible kink nucleation directions (and distances)

	double kinkDragCoefficient;				// Controls the mobility of kink segments (units: Pa * seconds).
	double maxKinkMotionDistance;			// Maximum distance kinks may move during one simulation step (Burgers vector units).
	double stress_total;
	SymmetricTensor2 externalStress;		// The external stress that drives the dislocation (units: Pa).

	double shearModulus;					// The elastic shear modulus (mu) of the material (units: Pa).
	double poissonRatio;					// The elastic Poisson ratio (nu) of the material.
	double coreWidth;						// The core width parameter of the non-singular dislocation theory (units: b).

	double temperature;						// The temperature parameter of the simulation (units: Kelvin).
	int numEventsPerSegment;				// The number of nucleation events to generate per screw segment.
	double nucleationAttemptFrequency;		// The prefactor for kink-pair nucleation events (units: 1/seconds/b).

	int numSimulationSteps;					// Number of simulation steps to perform.
	unsigned int randomSeed;				// The random number generator seed.

	int printEveryNTimesteps;				// Controls how often the current dislocation position is written to the log file.
	std::string positionOutputFilename;		// Output filename for dislocation position data.

	int outputEveryNTimesteps;				// Controls how often the dislocation configuration is written to the output file.
	std::string vtkOutputFilename;			// Output filename for dislocation snapshots.
	std::string pointDefectOutputFilename;	// Output filename for point defect cloud snapshots.
	double outputScalingFactor;				// Scaling applied to the Z coordinates when writing dislocation nodes to the VTK file.

	double peierlsStress;					// Peierls stress in formula for stress-dependent kink energy (units: Pa).

	double kinkWidth;						// The width of a single kink (units: b).

	double kpwidth_l0;						// Prefactor in formula for stress-dependent kink pair separation (units: b).
	double kpwidth_w0;						// Parameter in formula for stress-dependent kink pair separation (no units).
	double kpwidth_p;						// Exponent parameter in formula for stress-dependent kink pair separation (no units).
	double kpwidth_q;						// Exponent parameter in formula for stress-dependent kink pair separation (no units).

	double kpenergy_deltaH0;				// Prefactor in formula for stress-dependent kink energy (units: eV).
	double kpenergy_p;						// Exponent parameter in formula for stress-dependent kink energy (no units).
	double kpenergy_q;						// Exponent parameter in formula for stress-dependent kink energy (no units).

	double crossKinkLockDistance;			// If two kinks, separated by a distance in the screw direction less than this, feel an attractive
											// force then they are considered a compound object that moves in a synchronized fashion (units: b).

	bool enableLocalStress;					// Enables/disables the calculation of local stresses (in addition to the global applied stress).

	double nonschmid_a1;					// First inverse prefactor in non-Schmid formula for computing the RSS (no units).
	double nonschmid_a2;
	double nonschmid_a3;
	double nonschmid_tc;					// Second parameter in non-Schmid formula for computing the RSS (no units).

	double kinkDiffusivityCoefficient;		// The diffusivity of kink motion as linear function of temperature. (default: 0) (units: meter^2/(seconds*kelvin))
	double pointDefectConcentration;		// The fraction of lattice sites occupied by point defects.
	double pointDefectWindowRadius;			// Controls the size of the spatial window around the dislocation in which point defects are explicitly modeled (units: lattice constants).
};

/*
 * Stores the rate of a specific kink-pair nucleation event.
 */
struct KinkPairNucleationEvent
{
	/// The H-segment on which the nucleation occurs.
	SegmentHandle segment;

	/// The KP center of mass of the kink pair on the screw segment (in units of the Burgers vector).
	double kinkPairPosition;

	/// The separation of the kink pair (in units of the Burgers vector).
	double kinkPairWidth;

	/// The direction of the kink pair. This is an index into the global list of possible kink directions
	/// defined in the SimulationParameters structure.
	int kinkDirection;

	/// The rate of this event [1/s].
	double rate;
};

struct PointDefectsEventList
{
	double rate;
	PointDefect* p;
	double x;   //
	double y;
	double z;
	double x1;
	double y1;
	double z1;
	int t;
	double x_line;
	double y_line;
	double x1_line;
	double y1_line;
};


/**
 * The main simulation object.
 */
class Simulation : public UsingContext
{
public:

	/// The type of pseudo-random number generator used by the simulation.
	typedef boost::mt19937 random_generator_type;

public:

	/// Constructor.
	Simulation(Context& context);

	/// Returns a reference to the parameter set.
	const SimulationParameters& params() const { return _params; }

	/// Returns a reference to the parameter set.
	SimulationParameters& params() { return _params; }

	/// Initializes the simulation.
	void initialize();

	/// Performs a single iteration of the kMC algorithm.
	double performIteration();

	/// Prints a line to the log file.
	void printSimulationInfo();

	/// Dumps the current dislocation structure to a file.
	void outputDislocations();

	void outputPointDefects();

	/// Returns a reference to the dislocation network.
	DislocationNetwork& network() { return _dislocNet; }

	/// Returns a reference to the point defect cloud.
	PointDefects& pointDefects() { return _pointDefects; }

	/// Returns the global random number generator.
	random_generator_type& randomGenerator() { return _random; }

	/// Return the number of simulation iterations performed so far.
	int simulationStep() const { return _simulationStep; }

	/// Returns the simulated time.
	double simulationTime() const { return _simulationTime; }

	int number() const {return _number;}
	///9.6[MODIFICATION FOR MOVING SOLUTE

	double simulationPdefectMeanSquareMotionReduced() const { return _simulationPdefectMeanSquareMotionReduced; }
	///

	void modifySimulationPdefectMeanSquareMotionReduced(double& msmr) { _simulationPdefectMeanSquareMotionReduced = msmr; }

	///
	void modifySimulationPdefectMeanSquareMotion(double& msm) { _simulationPdefectMeanSquareMotion = msm; }

	///9.6 MODIFICATION FOR MOVING SOLUTE]
	double simulationPdefectMeanSquareMotion() const { return _simulationPdefectMeanSquareMotion; }


	/// Generates a list of all possible nucleation events and computes their rates.
	double generateNucleationEventList(std::vector<KinkPairNucleationEvent>& events);

	double generatePointDefectsEventList(std::vector<PointDefectsEventList>& pevents);

	double generateDislocationBindingEventList(std::vector<DislocationBindingEventList>& devents);

	/// Calculates the activation energy for a kink-pair nucleation event.
	bool calculateActivationEnergy(const KinkPairNucleationEvent& event, double& sustainableKinkSeparation, double& activationEnergy, double& frequencyPrefactor);

	bool calculateActivationEnergyPointDefects(const PointDefectsEventList& pevent, double& activationEnergy, double& frequencyPrefactor, const Point3& Pointp);

	/// Executes the selected kink-pair nucleation event by modifying the current dislocation configuration.
	bool executeNucleationEvent(const KinkPairNucleationEvent& event);

	bool executePointDefectsEvent(const PointDefectsEventList& pevent);

	//bool executeBindingEvent(const SegmentHandle& devent);



private:

	/// The input parameter set.
	SimulationParameters _params;

	/// The dislocation network.
	DislocationNetwork _dislocNet;

	/// The point defect cloud.
	PointDefects _pointDefects;

	/// The global random number generator.
	random_generator_type _random;

	/// The number of simulation iterations performed so far.
	int _simulationStep;

	/// The simulated time.
	double _simulationTime;

	///9.6;[MODIFICATION FOR MOVING SOLUTE
	double _simulationPdefectMeanSquareMotionReduced;

	///9.6;
	double _simulationPdefectMeanSquareMotion;

	///number of pointdefects
	int _number;

	///9.6version nucleation rate
	double _simulationNucleationRate;

	///9.6version diffusion rate ]
	double _simulationDiffusionRate;

	///
	double _simulationDetachmentRate;
	/// The log file to which the current dislocation position is written.
	std::ofstream _dislocPositionFile;
};

#endif // __SIMULATION_H
