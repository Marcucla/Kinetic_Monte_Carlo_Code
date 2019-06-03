#include "Simulation.h"

using namespace std;

/******************************************************************************
* Constructor.
*****************************************************************************/
Simulation::Simulation(Context& _context) : UsingContext(_context), _dislocNet(_context, this, _params),_pointDefects(_context, this, _params)
{
	_params.latticeParameter = 0;
	_params.lineLength = 0;
	_params.kinkDragCoefficient = 0;
	_params.maxKinkMotionDistance = 0;
	_params.externalStress = SymmetricTensor2(0);
	_params.numSimulationSteps = 0;
	_params.randomSeed = 1;
	_params.printEveryNTimesteps = 100;
	_params.outputEveryNTimesteps = 0;
	_params.vtkOutputFilename = "dislocation.vtk";
	_params.outputScalingFactor = 1;
	_params.shearModulus = 0;
	_params.poissonRatio = 0;
	_params.coreWidth = 0;
	_params.temperature = 0;
	_params.numEventsPerSegment = 0;
	_params.nucleationAttemptFrequency = 0;
	_params.peierlsStress = 0;
	_params.kinkWidth = 0;
	_params.kpwidth_l0 = 0;
	_params.kpwidth_w0 = 0;
	_params.kpwidth_p = 0;
	_params.kpwidth_q = 0;
	_params.kpenergy_deltaH0 = 0;
	_params.kpenergy_p = 0;
	_params.kpenergy_q = 0;
	_params.enableLocalStress = true;
	_params.crossKinkLockDistance = 0;
	_params.nonschmid_a1 = 1;
	_params.nonschmid_a2 = 0;
	_params.nonschmid_a3 = 0;
	_params.nonschmid_tc = 0;
	_params.kinkDiffusivityCoefficient = 0;
	_params.stress_total = 0;

	_params.pointDefectConcentration = 0;
	_params.pointDefectWindowRadius = 0;

	_simulationTime = 0;
	_simulationStep = 0;

	_simulationDiffusionRate = 0;
	_simulationNucleationRate = 0;
	_simulationPdefectMeanSquareMotion = 0;
	_simulationPdefectMeanSquareMotionReduced = 0;
	_number = 0;
}

/******************************************************************************
* Initializes the simulation.
* Computes dependent parameters from user-supplied parameters.
*****************************************************************************/
void Simulation::initialize()
{
	// Validate parameters.
	SIMULATION_ASSERT(params().latticeParameter > 0);
	SIMULATION_ASSERT(params().peierlsStress > 0);

	// Specify primitive cell of lattice - in this case bcc.
	// Burgers vector of screw dislocation coincides with third edge vector of unit cell.
	_params.unitCell.setColumn(0, Vector3(sqrt(6.0) / 3.0, 0, -sqrt(3.0) / 6.0));
	_params.unitCell.setColumn(1, Vector3(sqrt(6.0) / 6.0, sqrt(2.0) / 2.0, sqrt(3.0) / 6.0));
	_params.unitCell.setColumn(2, Vector3(0, 0, sqrt(3.0) / 2.0));

	// Verify volume of bcc primitive cell.
	SIMULATION_ASSERT(fabs(fabs(_params.unitCell.determinant()) - 0.5) < CAFLOAT_EPSILON);

	_params.peierlsStress = 1.07 * 1e9/(_params.pointDefectConcentration + 0.51);

	// Scale lattice cell with lattice parameter.
	_params.unitCell *= _params.latticeParameter;
	_params.inverseUnitCell = _params.unitCell.inverse();

	// Burgers vector:
	_params.b = NodalVector(0,0,1);
	_params.blength = _params.latticeParameter * sqrt(3.0) / 2.0;
	_params.b_spatial = Vector3(0, 0, _params.blength);

	// Line length:
	SIMULATION_ASSERT(_params.lineLength > 0);
	_params.pbcLength = _params.blength * _params.lineLength;

	// List of kink directions (expressed in terms of primitive cell vectors. all are <111>/2 type).
	_params.kinkDirections.resize(6);
	_params.kinkDirections[0].h = NodalVector(+1,  0, +1.0/3.0);
	_params.kinkDirections[1].h = NodalVector(-1,  0, -1.0/3.0);
	_params.kinkDirections[2].h = NodalVector( 0, +1, -1.0/3.0);
	_params.kinkDirections[3].h = NodalVector( 0, -1, +1.0/3.0);
	_params.kinkDirections[4].h = NodalVector(-1,  1, -2.0/3.0);
	_params.kinkDirections[5].h = NodalVector( 1, -1, +2.0/3.0);
	_params.kinkDirections[0].reverseDirection = 1;
	_params.kinkDirections[1].reverseDirection = 0;
	_params.kinkDirections[2].reverseDirection = 3;
	_params.kinkDirections[3].reverseDirection = 2;
	_params.kinkDirections[4].reverseDirection = 5;
	_params.kinkDirections[5].reverseDirection = 4;
	for(int i = 0; i < _params.kinkDirections.size(); i++) {
		_params.kinkDirections[i].h_spatial = _params.unitCell * _params.kinkDirections[i].h;
		_params.kinkDirections[i].h_mag = Length(_params.kinkDirections[i].h_spatial);
		_params.kinkDirections[i].theta = atan2(_params.kinkDirections[i].h_spatial.Y, _params.kinkDirections[i].h_spatial.X);
		std::cerr << i << " " << _params.kinkDirections[i].h_spatial << " " << rad2deg(_params.kinkDirections[i].theta) << std::endl;
		_params.kinkDirections[i].normal = Normalize(Vector3(-_params.kinkDirections[i].h_spatial.Y, _params.kinkDirections[i].h_spatial.X, 0));
		SIMULATION_ASSERT(fabs(_params.kinkDirections[i].h_mag - _params.latticeParameter*sqrt(6.0)/3.0) < 1e-12);
	}

	// Initialize random number generator.
	_random.seed(params().randomSeed);

	// Output simulation parameters to log file for reference.
	context().msgLogger(VERBOSITY_NORMAL) << "Program version:                  " << PROGRAM_VERSION_STRING << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "---------- Simulation parameters ------------" << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Lattice parameter [A]:               " << params().latticeParameter << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Unit cell [A]:                       " << params().unitCell.row(0) << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "                                     " << params().unitCell.row(1) << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "                                     " << params().unitCell.row(2) << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Burgers vector:                      " << params().b << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Burgers vector magnitude [A]:        " << params().blength << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Line length [b]:                     " << params().lineLength << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Periodicity length [A]:              " << params().pbcLength << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Temperature [K]:                     " << params().temperature << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink drag coefficient [Pa*s]  :      " << params().kinkDragCoefficient << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Maximum kink-motion distance [b]:    " << params().maxKinkMotionDistance << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Number of events per segment:        " << params().numEventsPerSegment << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink-pair attempt frequency [1/s/b]: " << params().nucleationAttemptFrequency << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "External total stress [Pa]:          " << params().stress_total << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "External stress [Pa]:                " << params().externalStress.row(0) << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "                                     " << params().externalStress.row(1) << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "                                     " << params().externalStress.row(2) << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Shear modulus [Pa]:                  " << params().shearModulus << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Poisson's ratio:                     " << params().poissonRatio << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Core width [b]:                      " << params().coreWidth << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink height [A]:                     " << params().kinkDirections[0].h_mag << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Peierls stress [Pa]:                 " << params().peierlsStress << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink width [b]:                      " << params().kinkWidth << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink pair separation l0 [b]:         " << params().kpwidth_l0 << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink pair separation w0:             " << params().kpwidth_w0 << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink pair separation p:              " << params().kpwidth_p << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink pair separation q:              " << params().kpwidth_q << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink pair energy Delta H0 [eV]:      " << params().kpenergy_deltaH0 << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink pair energy p:                  " << params().kpenergy_p << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink pair energy q:                  " << params().kpenergy_q << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Non-Schmid parameter a1:             " << params().nonschmid_a1 << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Non-Schmid parameter a2:             " << params().nonschmid_a2 << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Non-Schmid parameter a3:             " << params().nonschmid_a3 << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Non-Schmid parameter tc:             " << params().nonschmid_tc << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Kink-diffusivity coeff. [m^2/s/K]:   " << params().kinkDiffusivityCoefficient << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Enable local stress calculation:     " << params().enableLocalStress << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Cross-kink lock distance [b]:        " << params().crossKinkLockDistance << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Random number seed:                  " << params().randomSeed << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Number of simulation steps:          " << params().numSimulationSteps << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Print interval:                      " << params().printEveryNTimesteps << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "Output interval:                     " << params().outputEveryNTimesteps << endl;
	context().msgLogger(VERBOSITY_NORMAL) << "---------------------------------------------" << endl;

	// Set up initial dislocation configuration.
	network().initialize();

	pointDefects().initialize();

	pointDefects().updateWindow();

	// Open log file.
	if(params().positionOutputFilename.empty() == false) {
		_dislocPositionFile.open(params().positionOutputFilename.c_str());
		if(!_dislocPositionFile.is_open())
			context().error("Failed to open output log file for writing: %s", params().positionOutputFilename.c_str());

		_dislocPositionFile
			<< "# Program version:                     " << PROGRAM_VERSION_STRING << endl
			<< "# Line length [b]:                     " << params().lineLength << endl
			<< "# Temperature [K]:                     " << params().temperature << endl
			<< "# Sigma total [Pa]:                    " << params().stress_total << endl
			<< "# Sigma xz [Pa]:                       " << params().externalStress(0,2) << endl
			<< "# Sigma yz [Pa]:                       " << params().externalStress(1,2) << endl
			<< "# Random seed:                         " << params().randomSeed << endl
			<< "# Kink drag coefficient [Pa*s] :       " << params().kinkDragCoefficient << endl
			<< "# Maximum kink-motion distance [b]:    " << params().maxKinkMotionDistance << endl
			<< "# Number of events per segment:        " << params().numEventsPerSegment << endl
			<< "# Kink-pair attempt frequency [1/s/b]: " << params().nucleationAttemptFrequency << endl
			<< "# Cross-kink lock distance [b]:        " << params().crossKinkLockDistance << endl
			<< "# Enable local stress calculation:     " << params().enableLocalStress << endl
			<< "# Core width [b]:                      " << params().coreWidth << endl
			<< "# Peierls stress [Pa]:                 " << params().peierlsStress << endl
			<< "# Kink width [b]:                      " << params().kinkWidth << endl
			<< "# Shear modulus [Pa]:                  " << params().shearModulus << endl
			<< "# Poisson's ratio:                     " << params().poissonRatio << endl;

		_dislocPositionFile << "#" << endl << "# Columns:" << endl
				<< "# 1. simulation_step" << endl
				<< "# 2. time [s]" << endl
				<< "# 3. position_X [A]" << endl
				<< "# 4. position_Y [A]" << endl
				<< "# 5. line_width_X [A]" << endl
				<< "# 6. line_width_Y [A]" << endl
				<< "# 7. loop_count" << endl
				<< "# 8. segment_count" << endl
				<< "# 9. simulationDiffusionRate" << endl
				<< "# 10. simulationNucleationRate" << endl
				<< "# 10. simulationDetachmentRate" << endl
				<< "# 11. simulationPdefectMeanSquareMotion" << endl;
	}
}
