#include "../StandardIncludes.h"
#include "../Settings.h"
#include "../context/Context.h"
#include "../simulation/Simulation.h"

using namespace std;

int main(int argc, char* argv[])
{
	namespace po = boost::program_options;

	try {
		// Create global context object.
		Context context;

		// Create main simulation object.
		Simulation simulation(context);

		// Define command line and configuration options:

		string config_file;
		po::options_description generic_options("Generic options");
		generic_options.add_options()
			("help", "Print help message")
			("config", po::value<string>(&config_file), "Name of configuration file.")
		;

		po::options_description param_options("Simulation parameters");
		param_options.add_options()
			("lattice_param", po::value<double>(&simulation.params().latticeParameter), "BCC lattice parameter (units: Angstrom)")
			("line_length", po::value<int>(&simulation.params().lineLength), "Screw dislocation line length (units: b)")
			("kink_drag_coefficient", po::value<double>(&simulation.params().kinkDragCoefficient), "Controls the mobility of kink segments (units: Pa*seconds)")
			("max_kink_motion_distance", po::value<double>(&simulation.params().maxKinkMotionDistance), "Maximum distance a kink segment may move within one simulation step (units: b)")
			("random_seed", po::value<unsigned int>(&simulation.params().randomSeed), "Seed for random number generator")
			("simulation_steps", po::value<int>(&simulation.params().numSimulationSteps), "Number of simulation steps to perform")
			("stress_total", po::value<double>(&simulation.params().stress_total), "External total stress(units: Pa)")
			("stress_xx", po::value<double>(&simulation.params().externalStress(0,0)), "External stress tensor component (units: Pa)")
			("stress_yy", po::value<double>(&simulation.params().externalStress(1,1)), "External stress tensor component (units: Pa)")
			("stress_zz", po::value<double>(&simulation.params().externalStress(2,2)), "External stress tensor component (units: Pa)")
			("stress_xy", po::value<double>(&simulation.params().externalStress(0,1)), "External stress tensor component (units: Pa)")
			("stress_xz", po::value<double>(&simulation.params().externalStress(0,2)), "External stress tensor component (units: Pa)")
			("stress_yz", po::value<double>(&simulation.params().externalStress(1,2)), "External stress tensor component (units: Pa)")
			("shear_modulus", po::value<double>(&simulation.params().shearModulus), "Elastic shear modulus of the material (units: Pa)")
			("poisson_ratio", po::value<double>(&simulation.params().poissonRatio), "Poisson's ratio of the material")
			("core_width", po::value<double>(&simulation.params().coreWidth), "Core width parameter of the non-singular dislocation theory (units: b)")
			("temperature", po::value<double>(&simulation.params().temperature), "Temperature (units: Kelvin)")
			("num_events", po::value<int>(&simulation.params().numEventsPerSegment), "Number of nucleation events to generate per screw segment")
			("attempt_frequency", po::value<double>(&simulation.params().nucleationAttemptFrequency), "Prefactor for kink-pair nucleation events (units: 1/seconds/b)")
			("peierls_stress", po::value<double>(&simulation.params().peierlsStress), "Peierls stress in formula for stress-dependent kink energy (units: Pa)")
			("kink_width", po::value<double>(&simulation.params().kinkWidth), "The width of a single kink in the screw direction (units: b)")
			("kpenergy_deltaH0", po::value<double>(&simulation.params().kpenergy_deltaH0), "Prefactor in formula for stress-dependent kink energy (units: eV)")
			("kpenergy_p", po::value<double>(&simulation.params().kpenergy_p), "Exponent in formula for stress-dependent kink energy")
			("kpenergy_q", po::value<double>(&simulation.params().kpenergy_q), "Exponent in formula for stress-dependent kink energy")
			("kpwidth_l0", po::value<double>(&simulation.params().kpwidth_l0), "Prefactor in formula for stress-dependent kink pair separation (units: b)")
			("kpwidth_w0", po::value<double>(&simulation.params().kpwidth_w0), "Parameter in formula for stress-dependent kink pair separation")
			("kpwidth_p", po::value<double>(&simulation.params().kpwidth_p), "Exponent in formula for stress-dependent kink pair separation")
			("kpwidth_q", po::value<double>(&simulation.params().kpwidth_q), "Exponent in formula for stress-dependent kink pair separation")
			("enable_local_stress", po::value<bool>(&simulation.params().enableLocalStress), "Enables/disables the calculation of local stresses")
			("cross_kink_lock_distance", po::value<double>(&simulation.params().crossKinkLockDistance), "Distance at which attraction/repulsion of kink pairs is determined (units: b)")
			("nonschmid_a1", po::value<double>(&simulation.params().nonschmid_a1), "First inverse prefactor in non-Schmid formula for computing RSS (no units). Default value is 1 for Schmid law.")
			("nonschmid_a2", po::value<double>(&simulation.params().nonschmid_a2), "Second parameter in non-Schmid formula for computing RSS (no units). Default value is 0 for Schmid law.")
			("nonschmid_a3", po::value<double>(&simulation.params().nonschmid_a3), "First inverse prefactor in non-Schmid formula for computing RSS (no units). Default value is 1 for Schmid law.")
			("nonschmid_tc", po::value<double>(&simulation.params().nonschmid_tc), "Second parameter in non-Schmid formula for computing RSS (no units). Default value is 0 for Schmid law.")
			("kink_diffusivity_coeff", po::value<double>(&simulation.params().kinkDiffusivityCoefficient), "The diffusivity of kink motion as linear function of temperature (units: m^2/s/K)")
			("point_defect_concentration", po::value<double>(&simulation.params().pointDefectConcentration), "The fraction of lattice sites occupied by point defects (valid range: 0 - 1).")
			("point_defect_window_radius", po::value<double>(&simulation.params().pointDefectWindowRadius), "Controls the size of the spatial window around the dislocation in which point defects are explicitly modeled (units: lattice constants).")
		;

		po::options_description output_options("Output options");
		output_options.add_options()
			("verbosity", po::value<int>(), "Verbosity level (0-3)")
			("position_output_interval", po::value<int>(&simulation.params().printEveryNTimesteps), "Controls how often the dislocation position is written to the position file.")
			("position_file", po::value<string>(&simulation.params().positionOutputFilename), "The output file for dislocation positions.")
			("vtk_output_interval", po::value<int>(&simulation.params().outputEveryNTimesteps), "Controls how often the dislocation configuration is output.")
			("vtk_output_file", po::value<string>(&simulation.params().vtkOutputFilename), "The output filename for dislocation snapshots.")
			("point_output_file", po::value<string>(&simulation.params().pointDefectOutputFilename), "The output filename for point defect snapshots.")
			("output_scaling_factor", po::value<double>(&simulation.params().outputScalingFactor), "Scaling applied to the Z coordinates when writing dislocation nodes to the VTK file.")
		;

		po::options_description cmdline_options;
		po::options_description config_file_options;
		po::options_description visible_options;
		cmdline_options.add(generic_options).add(param_options).add(output_options);
		config_file_options.add(param_options).add(output_options);
		visible_options.add(generic_options).add(output_options).add(param_options);

		// Parse command line options.
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv). options(cmdline_options).run(), vm);
		po::notify(vm);

		// Read in configuration file.
		if(!config_file.empty()) {
			ifstream ifs(config_file.c_str());
			if(!ifs)
				context.error("Failed to open configuration file: %s", config_file.c_str());
			po::store(parse_config_file(ifs, config_file_options), vm);
			po::notify(vm);
		}

		// Set verbosity level.
		if(vm.count("verbosity"))
			context.setVerbosityLevel(vm["verbosity"].as<int>());

		if(vm.count("version")) {
			cout << "Kinetic Monte-Carlo simulation code" << endl;
			cout << "Version: " << PROGRAM_VERSION_STRING << endl;
			return 0;
		}

		if(vm.count("help") || argc == 1) {
			cout << endl << "Usage: DislocationKMC [options]" << endl;
			cout << visible_options << endl;
			return 0;
		}

		// Verify simulation parameters.
		if(!vm.count("lattice_param") || simulation.params().latticeParameter <= 0) context.error("Missing or invalid parameter: lattice_param");
		if(!vm.count("line_length") || simulation.params().lineLength <= 0) context.error("Missing or invalid parameter: line_length");
		if(!vm.count("kink_drag_coefficient") || simulation.params().kinkDragCoefficient <= 0) context.error("Missing or invalid parameter: kink_drag_coefficient");
		if(!vm.count("max_kink_motion_distance") || simulation.params().maxKinkMotionDistance < 0) context.error("Missing or invalid parameter: max_kink_motion_distance");
		if(!vm.count("random_seed")) context.error("Missing or invalid parameter: random_seed");
		if(!vm.count("simulation_steps") || simulation.params().numSimulationSteps < 0) context.error("Missing or invalid parameter: simulation_steps");

		// Setup simulation.
		simulation.initialize();

		// Run kMC simulation.
		context.msgLogger(VERBOSITY_NORMAL) << "Starting simulation run" << endl;
		simulation.printSimulationInfo();
		simulation.outputDislocations();
		for(int step = 1; step <= simulation.params().numSimulationSteps; step++) {
			context.msgLogger(VERBOSITY_HIGH) << "***************** Step " << step << " *****************" << endl;
			simulation.performIteration();
			simulation.printSimulationInfo();
			simulation.outputDislocations();
			simulation.outputPointDefects();
		}

		return 0;
	}
	catch(const std::bad_alloc& ex) {
		cerr << endl << "ERROR: Out of memory" << endl;
		return 1;
	}
	catch(const exception& ex) {
		cerr << endl << "ERROR: " << ex.what() << endl;
		return 1;
	}

	return 0;
}
