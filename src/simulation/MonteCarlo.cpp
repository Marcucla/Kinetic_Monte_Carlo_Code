#include "Simulation.h"

using namespace std;

/******************************************************************************
* Performs a single iteration of the kMC algorithm.
*****************************************************************************/
double Simulation::performIteration()
{
	// Assign velocities to kink segments and nodes.
	network().calculateKinkVelocities(randomGenerator());

	// Determine how long kinks can move until they reach the maximum migration distance set by the user.
	double maxMotionTime = CAFLOAT_MAX;
	if(params().maxKinkMotionDistance > 0) {
		maxMotionTime = network().minTimeForDistance(params().maxKinkMotionDistance);
	}

	// Determine how long kinks can move until they hit another kink segment.
	double freeMotionTime = network().freeMotionTime();


	// The kink migration time.
	double t_mig1 = min(maxMotionTime, freeMotionTime);

	//  Determine how long kinks can move until they hit solute atoms.
	double freeBindTime = network().minTimeForBind(t_mig1);

	// The kink migration time.
	double t_mig = min(t_mig1, freeBindTime);

	// Generate list of kink-pair nucleation events and compute rates.
	vector<KinkPairNucleationEvent> events;
	vector<PointDefectsEventList> pevents;
	vector<DislocationBindingEventList> devents;
	double totalRate = generateNucleationEventList(events);
	double totalRatePointDefects = generatePointDefectsEventList(pevents);
	double totalRateDislocation = generateDislocationBindingEventList(devents);
	_simulationNucleationRate = totalRate;
	_simulationDiffusionRate = totalRatePointDefects;
	_simulationDetachmentRate = totalRateDislocation;
	//context().msgLogger(VERBOSITY_NORMAL) << "Nucleation: " << totalRate << " bindingrate: " << totalRateDislocation << "diffusion: " << totalRatePointDefects << endl;
	//context().msgLogger(VERBOSITY_NORMAL) << "Nucleation: " << totalRate << endl;
	double totalRateall = totalRate + totalRatePointDefects + totalRateDislocation;
	double t_nuc;

	if(totalRateall > 0) {
		// Generate nucleation time from exponential distribution.
		boost::exponential_distribution<long double> expDist(totalRateall);
		boost::uniform_01<random_generator_type> uniform01(randomGenerator());
		t_nuc = expDist(uniform01);
	}
	else {
		if(t_mig == CAFLOAT_MAX)
			context().error("Total kink nucleation rate became zero and there are no kinks to move at timestep %i.", simulationStep());
		t_nuc = CAFLOAT_MAX;
	}

	context().msgLogger(VERBOSITY_HIGH) << "Free kink motion time: " << freeMotionTime <<
			"  Maximum kink motion time: " << maxMotionTime <<
			"  Kink migration time: " << t_mig <<
			"  Nucleation time: " << t_nuc <<
			"  Total rate: " << totalRate << endl;

	bool modifiedDislocationNetwork = false;
	bool nucleatedKinkPair = false;
	bool modifiedPointDefectsClouds = false;
	double deltaTime;

	//context().msgLogger(VERBOSITY_NORMAL) << "t_mig: " << t_mig << " t_nuc: " << t_nuc << endl;
	
	if(t_mig <= t_nuc) {
		// All kinks move with their current velocities for a time period t_mig.
		context().msgLogger(VERBOSITY_HIGH) << "Propagating kink over time period: " << t_mig << " [s]" << endl;
		network().propagateKinks(t_mig);
		//context().msgLogger(VERBOSITY_NORMAL) << "t_mig: " << t_mig << " t_nuc: " << t_nuc << endl;
		modifiedDislocationNetwork = true;
		deltaTime = t_mig;
	}
	else {
		// All kinks move with their current velocities for a time period t_nuc.
		if(t_mig != CAFLOAT_MAX) {
			context().msgLogger(VERBOSITY_HIGH) << "Propagating kink over time period: " << t_nuc << " [s]" << endl;
			network().propagateKinks(t_nuc);
			modifiedDislocationNetwork = true;
		}

		// Generate random number in range [0,totalRate).
		boost::uniform_real<> uniformDist(0, totalRateall);
		double r = uniformDist(randomGenerator());
		SIMULATION_ASSERT(r >= 0 && r < totalRateall);

		// Pick an event based on random number.
		double cumRate = 0;
		if (r <= totalRate){
		auto event = events.begin();
		for(;; ++event) {
			SIMULATION_ASSERT(event != events.end());
			cumRate += event->rate;
			if(cumRate >= r) break;
		}
		
		SIMULATION_ASSERT(cumRate > 0.0);

		context().msgLogger(VERBOSITY_HIGH) << "Executing nucleation event:" << endl;
		context().msgLogger(VERBOSITY_HIGH) << "  Kink direction: " << event->kinkDirection << endl;
		context().msgLogger(VERBOSITY_HIGH) << "  Kink pair position: " << event->kinkPairPosition << endl;
		context().msgLogger(VERBOSITY_HIGH) << "  Kink separation: " << event->kinkPairWidth << " [b]" << endl;

		// Then execute kMC event.
		if(executeNucleationEvent(*event))
			modifiedDislocationNetwork = true;

		// Return physical time of this kMC step.
		deltaTime = 1.0 / totalRateall;

		nucleatedKinkPair = true;
		}
		else if(r <= (totalRate + totalRatePointDefects))
		{
		cumRate = totalRate;
		auto pevent = pevents.begin();
		for(;; ++pevent) {
			SIMULATION_ASSERT(pevent != pevents.end());
			cumRate += pevent->rate;
			if(cumRate >= r) break;
		}
		if(executePointDefectsEvent(*pevent))
			modifiedPointDefectsClouds = true;

		//context().msgLogger(VERBOSITY_NORMAL) << "modified: " << modifiedPointDefectsClouds << endl;
		deltaTime = 1.0 / totalRateall;

		}

		else
		{
		cumRate = totalRate + totalRatePointDefects;
		auto devent = devents.begin();
		double kT = params().temperature * 8.6173324e-5;
		for(;; ++devent) {
			SIMULATION_ASSERT(devent != devents.end());
			cumRate += devent->rate;
			//cumRate += params().nucleationAttemptFrequency / params().numEventsPerSegment * exp(-0.175 / kT);
			if(cumRate >= r) break;
		}

		if(network().executeBindingEvent(*devent))
			modifiedDislocationNetwork = true;

		deltaTime = 1.0 / totalRateall;
	    }
	}

#ifdef DEBUG_SIMULATION
	// Make sure dislocation data structure is in good order.
	network().validate(false);
#endif

	// Find attracting cross-kinks that form a locked configuration.
	if(network().detectLockedCrossKinks())
		modifiedDislocationNetwork = true;

	// Remove redundant segments and handle dislocation reactions.
	if(modifiedDislocationNetwork) {
		network().cleanupNetwork();

#ifdef DEBUG_SIMULATION
		// Make sure dislocation data structure is in good order.
		network().validate(true);
#endif
	}
	// Update point defect window when dislocation has moved via kink-pair mechanism.
	//if(nucleatedKinkPair)
	pointDefects().updateWindow();


	_simulationStep++;
	_simulationTime += deltaTime;

	//context().msgLogger(VERBOSITY_NORMAL) << "deltaTime: " << deltaTime << endl;

	// Return physical time increment of this simulation step.
	return deltaTime;

}

