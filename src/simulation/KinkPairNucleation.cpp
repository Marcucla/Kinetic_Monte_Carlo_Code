#include "Simulation.h"

using namespace std;

/******************************************************************************
* Generates a list of all possible nucleation events and computes their rates.
*****************************************************************************/
double Simulation::generateNucleationEventList(vector<KinkPairNucleationEvent>& events)
{
	double totalRate = 0;

	// Boltzmann constant times temperature (in [eV] units).
	double kT = params().temperature * 8.6173324e-5;

	for(SegmentIterator segment(network()); !segment.atEnd(); ++segment) {

		// Kinks can only be nucleated on screw segments, and not on edge segments.
		if(segment->isScrew() == false)
			continue;

		double segmentWidth = segment->getHSegmentLength();

		// Kinks are not nucleated on degenerate (zero length) screw segments.
		if(segmentWidth <= CAFLOAT_EPSILON)
			continue;

		// We don't allow kink pairs to be nucleated on very short screw segments.
		if(segmentWidth <= params().kinkWidth)
			continue;

		// Determine interval within which new KPs can be nucleated.
		double nucleation_interval_start = segment->node1()->pos().Z + params().kinkWidth / 2;
		double nucleation_interval_end   = segment->node1()->pos().Z + segmentWidth - params().kinkWidth / 2;

		// Random number distribution used to sample the nucleation positions along the current segment.
		boost::uniform_real<> positionDistribution(nucleation_interval_start, nucleation_interval_end);

		for(int direction = 0; direction < params().kinkDirections.size(); direction++) {

			// Generate several nucleation events per screw segment.
			for(int i = 0; i < params().numEventsPerSegment; i++) {

				KinkPairNucleationEvent event;
				event.segment = *segment;
				event.kinkDirection = direction;

				// Sample kink pair positions along screw segment.
				event.kinkPairPosition = positionDistribution(randomGenerator());

				// Compute activation energy (in eV) of nucleation event and sustainable kink separation.
				double activationEnergy;
				double frequencyPrefactor;
				if(!calculateActivationEnergy(event, event.kinkPairWidth, activationEnergy, frequencyPrefactor))
					continue;

				SIMULATION_CHECK_VALUE(activationEnergy);
				SIMULATION_CHECK_VALUE(event.kinkPairWidth);
				SIMULATION_ASSERT(event.kinkPairWidth > 0);
				if(activationEnergy <= 0) {
					context().msgLogger(VERBOSITY_HIGH) << "Kink pair width: " << event.kinkPairWidth << endl;
					context().msgLogger(VERBOSITY_HIGH) << "Kink pair nucleation activation energy: " << activationEnergy << endl;
					Point3 nucleationSite = params().unitCell * NodalPosition(event.segment->node1()->pos().X, event.segment->node1()->pos().Y, event.kinkPairPosition);
					context().error("Activation energy for kink pair nucleation became non-positive. Nucleation site is at spatial position (%f %f %f). Kink direction is %i.",
							nucleationSite.X, nucleationSite.Y, fmod(nucleationSite.Z + params().pbcLength, params().pbcLength), event.kinkDirection);
				}

				// Compute attempt frequency.
				// This is proportional to the segment's length.
				// Normalize it by the number of MC events we generate per segment.
				double frequency = frequencyPrefactor * params().nucleationAttemptFrequency / params().numEventsPerSegment;
				SIMULATION_ASSERT(frequency >= 0);
				if(frequency == 0)
					continue;

				// Compute nucleation rate.
				event.rate = frequency * exp(-activationEnergy / kT);
				SIMULATION_ASSERT(event.rate >= 0);

#if 0
				context().msgLogger(VERBOSITY_HIGH) << "Generated nucleation event:" << endl;
				context().msgLogger(VERBOSITY_HIGH) << "  Kink direction: " << event.kinkDirection << endl;
				context().msgLogger(VERBOSITY_HIGH) << "  Kink separation: " << event.kinkPairWidth << " [b]" << endl;
				context().msgLogger(VERBOSITY_HIGH) << "  Activation energy: " << activationEnergy << " [eV]" << endl;
				context().msgLogger(VERBOSITY_HIGH) << "  Activation rate: " << event.rate << " [1/s]" << endl;
#endif

				// Insert event into catalog.
				totalRate += event.rate;
				events.push_back(event);
			}
		}
	}

	return totalRate;
}

/******************************************************************************
* Calculates the activation energy for a kink pair nucleation event.
*****************************************************************************/
bool Simulation::calculateActivationEnergy(const KinkPairNucleationEvent& event, double& sustainableKinkSeparation, double& activationEnergy, double& frequencyPrefactor)
{
	// Site of nucleation where kink pair should be nucleated.
	Point3 nucleationSite = params().unitCell * NodalPosition(event.segment->node1()->pos().X, event.segment->node1()->pos().Y, event.kinkPairPosition);

	double resolvedShearStress;
	double s;
	if(params().nonschmid_a1 == 1 && params().nonschmid_a2 == 0) {
		// Implementation of standard Schmid law:

		// Compute local resolved shear stress in kink direction.
		resolvedShearStress = network().calculateResolvedShearStress(nucleationSite, event.kinkDirection);
	}
	else {
		// Implementation of non-Schmid corrections:

		// Compute local stress tensor.
		SymmetricTensor2 stress = network().calculateLocalStress(nucleationSite);

		// MRSS plane direction.
		double theta_MRSS = atan2(-stress(0,2), stress(1,2));
		// MRSS value.
		double sigma_MRSS = sqrt(stress(0,2)*stress(0,2) + stress(1,2)*stress(1,2));

		// Kink plane direction.
		double theta_k = params().kinkDirections[event.kinkDirection].theta;

		double sigma1 = -stress(0,2)*sin(theta_k)+stress(1,2)*cos(theta_k);

		double sigma2 = -stress(0,2)*sin(theta_k+deg2rad(60))+stress(1,2)*cos(theta_k+deg2rad(60));

		double sigma3 = (stress(1,1)-stress(0,0))*sin(theta_k)*cos(theta_k)+sin(theta_k)*sin(theta_k)*stress(0,1)-cos(theta_k)*cos(theta_k)*stress(1,0);

		double sigma4 = (stress(1,1)-stress(0,0))*sin(theta_k+deg2rad(60))*cos(theta_k+deg2rad(60))+sin(theta_k+deg2rad(60))*sin(theta_k+deg2rad(60))*stress(0,1)-cos(theta_k+deg2rad(60))*cos(theta_k+deg2rad(60))*stress(1,0);


		// Compute signed angle between kink plane and MRSS plane.
		double chi = theta_MRSS - theta_k;

		// Compute local resolved shear stress in kink direction.
		s = (sigma1 + params().nonschmid_a1 * sigma2)/(params().nonschmid_tc * params().peierlsStress * (2/(1+exp(2*((params().nonschmid_a2 * sigma3 + params().nonschmid_a3 * sigma4)/(params().nonschmid_tc * params().peierlsStress))))));
	}

	// Cannot nucleate sustainable kink pair if stress is non-positive.
	if(s <= 0)
		return false;

	// Check if local stress exceeds Peierls stress.
	if(s >= 1.0) {
		context().msgLogger(VERBOSITY_NORMAL) << "WARNING: Local resolved shear stress exceeds Peierls stress at spatial position (" <<
				nucleationSite.X << " " << nucleationSite.Y << " " << (fmod(nucleationSite.Z + params().pbcLength, params().pbcLength) * params().outputScalingFactor) << "): " <<
				(params().peierlsStress * s * 1e-6) << " MPa" << endl;

		// Assume a small/negligible activation barrier for the nucleation process to avoid failure of the kMC algorithm.
		activationEnergy = params().kpenergy_deltaH0 * 1e-12;

		// Use a kink separation corresponding to s==0.4, where it's at the minimum.
		sustainableKinkSeparation = params().kpwidth_l0 * (pow(0.4, -params().kpwidth_p) + params().kpwidth_w0) * pow(1.0 - 0.4, -params().kpwidth_q);

		return true;
	}

	// Calculate kink pair separation.
	sustainableKinkSeparation = params().kpwidth_l0 * (pow(s, -params().kpwidth_p) + params().kpwidth_w0) * pow(1.0 - s, -params().kpwidth_q);
	SIMULATION_ASSERT(sustainableKinkSeparation > 0);
	double extraq1 = 0;
	Point3 p1 = event.segment->node1()->pos();
	Point3 p2 = event.segment->node2()->pos();
    p1.Z = event.kinkPairPosition - params().kinkWidth / 2; 
	p2.Z = event.kinkPairPosition + params().kinkWidth / 2;
	
	if(pointDefects().isSoluteOnTheScrewDislocation(p1,p2))
		extraq1 = 0.25;

	activationEnergy = (params().kpenergy_deltaH0 + extraq1)* pow(1.0 - pow(s, params().kpenergy_p), params().kpenergy_q);//+0.25

	
	frequencyPrefactor = event.segment->getHSegmentLength() - sustainableKinkSeparation;
	if(frequencyPrefactor < 0) {
		if(s <= 0.4)
			frequencyPrefactor = 0;
		else
			frequencyPrefactor = 1;
	}

#if 0
	context().msgLogger(VERBOSITY_HIGH) << "Kink pair activation energy:" << endl;
	context().msgLogger(VERBOSITY_HIGH) << "  Kink direction: " << event.kinkDirection << endl;
	context().msgLogger(VERBOSITY_HIGH) << "  Segment: " << event.segment->node1()->pos() << " --- " << event.segment->node2()->pos() << endl;
	context().msgLogger(VERBOSITY_HIGH) << "  Resolved shear stress: " << resolvedShearStress << " [Pa]" << endl;
	context().msgLogger(VERBOSITY_HIGH) << "  Kink pair separation: " << sustainableKinkSeparation << " [b]" << endl;
	context().msgLogger(VERBOSITY_HIGH) << "  Delta H: " << activationEnergy << " eV" << endl;
#endif

	return true;
}

/******************************************************************************
* Executes the selected kink-pair nucleation event by modifying the current
* dislocation configuration.
*****************************************************************************/
bool Simulation::executeNucleationEvent(const KinkPairNucleationEvent& event)
{
	SIMULATION_ASSERT(event.kinkPairWidth > 0);

	if(network().insertKinkPair(event.segment, event.kinkPairPosition, event.kinkDirection, event.kinkPairWidth))
		return true; // Indicate that the dislocation configuration has changed.

	return false;	// No change was made to the current dislocation configuration.
}

bool Simulation::executePointDefectsEvent(const PointDefectsEventList& pevent)
{
	double x = pevent.x;
	double y = pevent.y;
	double z = pevent.z;
	double x1 = pevent.x1;
	double y1 = pevent.y1;
	double z1 = pevent.z1;
	int t = pevent.t;
	double x_line=pevent.x_line;
	double y_line=pevent.y_line;
	double x1_line=pevent.x1_line;
	double y1_line=pevent.y1_line;
	/*static const int sectorToLattice[8][3] = {
			{0,0,1},
			{0,0,-1},
			{1,0,0},
			{-1,0,0},
			{0,1,0},
			{0,-1,0},
			{1,1,1},
			{-1,-1,-1}
	};
	int x1 = x + sectorToLattice[pevent.n][0];
	int y1 = y + sectorToLattice[pevent.n][1];
	int z1 = z + sectorToLattice[pevent.n][2];*/
if(pointDefects().migration(pevent.p, x, y, z, x_line, y_line,x1, y1, z1,x1_line, y1_line,t))
	return true;

return false;
}

double Simulation::generatePointDefectsEventList(vector<PointDefectsEventList>& pevents)
{
	double totalRate = 0;
	_number = 0;
	double meanSquareMotion = 0;
	//double _number_of_atoms = 0;
	double kT = params().temperature * 8.6173324e-5;
	static const double sectorToLattice[6][4][5] = {
	{
			{-1/4.0,0,1/3.0,0,0},// 0->1
			{-1/4.0,1/4.0,-1/3.0,0,0},// 0->2
			{1/4.0,1/4.0,0,1/3.0,1/3.0},// 0->3
			{1/4.0,-1/2.0,0,1/3.0,-2/3.0}// 0->4
	},
	{
			{1/4.0,0,-1/3.0,0,0},//1->0
			{0,1/4.0,1/3.0,0,0},//1->2
			{1/4.0,-1/2.0,0,1/3.0,-2/3.0},//1->5
			{-1/2.0,1/4.0,0,-2/3.0,1/3.0}//1->3
	},
	{
			{1/4.0,-1/4.0,1/3.0,0,0},//2->0
			{0,-1/4.0,-1/3.0,0,0},//2->1
			{-1/2.0,1/4.0,0,-2/3.0,1/3.0},//2->4
			{1/4.0,1/4.0,0,1/3.0,1/3.0}//2->5
	},
	{
			{0,1/4.0,1/3.0,0,0},//3->4
			{-1/4.0,1/4.0,-1/3.0,0,0},//3->5
			{-1/4.0,-1/4.0,0,-1/3.0,-1/3.0},//3->0
			{1/2.0,-1/4.0,0,2/3.0,-1/3.0}//3->1
	},
	{
			{0,-1/4.0,-1/3.0,0,0},//4->3
			{-1/4.0,0,1/3.0,0,0},//4->5
			{-1/4.0,1/2.0,0,-1/3.0,2/3.0},//4->0
			{1/2.0,-1/4.0,0,2/3.0,-1/3.0}//4->2
	},
	{
			{1/4.0,-1/4.0,1/3.0,0,0},//5->3
			{1/4.0,0,-1/3.0,0,0},//5->4
			{-1/4.0,1/2.0,0,-1/3.0,2/3.0},//5->1
			{-1/4.0,-1/4.0,0,-1/3.0,-1/3.0}//5->2
	}};
	static const int sectorToLattice2[6][4] = {
			{1,2,3,4},
			{0,2,5,3},
			{0,1,4,5},
			{4,5,0,1},
			{3,5,0,2},
			{3,4,1,2}
	};
	for(auto row = pointDefects().defects().begin(); row != pointDefects().defects().end(); ++row) {
		SIMULATION_ASSERT(row != pointDefects().defects().end());
		double x_line =row->first.first;
		double y_line =row->first.second;
		//context().msgLogger(VERBOSITY_NORMAL) << " x_line: " << x_line << " y_line: " << y_line<< endl;
		PointDefect* head = row->second;
		for(PointDefect* pd = head; pd != nullptr; pd = pd->next) {
		//_number_of_atoms+=1;
		SIMULATION_ASSERT(pd != nullptr);
		double x = pd->x;
		double y = pd->y;
		double z = pd->position;
		int t = pd->p_type;
			for(int sector = 0; sector < 4; sector++) {
				bool isPointOccupied = false;
				double x1 = x + sectorToLattice[t][sector][0];
				double y1 = y + sectorToLattice[t][sector][1];
				double z1 = z + sectorToLattice[t][sector][2];
				double x1_line =x_line +sectorToLattice[t][sector][3];
				double y1_line =y_line +sectorToLattice[t][sector][4];
				int t1 = sectorToLattice2[t][sector];
				for(auto row2 = pointDefects().defects().begin(); row2 != pointDefects().defects().end(); ++row2) {
					PointDefect* head2 = row2->second;
					for(PointDefect* pt = head2; pt != nullptr; pt = pt->next) {
						double x2 = pt->x;
						double y2 = pt->y;
						double z2 = pt->position;
						if((x1==x2)&&(y1==y2)&&(z1==z2)) {
							isPointOccupied = true;
							break;
						}
					}
				
				if(!isPointOccupied){
				PointDefectsEventList pevent;
				double activationEnergy = 0;
				double frequencyPrefactor = 0;
				Point3 Pointp = pointDefects().getWorldPosition(x1, y1, z1);//modify because we need to bond with dislocation in lattice space
				pevent.p = pd;
				pevent.x = x;
				pevent.y = y;
				pevent.z = z;
				pevent.x1 = x1;
				pevent.y1 = y1;
				pevent.z1 = z1;
				pevent.t = t1;
				pevent.x_line = x_line;
				pevent.y_line = y_line;
				pevent.x1_line = x1_line;
				pevent.y1_line = y1_line;
				if(!calculateActivationEnergyPointDefects(pevent, activationEnergy, frequencyPrefactor,Pointp))
					continue;
				pevent.rate = frequencyPrefactor * exp(-activationEnergy / kT);
				//context().msgLogger(VERBOSITY_NORMAL) << "pevent.rate: " << pevent.rate << " pevent.x_line " << pevent.x_line<< "pevent.y_line " << pevent.y_line<< endl;
				SIMULATION_ASSERT(pevent.rate >= 0);
				totalRate += pevent.rate;
				pevents.push_back(pevent);
				}
				}
			}
		_number += 1;
		Point3 pOrginal = pd->p0;
		Point3 pNow = pd->p1;
		meanSquareMotion += (pNow.X - pOrginal.X)*(pNow.X - pOrginal.X) + (pNow.Y - pOrginal.Y)*(pNow.Y - pOrginal.Y) + (pNow.Z - pOrginal.Z)*(pNow.Z - pOrginal.Z);
		/// keep record of deminished pointdefects;
		_simulationPdefectMeanSquareMotion = (meanSquareMotion + _simulationPdefectMeanSquareMotionReduced)/_number;		
		}
	}
//context().msgLogger(VERBOSITY_NORMAL) << " _number_of_atoms: " << _number_of_atoms << endl;
return totalRate;
}
//<SegmentHandle>&?
double Simulation::generateDislocationBindingEventList(vector<DislocationBindingEventList>& devents){
	double totalRate = 0;
	double kT = params().temperature * 8.6173324e-5;
	for(SegmentIterator segment(network()); !segment.atEnd(); ++segment) {
		if(segment->isKink() == false){
			continue;
		}
		Point3 p1 = segment->node1()->pos();
		Point3 p2 = segment->node2()->pos();
		if (pointDefects().isSoluteOnTheDislocation(p1,p2))
		{
			DislocationBindingEventList devent1;//Segmenthandle devent = *segment;
			devent1.p = *segment;
			/*Point3 nucleationSite = params().unitCell * NodalPosition((segment->node1()->pos().X + segment->node2()->pos().X)/2, (segment->node1()->pos().Y + segment->node2()->pos().Y)/2, (segment->node1()->pos().Z + segment->node2()->pos().Z)/2);
			double resolvedShearStress;
			//resolvedShearStress = network().calculateResolvedShearStress(nucleationSite, segment->kinkDirection());
			SymmetricTensor2 stress = network().calculateLocalStress(nucleationSite);

			// MRSS plane direction.
			double theta_MRSS = atan2(-stress(0,2), stress(1,2));
			// MRSS value.
			double sigma_MRSS = sqrt(stress(0,2)*stress(0,2) + stress(1,2)*stress(1,2));

			// Kink plane direction.
			double theta_k = params().kinkDirections[segment->kinkDirection()].theta;

			// Compute signed angle between kink plane and MRSS plane.
			double chi = theta_MRSS - theta_k;

			// Compute local resolved shear stress in kink direction.
			resolvedShearStress = sigma_MRSS * (cos(chi) + params().nonschmid_a2 * cos(chi + deg2rad(60))) / params().nonschmid_a1;*/
			//double extraenergy;
			//extraenergy = 0.275 * params().stress_total/1e9;
			double extraenergy;
			extraenergy = fabs(0.275 * network().calculateLocalStressKink(segment->node1()->outSegment())/1e9);
			//context().msgLogger(VERBOSITY_NORMAL) << "extraenergy: " << extraenergy << endl;
			devent1.rate = 1.5 * 10e12 * exp(-(1.8 - extraenergy) / kT);
			totalRate += devent1.rate;
			devent1.t = 1;
			//totalRate += params().nucleationAttemptFrequency / params().numEventsPerSegment * exp(-0.175 / kT);
			devents.push_back(devent1);

			DislocationBindingEventList devent2;
			devent2.p = *segment;
			devent2.rate = 1.5 * 10e12 * exp(-(1.8 + extraenergy) / kT);
			totalRate += devent2.rate;
			devent2.t = -1;
			devents.push_back(devent2);
		}

	}
	return totalRate;
}

/*bool Simulation::executeBindingEvent(const SegmentHandle& devent)
{
	double z = 1;
	Point3 p1 = params().unitCell * devent->node1()->pos();
	Point3 p2 = params().unitCell * devent->node2()->pos();
	Point3 p = (p1 + p2)/2;
	SymmetricTensor2 stress = network().calculateLocalStress(p);
	if (stress(2,2) <= 0)
		z = -1;
	NodeHandle inCornerNode;
	NodeHandle outCornerNode;
	if (devent->inSegment()->isKink()){
		inCornerNode = _nodePool.construct(NodalPosition(devent->node1()->pos().X, devent->node1()->pos().Y, devent->node1()->pos().Z), 0);
		inCornerNode->setNext(devent->node1());
		inCornerNode->setPrev(devent->inSegment()->node1());
		devent->inSegment()->node1()->setNext(inCornerNode);
		devent->node1()->setPrev(inCornerNode);
	}
	if (devent->outSegment()->isKink()){
		outCornerNode = _nodePool.construct(NodalPosition(devent->node2()->pos().X, devent->node2()->pos().Y, devent->node2()->pos().Z), 0);
		outCornerNode->setNext(devent->node2());
		outCornerNode->setPrev(devent->node1());
		devent->node2()->setPrev(outCornerNode);
		devent->node1()->setNext(outCornerNode);
	}
	network().displaceKink(devent, z);

	 // Indicate that the dislocation configuration has changed.

	return true;	// No change was made to the current dislocation configuration.
}
*/

bool Simulation::calculateActivationEnergyPointDefects(const PointDefectsEventList& pevent, double& activationEnergy, double& frequencyPrefactor, const Point3& Pointp){
frequencyPrefactor = 6.4*10e12;
/*
double B = 309;
double U = 163;
double E = 415;
double v = 0.276;
double deltaV;
double deltaE;
double W1;
double W2;
SymmetricTensor2 stress = network().calculateLocalStress(Pointp);
SymmetricTensor2 strain;
strain(0,0) = stress(0,0) - v * stress(1,1) - v * stress(2,2);
strain(1,1) = stress(1,1) - v * stress(0,0) - v * stress(2,2);
strain(2,2) = stress(2,2) - v * stress(1,1) - v * stress(0,0);
strain(1,2) = 2 * (1 + v) * stress(1,2);
strain(2,1) = 2 * (1 + v) * stress(1,2);
strain(1,0) = 2 * (1 + v) * stress(0,1);
strain(0,1) = 2 * (1 + v) * stress(1,0);
strain(0,2) = 2 * (1 + v) * stress(0,2);
strain(2,0) = 2 * (1 + v) * stress(2,0);

double G = 2/U * (stress(0,1) + stress(0,2) + stress(1,2)) / 10e9;
SIMULATION_ASSERT(G <= 0.4);
double Er = 1/(3 * B) * (stress(0,0) + stress(1,1) + stress(2,2)) / 10e9;
SIMULATION_ASSERT(Er <= 0.4);
deltaV = - 2978 * G * G * G + 80 * G * G - 2.35 * G - 8.867 * Er * Er + 3.135 * Er + 1.0;
deltaE = - 4.77 * Er * Er - 0.219 * Er - 522 * G * G * G + 32.8 * G * G + 0.77 * G + 1.0;


activationEnergy = deltaV * 3.64 + deltaE * 1.54;

bool iskinkbool = false;
bool iskinkorg = false;
double activep = 0;
double activeorg = 0;

Point3 Pointorg = pointDefects().getWorldPosition(pevent.x, pevent.y, pevent.z);

 if (network().isSoluteOnTheDislocation(Pointp, iskinkbool))
 {
 	if (iskinkbool)
 		activep = 0.45;
 	else
 		activep = 0.25;
 }
if (network().isSoluteOnTheDislocation(Pointorg, iskinkorg))
 {
 	if (iskinkbool)+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 		activeorg = 0.45;
 	else
 		activeorg = 0.25;
 }

activationEnergy = activationEnergy + activeorg - activep;

*/
/*
double B = 309;
double U = 163;
double E = 415;
double v = 0.276;
double deltaV;
double deltaE;
double W1;
double W2;
SymmetricTensor2 stress = network().calculateLocalStress(Pointp);
SymmetricTensor2 strain;
strain(0,0) = stress(0,0) - v * stress(1,1) - v * stress(2,2);
strain(1,1) = stress(1,1) - v * stress(0,0) - v * stress(2,2);
strain(2,2) = stress(2,2) - v * stress(1,1) - v * stress(0,0);
strain(1,2) = 2 * (1 + v) * stress(1,2);
strain(2,1) = 2 * (1 + v) * stress(1,2);
strain(1,0) = 2 * (1 + v) * stress(0,1);
strain(0,1) = 2 * (1 + v) * stress(1,0);
strain(0,2) = 2 * (1 + v) * stress(0,2);
strain(2,0) = 2 * (1 + v) * stress(2,0);

double G = 2/U * (stress(0,1) + stress(0,2) + stress(1,2)) / 10e9;
SIMULATION_ASSERT(G <= 0.4);
double Er = 1/(3 * B) * (stress(0,0) + stress(1,1) + stress(2,2)) / 10e9;
SIMULATION_ASSERT(Er <= 0.4);
deltaV = - 2978 * G * G * G + 80 * G * G - 2.35 * G - 8.867 * Er * Er + 3.135 * Er + 1.0;
deltaE = - 4.77 * Er * Er - 0.219 * Er - 522 * G * G * G + 32.8 * G * G + 0.77 * G + 1.0;


// first test: activationEnergy = deltaV * 0.35 + deltaE * 0.40;
//activationEnergy = deltaV * 0.30 + deltaE * 0.35;
*/
activationEnergy = 0.2;
bool iskinkbool = false;
bool iskinkorg = false;
double activep = 0;
double activeorg = 0;

Point3 Pointorg = pointDefects().getWorldPosition(pevent.x_line, pevent.y_line, pevent.z);
/*
 if (network().isSoluteOnTheDislocation(Pointp, iskinkbool))
 {
 	if (iskinkbool)
 		activep = 0.25;
 	else
 		activep = 0.25;
 }
 */
/*
pevent.p = pd;
pevent.x = x;
pevent.y = y;
pevent.z = z;
pevent.x1 = x1;
pevent.y1 = y1;
pevent.z1 = z1;
pevent.t = t1;
pevent.x_line = x_line;
pevent.y_line = y_line;
pevent.x1_line = x1_line;
pevent.y1_line = y1_line;
*/
if (network().isSoluteOnTheDislocation(Pointorg, iskinkorg))
 {
 	/*
 	if (iskinkbool)
 		activeorg = 1.8;
 	else
 		activeorg = 1.8;*/
 	//context().msgLogger(VERBOSITY_NORMAL) << " whether bond: " << "YES" << endl;
 	activeorg = 1.8;


 }

activationEnergy = activationEnergy + activeorg;
//context().msgLogger(VERBOSITY_NORMAL) << " activationEnergy " << activationEnergy << endl;

//activationEnergy = 3.64 + 1.54;
//activationEnergy = 0.8;
return true;
}





















