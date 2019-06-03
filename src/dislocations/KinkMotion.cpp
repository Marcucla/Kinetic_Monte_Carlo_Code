#include "Dislocations.h"
#include "../simulation/Simulation.h"

#include <iomanip>

using namespace std;

/*******************************************************************************
* Assigns velocities to the kink segments.
********************************************************************************/
void DislocationNetwork::calculateKinkVelocities(boost::mt19937& rng)
{
	//bool isKinkDiffusionActive = 0;
	bool isKinkDiffusionActive = params().kinkDiffusivityCoefficient != 0;
	boost::uniform_int<int> diffusionDirRNG(0,1);

	// Compute velocity of each kink segment.
	for(SegmentIterator segmentIterator(*this); !segmentIterator.atEnd(); ++segmentIterator) {
		segmentIterator->_kinkDiffusionDir = 0;
		if(segmentIterator->isKink()) {
			NodeHandle nodep1 = segmentIterator->node1();
			NodeHandle nodep2 = segmentIterator->node2();
			Point3 nodep1p = nodep1->pos();
			Point3 nodep2p = nodep2->pos();
			if (simulation().pointDefects().isSoluteOnTheDislocation(nodep1p,nodep2p)){
				segmentIterator->setKinkVelocity(0);
				//context().msgLogger(VERBOSITY_NORMAL) << "yes " << endl;
				segmentIterator->node1()->_nodeDiffusionDir = 0;
				segmentIterator->node2()->_nodeDiffusionDir = 0;
			}
			else{
				segmentIterator->setKinkVelocity(calculateKinkVelocity(*segmentIterator));
				//context().msgLogger(VERBOSITY_NORMAL) << segmentIterator->kinkVelocity() << endl;
				if(isKinkDiffusionActive)
				segmentIterator->_kinkDiffusionDir = diffusionDirRNG(rng) ? +1 : -1;
			}

			// Determine random diffusion direction.
			//if((isKinkDiffusionActive)&&(!simulation().pointDefects().isSoluteOnTheDislocation(nodep1p,nodep2p)))
				//segmentIterator->_kinkDiffusionDir = diffusionDirRNG(rng) ? +1 : -1;
		}
	}

	// Find groups of locked cross-kinks and let the move in sync with a net velocity.
	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		SegmentHandle headSegment = loop->segments();
		do {
			if(headSegment->isScrew()) break;
			headSegment = headSegment->nextSegment();
		}
		while(headSegment != loop->segments());

		SegmentHandle segment1 = headSegment;
		do {
			if(segment1->isKink()) {
				// Compute average velocity of all kinks belonging to one locked structure.
				bool isbindwithpointdefects = false;
				double averageVelocity = segment1->kinkVelocity();
				if(averageVelocity == 0)
					isbindwithpointdefects = true;
				int kinkCount = 1;
				SegmentHandle segment2 = segment1->nextSegment();
				while(segment2->isKink() && segment2 != segment1) {
					averageVelocity += segment2->kinkVelocity();
					if(segment2->kinkVelocity()== 0)
					isbindwithpointdefects = true;	
					kinkCount++;
					segment2 = segment2->nextSegment();
				}
				if(kinkCount > 1) {
					// Set the velocity of all kinks to the same average velocity.
					// Disable diffusional motion for cross kink segments.
					averageVelocity /= kinkCount;
					if(isbindwithpointdefects)
						averageVelocity = 0;
					SegmentHandle s = segment1;
					do {
						s->setKinkVelocity(averageVelocity);
						s->_kinkDiffusionDir = 0;
						s = s->nextSegment();
					}
					while(s != segment2);
				}
				if(segment2 == segment1) break;
				segment1 = segment2;
			}
			else segment1 = segment1->nextSegment();
		}
		while(segment1 != headSegment);
	}

	// Transfer kink segment velocities to adjacent nodes.
	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		NodeHandle node = loop->nodes();
		do {
			SegmentHandle inSegment = node->inSegment();
			SegmentHandle outSegment = node->outSegment();
			if(inSegment->isKink() == true && outSegment->isKink() == false) {
				node->setNodeVelocity(inSegment->kinkVelocity());
				node->_nodeDiffusionDir = inSegment->_kinkDiffusionDir;
			}
			else if(inSegment->isKink() == false && outSegment->isKink() == true) {
				node->setNodeVelocity(outSegment->kinkVelocity());
				node->_nodeDiffusionDir = outSegment->_kinkDiffusionDir;
			}
			else if(inSegment->isKink() == true && outSegment->isKink() == true) {
				// Check if two adjacent kink segments (=cross-kink) have the same velocity.
				SIMULATION_ASSERT(fabs(inSegment->kinkVelocity() - outSegment->kinkVelocity()) < NODAL_VELOCITY_EPSILON);
				node->setNodeVelocity(outSegment->kinkVelocity());
				node->_nodeDiffusionDir = outSegment->_kinkDiffusionDir;
			}
			else {
				node->setNodeVelocity(0);
				node->_nodeDiffusionDir = 0;
			}
			node = node->nextNode();
		}
		while(node != loop->nodes());
	}

#ifdef DEBUG_SIMULATION
	// Verify nodal velocities.
	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		SegmentHandle segment = loop->segments();
		do {
			if(segment->isKink()) {
				SIMULATION_ASSERT(segment->node1()->nodeVelocity() == segment->node2()->nodeVelocity());
				SIMULATION_ASSERT(segment->node1()->_nodeDiffusionDir == segment->node2()->_nodeDiffusionDir);
			}
			segment = segment->nextSegment();
		}
		while(segment != loop->segments());
	}
#endif
}

/*******************************************************************************
* Calculates the velocity from the force acting on a kink segment.
* The returned velocity if in units of b per second.
********************************************************************************/
double DislocationNetwork::calculateKinkVelocity(SegmentHandle kinkSegment)
{
	// Compute resolved stress acting on segment.
	double stress = -calculateResolvedShearStress(getSegmentCenter(kinkSegment), kinkSegment->kinkDirection());

	// Velocity is tau/drag (units [b/sec]).
	return stress / params().kinkDragCoefficient;
}
double DislocationNetwork::calculateLocalStressKink(SegmentHandle kinkSegment)
{
	double stress = -calculateResolvedShearStress(getSegmentCenter(kinkSegment), kinkSegment->kinkDirection());
	return stress;
}
/*******************************************************************************
* Determines the earliest time until two nodes will collide based on their
* current velocities.
********************************************************************************/
double DislocationNetwork::freeMotionTime()
{
	double D = (params().kinkDiffusivityCoefficient * params().temperature) * 1e20 /* Angstrom^2/m^2 */ / (params().blength * params().blength);

	// Nested loop over all pairs of nodes.
	double minTime = CAFLOAT_MAX;
	for(auto loop1 = dislocationLoops().begin(); loop1 != dislocationLoops().end(); ++loop1) {
		NodeHandle node1 = loop1->segments();
		do {
			for(auto loop2 = loop1; loop2 != dislocationLoops().end(); ++loop2) {
				NodeHandle node2 = loop2->nodes();
				if(loop2 == loop1) {
					node2 = node1->nextNode();
					if(node2 == loop1->nodes()) continue;
				}
				do {
					// Is second node on same atomic row as first node?
					// Do not care about nodes which are already on top of each other.
					if(node2->pos().X == node1->pos().X && node2->pos().Y == node1->pos().Y && positionCompare(node2->pos().Z, node1->pos().Z) == false) {
						SIMULATION_ASSERT(node2 != node1);

						// Compute collision time.
						double collisionTime = CAFLOAT_MAX;
						double relativeVelocity = node2->nodeVelocity() - node1->nodeVelocity();
						if(relativeVelocity != 0) {
							double relativePosition = node2->pos().Z - node1->pos().Z;
							// Take into account periodic boundary.
							while(relativePosition >=  params().lineLength) relativePosition -= params().lineLength;
							while(relativePosition < 0) relativePosition += params().lineLength;
							SIMULATION_ASSERT(relativePosition >= 0 && relativePosition < params().pbcLength);

							if(D == 0 || node1->_nodeDiffusionDir - node2->_nodeDiffusionDir == 0) {
								if(relativeVelocity > 0) {
									collisionTime = ((double)params().lineLength - relativePosition) / relativeVelocity;
								}
								else if(relativeVelocity < 0) {
									collisionTime = relativePosition / -relativeVelocity;
								}
								if(collisionTime >= 0 && collisionTime < minTime) minTime = collisionTime;
							}
							else {
								double relativeDiffusionDir = node2->_nodeDiffusionDir - node1->_nodeDiffusionDir;
								double asq = relativeDiffusionDir * relativeDiffusionDir;
								double denom = 1.0 * 2.0 * relativeVelocity * relativeVelocity;
								double term2 = -4.0 * relativeVelocity * relativePosition + D * asq;
								if(term2 >= 0) {
									double term1 = D * asq - 2.0 * relativeVelocity * relativePosition;
									double term3 = relativeDiffusionDir * sqrt(D*term2);
									collisionTime = (term1 + term3) / denom;
									if(collisionTime >= 0 && collisionTime < minTime) minTime = collisionTime;
									collisionTime = (term1 - term3) / denom;
									if(collisionTime >= 0 && collisionTime < minTime) minTime = collisionTime;
								}
								relativePosition = (double)params().lineLength - relativePosition;
								term2 = 4.0 * relativeVelocity * relativePosition + D * asq;
								if(term2 >= 0) {
									double term1 = D * asq + 2.0 * relativeVelocity * relativePosition;
									double term3 = relativeDiffusionDir * sqrt(D*term2);
									collisionTime = (term1 + term3) / denom;
									if(collisionTime >= 0 && collisionTime < minTime) minTime = collisionTime;
									collisionTime = (term1 - term3) / denom;
									if(collisionTime >= 0 && collisionTime < minTime) minTime = collisionTime;
								}
							}
						}
					}
					node2 = node2->nextNode();
				}
				while(node2 != loop2->nodes());
			}
			node1 = node1->nextNode();
		}
		while(node1 != loop1->nodes());
	}

	return minTime;
}

/*******************************************************************************
* Determines the time the fastest kink segment needs to bind solute atoms
* maximum distance (in Burgers vector units).
********************************************************************************/
double DislocationNetwork::minTimeForBind(double t_mig)
{
	double D = (params().kinkDiffusivityCoefficient * params().temperature) * 1e20 /* Angstrom^2/m^2 */ / (params().blength * params().blength);

	double minTime = t_mig;
	//context().msgLogger(VERBOSITY_NORMAL)  << "CAFLOAT_MAX" << CAFLOAT_MAX  << endl;
	double timep1 = minTime;
	NodeHandle tempp1;
	NodeHandle tempp2;
	double movement1Z;
	for(SegmentIterator segmentIterator(*this); !segmentIterator.atEnd(); ++segmentIterator) 
	{
		if(segmentIterator->isKink() == false)
			continue;
		NodeHandle nodep1 = segmentIterator->node1();
		NodeHandle nodep2 = segmentIterator->node2();
		double nodep1velocity = nodep1->nodeVelocity();
		double nodep2velocity = nodep2->nodeVelocity();
		if (nodep1velocity == 0)
			continue;
		if (nodep2velocity == 0)
			continue;

			// Add diffusional displacement.

		double nodep1displacement = nodep1velocity * timep1;

		double diffusional_displacement = sqrt(timep1 * D);
		
		//context().msgLogger(VERBOSITY_NORMAL)  << "nodep1velocity" << nodep1velocity  << endl;
		
		//context().msgLogger(VERBOSITY_NORMAL)  << "timep1" << timep1  << endl;
		//context().msgLogger(VERBOSITY_NORMAL)  << "nodep1displacement" << nodep1displacement  << endl;
		//context().msgLogger(VERBOSITY_NORMAL)  << "diffusional_displacement" << diffusional_displacement  << endl;
		//SIMULATION_ASSERT(nodep1->_nodeDiffusionDir != 0);
		
		if(nodep1->_nodeDiffusionDir == 1)
			nodep1displacement += diffusional_displacement;
		else if(nodep1->_nodeDiffusionDir == -1)
			nodep1displacement -= diffusional_displacement;
		

		Point3 nodep1p = nodep1->pos();
		Point3 nodep2p = nodep2->pos();
		double temp1Z = simulation().pointDefects().bindPointDefects(nodep1p,nodep2p,nodep1displacement);

		double sign = 0;

		if(nodep1->_nodeDiffusionDir * nodep1velocity > 0)
			sign = 1;
		else if (nodep1->_nodeDiffusionDir * nodep1velocity < 0)
			sign = -1;

		//SIMULATION_ASSERT(sign != 0);

		if(fabs(temp1Z)<fabs(nodep1displacement))
		{
			double t;
			double vk = fabs(nodep1velocity);
			double distance = fabs(temp1Z);
			if(nodep1->_nodeDiffusionDir != 0)
				t = fabs(D + 2.0 * vk * distance - sign * sqrt(D*(4.0 * vk * distance + D))) / (2.0 * vk * vk);
			else
				t = distance / vk;

/*		timep1 = fabs(temp1Z)/fabs(nodep1displacement) * timep1;
		tempp1 = nodep1;
		tempp2 = nodep2;
		movement1Z = temp1Z;
*/
			//SIMULATION_ASSERT(t < timep1);
			if (t < timep1)
			timep1 = t;

			// test
		/*
		double test = nodep1velocity * timep1;

		double test2 = sqrt(timep1 * D);

		if(nodep1->_nodeDiffusionDir == 1)
			test += test2;
		else if(nodep1->_nodeDiffusionDir == -1)
			test -= test2;

		SIMULATION_ASSERT(fabs(test - temp1Z) < 0.01);
		*/
		}
	}
/*	if((timep1 < minTime)&&(timep1 < maxMotionTime))
	{
	minTime = timep1;
	SegmentHandle kink = tempp1->outSegment();
	SegmentHandle segmentt = kink;
	kink->setKinkVelocity(0);
	double movement2Z = -movement1Z;
	displaceKink(kink,movement2Z);
	tempp1->setNodeVelocity(0);
	tempp2->setNodeVelocity(0);
	tempp1->_nodeDiffusionDir = 0;
	tempp2->_nodeDiffusionDir = 0;
	kink->_kinkDiffusionDir = 0;
	do {
		segmentt->setKinkVelocity(0);
		displaceKink(segmentt,movement1Z);
		segmentt->_kinkDiffusionDir = 0;
		segmentt->node1()->setNodeVelocity(0);
		segmentt->node1()->_nodeDiffusionDir = 0;
		segmentt->node2()->setNodeVelocity(0);
		segmentt->node2()->_nodeDiffusionDir = 0;
		segmentt = segmentt->nextSegment();
	}
	while(segmentt->isKink());
	
	segmentt = kink;
	do {
		segmentt->setKinkVelocity(0);
		displaceKink(segmentt,movement1Z);
		segmentt->_kinkDiffusionDir = 0;
		segmentt->node1()->setNodeVelocity(0);
		segmentt->node1()->_nodeDiffusionDir = 0;
		segmentt->node2()->setNodeVelocity(0);
		segmentt->node2()->_nodeDiffusionDir = 0;
		segmentt = segmentt->prevSegment();
	}
	while(segmentt->isKink());

}

*/
	minTime = timep1;
	//context().msgLogger(VERBOSITY_NORMAL)  << "minTime " << minTime << endl;
	return minTime;

}

/*******************************************************************************
* Determines the time the fastest kink segment needs to travel the given
* maximum distance (in Burgers vector units).
********************************************************************************/

double DislocationNetwork::minTimeForDistance(double distance)
{
	SIMULATION_ASSERT(distance >= 0);

	double minTime = CAFLOAT_MAX;
	double D = (params().kinkDiffusivityCoefficient * params().temperature) * 1e20 /* Angstrom^2/m^2 */ / (params().blength * params().blength);

	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		NodeHandle node = loop->nodes();
		do {
			if(node->nodeVelocity() != 0) {
				double vk = fabs(node->nodeVelocity());
				double t;
				if(node->_nodeDiffusionDir != 0)
					t = fabs(D + 2.0 * vk * distance - sqrt(D*(4.0 * vk * distance + D))) / (2.0 * vk * vk);
				else
					t = distance / vk;
				if(t < minTime) minTime = t;
			}
			node = node->nextNode();
		}
		while(node != loop->nodes());
	}
	return minTime;
}

/*******************************************************************************
* Advances the simulation by the given amount and lets kinks move.
********************************************************************************/
void DislocationNetwork::propagateKinks(double delta_t)
{
	double D = (params().kinkDiffusivityCoefficient * params().temperature) * 1e20 /* Angstrom^2/m^2 */ / (params().blength * params().blength);
	double diffusional_displacement = sqrt(delta_t * D);

	if(diffusional_displacement >= params().pbcLength * 0.5 / params().blength)
		context().error("Diffusional displacements are too large at this timestep and exceed half of the dislocation length. Delta t: %f, Delta x: %f, PBC length: %f", delta_t, diffusional_displacement, params().pbcLength / params().blength);

	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		NodeHandle node = loop->nodes();
		do {
			double velocity = node->nodeVelocity();

			// Shift node along atomic row.
			double deltaZ = velocity * delta_t;

			// Add diffusional displacement.
			if(diffusional_displacement != 0) {
				if(node->_nodeDiffusionDir == 1)
					deltaZ += diffusional_displacement;
				else if(node->_nodeDiffusionDir == -1)
					deltaZ -= diffusional_displacement;
			}

			if(deltaZ != 0) {
				//context().msgLogger(VERBOSITY_HIGH) << "Moving node at " << node->pos() << " with velocity " << velocity << " by " << deltaZ << endl;

				node->setPos(NodalPosition(node->pos().X, node->pos().Y, fmod(node->pos().Z + deltaZ, params().lineLength)));

				// Adjust line vectors of preceding and following segment.
				SegmentHandle inSegment = node->inSegment();
				SegmentHandle outSegment = node->outSegment();
				if(inSegment->isScrew()) {
					inSegment->setScrewLineVector(inSegment->getHSegmentLength() + deltaZ);
				}
				else {
					inSegment->setKinkLineVector(inSegment->kinkDirection(), inSegment->lineVector() + NodalVector(0, 0, deltaZ));
				}
				if(outSegment->isScrew()) {
					outSegment->setScrewLineVector(outSegment->getHSegmentLength() - deltaZ);
				}
				else {
					outSegment->setKinkLineVector(outSegment->kinkDirection(), outSegment->lineVector() - NodalVector(0, 0, deltaZ));
				}
			}
			node = node->nextNode();
		}
		while(node != loop->nodes());
	}

	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		NodeHandle segment = loop->segments();
		do {
			if(segment->isScrew()) {
				if(fabs(segment->getHSegmentLength()) > (double)params().lineLength + NODAL_POS_EPSILON)
					context().error("Length of screw segment exceeded box size (periodicity length): %f", segment->getHSegmentLength());
			}
			segment = segment->nextSegment();
		}
		while(segment != loop->segments());
	}
}

/*******************************************************************************
* Moves a kink segment along the screw direction by the given amount (in units of b).
* The caller of this method has to make sure that the displacement will not lead to
* an invalid dislocation configuration.
********************************************************************************/
void DislocationNetwork::displaceKink(SegmentHandle kink, double deltaZ)
{
	SIMULATION_ASSERT(kink->isKink());

	if(deltaZ == 0)
		return;	// Nothing to do.

	SegmentHandle segment;

	// Shift nodes (including those of locked cross-kinks).
	segment = kink;
	do {
		NodeHandle node = segment->node2();
		node->setPos(NodalPosition(node->pos().X, node->pos().Y, fmod(node->pos().Z + deltaZ, params().lineLength)));
		segment = segment->nextSegment();
	}
	while(segment->isKink());

	// Adjust line vector of following segment.
	SIMULATION_ASSERT(segment->isScrew());
	segment->setScrewLineVector(segment->getHSegmentLength() - deltaZ);

	// Shift nodes (including those of locked cross-kinks).
	segment = kink;
	do {
		NodeHandle node = segment->node1();
		node->setPos(NodalPosition(node->pos().X, node->pos().Y, fmod(node->pos().Z + deltaZ, params().lineLength)));
		segment = segment->prevSegment();
	}
	while(segment->isKink());

	// Adjust line vector of preceding segment.
	SIMULATION_ASSERT(segment->isScrew());
	segment->setScrewLineVector(segment->getHSegmentLength() + deltaZ);
}

bool DislocationNetwork::isSoluteOnTheDislocation(const Point3& p, bool& iskinkbool)
{		
		bool Is = false;
		double small=0;
		double large=0;
		Point3 p4 = params().inverseUnitCell * p;
		for(SegmentIterator segmentIterator(*this); !segmentIterator.atEnd(); ++segmentIterator) {
		//Vector3 l = getSegmentVectorWorld(*segmentIterator);
		Point3 p1 = segmentIterator->node1()->pos();
		Point3 p2 = segmentIterator->node2()->pos();

		//context().msgLogger(VERBOSITY_NORMAL)  << "p1:x " << p1.X << "p1:y " << p1.Y<< "p1:z " << p1.Z << "p2:z " << p2.Z << "p4:x " << p4.X << "p4:y " << p4.Y<< "p1:z " << p4.Z<< endl;
		if(segmentIterator->isKink()){
			//if((fabs(p1.X-p4.X)<=0.001)&&(fabs(p1.Y-p4.Y)<=0.001))
				//Is = true;
			//if((fabs(p2.X-p4.X)<=0.001)&&(fabs(p2.Y-p4.Y)<=0.001))
				//Is = true;

			if((fabs(p1.X-p4.X)<=0.1)&&(fabs(p1.Y-p4.Y)<=0.1)&&((fabs(p1.Z-p4.Z)<=0.001)||(fabs(p1.Z-p4.Z+params().lineLength)<=0.001)||(fabs(p1.Z-p4.Z-params().lineLength)<=0.001)))
				Is = true;
			if((fabs(p2.X-p4.X)<=0.1)&&(fabs(p2.Y-p4.Y)<=0.1)&&((fabs(p2.Z-p4.Z)<=0.001)||(fabs(p2.Z-p4.Z+params().lineLength)<=0.001)||(fabs(p2.Z-p4.Z-params().lineLength)<=0.001)))
				Is = true;
		}
		else{
			if(p1.Z>p2.Z){
			if((fabs(p1.X - p4.X)<=0.1)&&(fabs(p1.Y - p4.Y)<=0.1)&&(((p1.Z <= p4.Z)&&((p2.Z+params().lineLength) >= p4.Z))||((p1.Z <= (p4.Z + params().lineLength))&&((p2.Z+params().lineLength) >= (p4.Z + params().lineLength)))||((p1.Z <= (p4.Z - params().lineLength))&&((p2.Z+params().lineLength) >= (p4.Z - params().lineLength)))))
				Is = true;				
			}
			else{
			if((fabs(p1.X - p4.X)<=0.1)&&(fabs(p1.Y - p4.Y)<=0.1)&&(((p1.Z <= p4.Z)&&(p2.Z >= p4.Z))||((p1.Z <= (p4.Z + params().lineLength))&&(p2.Z >= (p4.Z + params().lineLength)))||((p1.Z <= (p4.Z - params().lineLength))&&(p2.Z >= (p4.Z - params().lineLength)))))
				Is = true;
			}

			//if((fabs(p1.X - p4.X)<=0.001)&&(fabs(p1.Y - p4.Y)<=0.001))
			//if((fabs(p1.X - p4.X)<=0.001)&&(fabs(p1.Y - p4.Y)<=0.001)&&(((p1.Z <= p4.Z)&&(p2.Z >= p4.Z))||((p1.Z <= (p4.Z + params().lineLength))&&(p2.Z >= (p4.Z + params().lineLength)))||((p1.Z <= (p4.Z - params().lineLength))&&(p2.Z >= (p4.Z - params().lineLength)))))
				//Is = true;
		} 
		//Vector3 l1 = p - pd;
		/*
		if(fabs(sqrt((l.X * l.X +l.Y * l.Y + l.Z * l.Z)*(l1.X * l1.X +l1.Y * l1.Y + l1.Z * l1.Z)) - (l.X * l1.X + l.Y * l1.Y + l.Z * l1.Z)) <= 0.01)
			{
				if((l.X * l.X +l.Y * l.Y + l.Z * l.Z) >= (l1.X * l1.X +l1.Y * l1.Y + l1.Z * l1.Z))
				{
					Is = true;
					if(segmentIterator->isKink())
					iskinkbool = true;
			}
			}
		*/

	}
	if (Is==true){
		//context().msgLogger(VERBOSITY_NORMAL)  << "Is: " << Is << endl;
	}
		return Is;
}


