#include "Dislocations.h"
#include "../simulation/Simulation.h"

using namespace std;

/*******************************************************************************
* Finds attracting cross-kinks that form a locked configuration.
********************************************************************************/
bool DislocationNetwork::detectLockedCrossKinks()
{
	if(params().crossKinkLockDistance == 0)
		return false;

	SIMULATION_ASSERT(params().enableLocalStress);

	bool modifiedNetwork = false;
	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		SegmentHandle segment = loop->segments();
		do {
			if(segment->isScrew() && segment->prevSegment()->isKink() && segment->nextSegment()->isKink()) {
				if(segment->getHSegmentLength() <= params().crossKinkLockDistance) {
					if(areCrossKinksLocked(segment)) {
						segment = lockCrossKink(segment, *loop);
						modifiedNetwork = true;
						continue;
					}
				}
			}
			segment = segment->nextSegment();
		}
		while(segment != loop->segments());
		eliminateSelfSegment(loop->segments());
	}

	return modifiedNetwork;
}

/*******************************************************************************
* Determines whether a cross-kink pair wants to stay locked or rather separate.
* This method takes the short screw segment that connects the two kinks as input.
********************************************************************************/
bool DislocationNetwork::areCrossKinksLocked(SegmentHandle screwSegment)
{
	SegmentHandle kink1 = screwSegment->prevSegment();
	SegmentHandle kink2 = screwSegment->nextSegment();

	SIMULATION_ASSERT(screwSegment->isScrew());
	SIMULATION_ASSERT(kink1->isKink() && kink2->isKink());

	// Before calculating the stresses, we need to separate the two kinks by extending the screw segment.
	double extensionLength = params().crossKinkLockDistance - screwSegment->getHSegmentLength();
	if(extensionLength > 0) {
		SegmentHandle nextScrew1 = kink1->prevSegment();
		SegmentHandle nextScrew2 = kink2->nextSegment();

		// If both kinks are already part of cross-kink locks, then they form a lock too.
		if(nextScrew1->isScrew() == false && nextScrew2->isScrew() == false)
			return true;

		double space1 = nextScrew1->isScrew() ? nextScrew1->getHSegmentLength() : 0.0;
		double space2 = nextScrew2->isScrew() ? nextScrew2->getHSegmentLength() : 0.0;
		if(space1 + space2 < extensionLength)
			return true;	// There is not enough space to unlock this cross-kink pair.

		double displacement1 = max(min(extensionLength * 0.5, space1), extensionLength - space2);
		double displacement2 = max(min(extensionLength * 0.5, space2), extensionLength - space1);
		SIMULATION_ASSERT(fabs(displacement1 + displacement2 - extensionLength) <= NODAL_POS_EPSILON);
		SIMULATION_ASSERT(displacement1 <= space1);
		SIMULATION_ASSERT(displacement2 <= space2);
		/*NodeHandle no1dep1 = kink1->node1();
		NodeHandle no1dep2 = kink1->node2();
		Point3 no1dep1p = no1dep1->pos();
		Point3 no1dep2p = no1dep2->pos();
		NodeHandle no2dep1 = kink2->node1();
		NodeHandle no2dep2 = kink2->node2();
		Point3 no2dep1p = no2dep1->pos();
		Point3 no2dep2p = no2dep2->pos();
		if(simulation().pointDefects().isSoluteOnTheDislocation(no2dep1p,no2dep2p))
		{
		displaceKink(kink2,  (displacement2 + displacement1));
		}
		else{
			if(simulation().pointDefects().isSoluteOnTheDislocation(no2dep1p,no2dep2p))
			{
			displaceKink(kink1, -(displacement1 + displacement2));
			}
			else{
			displaceKink(kink1, -displacement1);
			displaceKink(kink2,  displacement2);
			}
		}
		*/
		displaceKink(kink1, -displacement1);
		displaceKink(kink2,  displacement2);
		SIMULATION_ASSERT(fabs(screwSegment->getHSegmentLength() - params().crossKinkLockDistance) <= NODAL_POS_EPSILON);
	}

	// Compute resolved stress acting on kink 1.
	double stress1 = -calculateResolvedShearStress(getSegmentCenter(kink1), kink1->kinkDirection());

	// Compute resolved stress acting on kink 2.
	double stress2 = -calculateResolvedShearStress(getSegmentCenter(kink2), kink2->kinkDirection());

	// A positive difference in the forces means net attraction.
	return (stress1 - stress2) > 0;
}

/*******************************************************************************
* Collapses a short screw segment such that the adjacent kinks form a
* locked cross-kink.
********************************************************************************/
SegmentHandle DislocationNetwork::lockCrossKink(SegmentHandle screwSegment, DislocationLoop& loop)
{
	SegmentHandle kink1 = screwSegment->prevSegment();
	SegmentHandle kink2 = screwSegment->nextSegment();

	//context().msgLogger(VERBOSITY_HIGH) << "Locking cross kink at screw segment" << screwSegment->node1()->pos() << " - " << screwSegment->node2()->pos() << endl;

	SIMULATION_ASSERT(screwSegment->isScrew());
	SIMULATION_ASSERT(kink1->isKink() && kink2->isKink());

	double separation = screwSegment->getHSegmentLength();
	if(separation > 0) {
		SegmentHandle nextScrew1 = kink1->prevSegment();
		SegmentHandle nextScrew2 = kink2->nextSegment();

		if(nextScrew1->isScrew()) {
			displaceKink(kink1, +separation);
		}
		else if(nextScrew2->isScrew()) {
			displaceKink(kink2, -separation);
		}
		else {
			displaceKink(kink1, +separation * 0.5);
			displaceKink(kink2, -separation * 0.5);
		}
	}
	SIMULATION_ASSERT(fabs(screwSegment->getHSegmentLength()) <= NODAL_POS_EPSILON);

	// Finally remove degenerate screw segment from dislocation.
	if(kink1->kinkDirection() != params().kinkDirections[kink2->kinkDirection()].reverseDirection)
		return deleteSegment(screwSegment, loop);
	else {
		// Two opposite kinks will directly annihilate.
		deleteSegment(screwSegment, loop);
		deleteSegment(kink1, loop);
		return deleteSegment(kink2, loop);
	}
}

