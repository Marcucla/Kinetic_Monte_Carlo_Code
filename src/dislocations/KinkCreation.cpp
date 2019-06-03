#include "Dislocations.h"
#include "../simulation/Simulation.h"

using namespace std;

/*******************************************************************************
* Inserts a kink pair into a screw segment.
********************************************************************************/
bool DislocationNetwork::insertKinkPair(SegmentHandle segment, double kinkPairPos, int kinkDirectionIndex, double kinkPairWidth)
{
	SIMULATION_ASSERT(segment->isScrew());	// Segment must be an H-segment.
	SIMULATION_ASSERT(kinkPairWidth > 0);	// Width must be positive.
	SIMULATION_ASSERT(kinkPairWidth < params().lineLength);	// Kink separation cannot be larger than the simulation domain.

	NodeHandle node1 = segment->node1();
	NodeHandle node2 = segment->node2();

	// Compute position of left kink on screw segment.
	double leftPos = (kinkPairPos - 0.5 * kinkPairWidth);

	// Compute position of right kink on screw segment.
	double rightPos = leftPos + kinkPairWidth;

	// Special step for the initial straight dislocation line.
	if(segment->isStraightScrewLine()) {
		// Move the two existing nodes to the kink positions such that this segment can fully accommodate the kink pair.
		node1->setPos(NodalPosition(node1->pos().X, node1->pos().Y, leftPos));
		node2->setPos(NodalPosition(node2->pos().X, node2->pos().Y, rightPos));
		segment->setScrewLineVector(kinkPairWidth);
		segment->nextSegment()->setScrewLineVector(params().lineLength - kinkPairWidth);
	}

	double segmentWidth = segment->getHSegmentLength();
	if(segmentWidth <= 0)
		return false;	// Screw segment no longer exists.
	if(rightPos <= node1->pos().Z)
		return false;	// Kink pair is position outside the screw segment.
	if(leftPos >= node1->pos().Z + segmentWidth)
		return false;	// Kink pair is position outside the screw segment.

	// Clamp kink pair dimensions to screw segment.
	if(leftPos < node1->pos().Z) leftPos = node1->pos().Z;
	if(rightPos > node1->pos().Z + segmentWidth) rightPos = node1->pos().Z + segmentWidth;
	if(rightPos <= leftPos)
		return false;

	const NodalVector& kinkVector = params().kinkDirections[kinkDirectionIndex].h;

	NodeHandle firstCornerNode;
	if(!positionCompare(leftPos, node1->pos().Z)) {
		// Insert new corner node into H-segment.
		firstCornerNode = _nodePool.construct(NodalPosition(node1->pos().X, node1->pos().Y, leftPos), kinkVector, kinkDirectionIndex);
		firstCornerNode->setPrev(node1);
		node1->setNext(firstCornerNode);
		segment->setScrewLineVector(leftPos - node1->pos().Z);
	}
	else {
		firstCornerNode = node1;
		segment->setKinkLineVector(kinkDirectionIndex, kinkVector);
	}

	NodeHandle secondCornerNode;
	if(!positionCompare(rightPos, node2->pos().Z)) {
		// Insert new corner node into H-segment.
		secondCornerNode = _nodePool.construct(NodalPosition(node2->pos().X, node2->pos().Y, rightPos), node1->pos().Z + segmentWidth - rightPos);
		secondCornerNode->setNext(node2);
		node2->setPrev(secondCornerNode);
	}
	else secondCornerNode = node2;

	//context().msgLogger() << "Inserting kink pair in direction " << kinkDirectionIndex << " from " << firstCornerNode->pos() << " to " << secondCornerNode->pos() << endl;

	// Insert new segment between kinks.
	NodeHandle kinkNode1 = _nodePool.construct(firstCornerNode->pos() + kinkVector, rightPos - leftPos);
	NodeHandle kinkNode2 = _nodePool.construct(secondCornerNode->pos() + kinkVector, -kinkVector, params().kinkDirections[kinkDirectionIndex].reverseDirection);

	// Update doubly-linked list of nodes.
	SIMULATION_ASSERT(firstCornerNode->isKink());
	firstCornerNode->setNext(kinkNode1);
	kinkNode1->setPrev(firstCornerNode);
	kinkNode1->setNext(kinkNode2);
	kinkNode2->setPrev(kinkNode1);
	kinkNode2->setNext(secondCornerNode);
	secondCornerNode->setPrev(kinkNode2);

	return true;
}

bool DislocationNetwork::executeBindingEvent(DislocationBindingEventList& devent)
{
	double z = 1;
	//double D = (params().kinkDiffusivityCoefficient * params().temperature) * 1e20 /* Angstrom^2/m^2 */ / (params().blength * params().blength);
	/*Point3 p1 = params().unitCell * devent->node1()->pos();
	Point3 p2 = params().unitCell * devent->node2()->pos();
	Point3 p = (p1 + p2)/2;
	SymmetricTensor2 stress = calculateLocalStress(p);*/
	if ((calculateKinkVelocity(devent.p) <= 0)&& devent.t==1)
		z = -1;
	if((calculateKinkVelocity(devent.p)>= 0)&& devent.t==-1)
		z = -1;
	/*NodeHandle inCornerNode;
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
		outCornerNode->setNext(devent->nextSegment()->node1());
		outCornerNode->setPrev(devent->node2());
		devent->nextSegment()->node1()->setPrev(outCornerNode);
		devent->node2()->setNext(outCornerNode);
	}*/
	
	//double diffusional_displacement = sqrt(deltaTime * D);

	//if(diffusional_displacement >= params().pbcLength * 0.5 / params().blength)
	//context().error("Diffusional displacements are too large at this timestep and exceed half of the dislocation length. Delta t: %f, Delta x: %f, PBC length: %f", deltaTime, diffusional_displacement, params().pbcLength / params().blength);
	
	SegmentHandle segment;
	segment = devent.p;
	//segment->_kinkDiffusionDir = z;
	displaceKink(devent.p, z);
	//double deltaZ = z * diffusional_displacement;
	/*
	NodeHandle node1 = segment->node1();
	NodeHandle node2 = segment->node2();
	//context().msgLogger(VERBOSITY_HIGH) << "Moving node at " << node->pos() << " with velocity " << velocity << " by " << deltaZ << endl;

	node1->setPos(NodalPosition(node1->pos().X, node1->pos().Y, fmod(node1->pos().Z + deltaZ, params().lineLength)));
	// Adjust line vectors of preceding and following segment.
	SegmentHandle inSegment1 = node1->inSegment();
	SegmentHandle outSegment1 = node1->outSegment();
	if(inSegment1->isScrew()) {
	inSegment1->setScrewLineVector(inSegment1->getHSegmentLength() + deltaZ);
	}
	else {
	inSegment1->setKinkLineVector(inSegment1->kinkDirection(), inSegment1->lineVector() + NodalVector(0, 0, deltaZ));
	}
	if(outSegment1->isScrew()) {
	outSegment1->setScrewLineVector(outSegment1->getHSegmentLength() - deltaZ);
	}
	else {
	outSegment1->setKinkLineVector(outSegment1->kinkDirection(), outSegment1->lineVector() - NodalVector(0, 0, deltaZ));
	}


	
	node2->setPos(NodalPosition(node2->pos().X, node2->pos().Y, fmod(node2->pos().Z + deltaZ, params().lineLength)));
	// Adjust line vectors of preceding and following segment.
	SegmentHandle inSegment2 = node2->inSegment();
	SegmentHandle outSegment2 = node2->outSegment();
	if(inSegment2->isScrew()) {
	inSegment2->setScrewLineVector(inSegment2->getHSegmentLength() + deltaZ);
	}
	else {
	inSegment2->setKinkLineVector(inSegment2->kinkDirection(), inSegment2->lineVector() + NodalVector(0, 0, deltaZ));
	}
	if(outSegment2->isScrew()) {
	outSegment2->setScrewLineVector(outSegment2->getHSegmentLength() - deltaZ);
	}
	else {
	outSegment2->setKinkLineVector(outSegment2->kinkDirection(), outSegment2->lineVector() - NodalVector(0, 0, deltaZ));
	}
	*/
	 // Indicate that the dislocation configuration has changed.

	return true;	// No change was made to the current dislocation configuration.
}

