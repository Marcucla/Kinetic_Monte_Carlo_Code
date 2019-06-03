#include "Dislocations.h"
#include "../simulation/Simulation.h"

using namespace std;

/*******************************************************************************
* Removes redundant segments and nodes and makes topological changes at
* multi-nodes if indicated by the curvature criterion.
*
* Also removes annihilated segments and kills degenerate loops completely.
********************************************************************************/
void DislocationNetwork::cleanupNetwork()
{
	context().msgLogger(VERBOSITY_HIGH) << "Cleaning up dislocation network." << endl;

	// Find nodes which are located right on top of a screw segment.
	// Then insert an extra node into the screw segment to create
	// a proper multi-node (consisting of two or more 2-nodes at the same location).
	for(auto loop1 = dislocationLoops().begin(); loop1 != dislocationLoops().end(); ++loop1) {
		SegmentHandle segment = loop1->segments();
		do {
			if(segment->isScrew()) {
				SIMULATION_ASSERT(fabs(segment->getHSegmentLength()) - NODAL_POS_EPSILON <= (double)params().lineLength);
				for(auto loop2 = dislocationLoops().begin(); loop2 != dislocationLoops().end(); ++loop2) {
					NodeHandle node = loop2->nodes();
					do {
						NodeHandle startNode = segment->node1();
						NodeHandle endNode = segment->node2();
						if(node->pos().X == startNode->pos().X && node->pos().Y == startNode->pos().Y && node != startNode && node != endNode) {
							double deltaZ = node->pos().Z - startNode->pos().Z;
							bool isInside = false;
							if(segment->getHSegmentLength() > 0) {
								while(deltaZ > +params().lineLength) deltaZ -= params().lineLength;
								while(deltaZ < 0) deltaZ += params().lineLength;
								isInside = (deltaZ > NODAL_POS_EPSILON) && (deltaZ < segment->getHSegmentLength() - NODAL_POS_EPSILON);
							}
							else if(segment->getHSegmentLength() < 0) {
								while(deltaZ > 0) deltaZ -= params().lineLength;
								while(deltaZ < -params().lineLength) deltaZ += params().lineLength;
								isInside = (deltaZ < -NODAL_POS_EPSILON) && (deltaZ > segment->getHSegmentLength() + NODAL_POS_EPSILON);
							}
							if(isInside) {
								// Insert extra node into screw segment.
								NodeHandle newNode = _nodePool.construct(startNode->pos() + NodalVector(0,0,deltaZ), segment->getHSegmentLength() - deltaZ);
								context().msgLogger(VERBOSITY_HIGH) << "Inserting extra node into screw segment [" << startNode->pos() << "," << endNode->pos() << "] at " << newNode->pos() << "  deltaZ=" << deltaZ << " linevec=" << segment->lineVector() << endl;
								SIMULATION_ASSERT(positionCompare(newNode->pos(), node->pos()));
								newNode->setPrev(startNode);
								newNode->setNext(endNode);
								startNode->setNext(newNode);
								endNode->setPrev(newNode);
								startNode->setScrewLineVector(deltaZ);
							}
						}
						node = node->nextNode();
					}
					while(node != loop2->nodes());
				}
			}
			segment = segment->nextSegment();
		}
		while(segment != loop1->segments());
	}

	// Remove zero length screw segments.
	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		SegmentHandle segment = loop->segments();
		do {
			if(segment->isScrew() && fabs(segment->getHSegmentLength()) <= NODAL_POS_EPSILON) {
				SIMULATION_ASSERT(positionCompare(segment->node1()->pos(), segment->node2()->pos()));
				context().msgLogger(VERBOSITY_HIGH) << "Removed zero length H-segment " << segment->node1()->pos() << (params().unitCell * segment->node1()->pos()) << endl;
				segment = deleteSegment(segment, *loop);
			}
			else segment = segment->nextSegment();
		}
		while(segment != loop->segments());
	}

	// Doubly-nested loop over all nodes.
	// Find pairs of nodes located on the same location and reconnect the four incident segments based on curvature criterion.
	for(auto loop1 = dislocationLoops().begin(); loop1 != dislocationLoops().end(); ++loop1) {
		NodeHandle node1 = loop1->nodes();
		do {
			for(auto loop2 = loop1; loop2 != dislocationLoops().end(); ++loop2) {
				NodeHandle node2 = loop2->nodes();
				if(loop2 == loop1) {
					node2 = node1->nextNode();
					if(node2 == loop1->nodes()) continue;
				}
				do {
					SIMULATION_ASSERT(node1 != node2);
					if(positionCompare(node1->pos(), node2->pos())) {
						if(shouldFlipJunction(node1, node2)) {

							// Re-connect segments.
							NodeHandle prev1 = node1->prevNode();
							NodeHandle prev2 = node2->prevNode();
							SIMULATION_ASSERT(prev1 != prev2);
							prev1->setNext(node2);
							node2->setPrev(prev1);
							prev2->setNext(node1);
							node1->setPrev(prev2);
							NodalPosition oldPos = node1->pos();
							node1->setPos(node2->pos());
							node2->setPos(oldPos);
							if(loop2 == loop1) {
								// Create two separate loops if the 4-junction was part of a single loop.
								context().msgLogger(VERBOSITY_HIGH) << "Flipping 4-junction at " << node1->pos() << (params().unitCell * node1->pos()) << " -> Creating new loop." << endl;
								if(node1 == loop1->nodes()) {
									eliminateSelfSegment(node2);
									loop1->setHead(node2);
								}
								createLoop(node1);
								SIMULATION_ASSERT(node2 != node1->nextNode());
								node1 = node2;
							}
							else {
								// Coalesce loops if the 4-junction belonged to two separate loops.
								context().msgLogger(VERBOSITY_HIGH) << "Flipping 4-junction at " << node1->pos() << (params().unitCell * node1->pos()) << " -> Joining loops." << endl;
								auto prevLoop = loop2; --prevLoop;
								dislocationLoops().erase(loop2);
								loop2 = prevLoop;
								node1 = node2;
								break;
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

	// Remove zero length screw segments.
	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		SegmentHandle segment = loop->segments();
		do {
			if(segment->isScrew() && fabs(segment->getHSegmentLength()) <= NODAL_POS_EPSILON && segment->nextSegment() != segment) {
				SIMULATION_ASSERT(positionCompare(segment->node1()->pos(), segment->node2()->pos()));
				context().msgLogger(VERBOSITY_HIGH) << "Removed zero length H-segment " << segment->node1()->pos() << (params().unitCell * segment->node1()->pos()) << endl;
				segment = deleteSegment(segment, *loop);
			}
			else segment = segment->nextSegment();
		}
		while(segment != loop->segments());
	}

	// Remove redundant nodes on H-segments, i.e., nodes which are adjacent to two screw segments.
	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		SegmentHandle segment = loop->segments();
		do {
			// Find two consecutive H-segments. Remove shared node and coalesce them into a single segment.
			if(segment->isScrew()) {
				SegmentHandle nextSegment = segment->nextSegment();
				SIMULATION_ASSERT(nextSegment != segment);
				if(nextSegment->isScrew() && segment->getHSegmentLength() * nextSegment->getHSegmentLength() >= 0 && nextSegment->nextSegment() != segment) {
					NodeHandle nodeToRemove = nextSegment->node1();
					context().msgLogger(VERBOSITY_HIGH) << "Removing extra H-segment node at " << nodeToRemove->pos() << (params().unitCell * nodeToRemove->pos()) << endl;

					SIMULATION_ASSERT(segment->node1()->pos().X == nodeToRemove->pos().X && segment->node1()->pos().Y == nodeToRemove->pos().Y);
					SIMULATION_ASSERT(nextSegment->node2()->pos().X == nodeToRemove->pos().X && nextSegment->node2()->pos().Y == nodeToRemove->pos().Y);
					SIMULATION_ASSERT(nextSegment == nodeToRemove);

					// Extend first screw segment to cover the second screw segment.
					segment->setScrewLineVector(segment->getHSegmentLength() + nextSegment->getHSegmentLength());

					// Delete second screw segment.
					deleteSegment(nextSegment, *loop);

					SIMULATION_ASSERT(positionCompare(segment->node2()->pos().Z - segment->node1()->pos().Z, segment->lineVector().Z));
					if(fabs(segment->getHSegmentLength()) > (double)params().lineLength + NODAL_POS_EPSILON)
						context().error("Length of merged screw segment exceeds box size (periodicity length): %f [b]", segment->getHSegmentLength());
				}
			}
			segment = segment->nextSegment();
		}
		while(segment != loop->segments());
	}

	// Kill degenerate loops that consist of only two segments.
	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ) {
		SegmentHandle segment1 = loop->segments();
		SegmentHandle segment2 = segment1->nextSegment();
		if(segment2->nextSegment() == segment1 &&
				segment1->lineVector().X + segment2->lineVector().X == 0 &&
				segment1->lineVector().Y + segment2->lineVector().Y == 0 &&
				fabs(segment1->lineVector().Z + segment2->lineVector().Z) <= NODAL_POS_EPSILON) {
			_nodePool.destroy(segment1);
			if(segment2 != segment1)
				_nodePool.destroy(segment2);
			loop = dislocationLoops().erase(loop);
		}
		else ++loop;
	}
}

/*******************************************************************************
* Deletes the given segment from the given loop.
* Returns the segment that follows the deleted segment in the loop.
********************************************************************************/
SegmentHandle DislocationNetwork::deleteSegment(SegmentHandle segment, DislocationLoop& loop)
{
	SegmentHandle nextSegment = segment->nextSegment();
	segment->node1()->prevNode()->setNext(segment->node2());
	segment->node2()->setPrev(segment->node1()->prevNode());
	if(segment == loop.segments()) {
		eliminateSelfSegment(nextSegment);
		SIMULATION_ASSERT(nextSegment->nextSegment() != nextSegment);
		loop.setHead(nextSegment);
	}
	_nodePool.destroy(segment);

	return nextSegment;
}

/*******************************************************************************
* Determines whether the given 4-junction should be flipped based on a curvature
* criterion.
********************************************************************************/
bool DislocationNetwork::shouldFlipJunction(NodeHandle node1, NodeHandle node2)
{
	Vector3 seg1in  = getSegmentVectorWorld(node1->inSegment());
	Vector3 seg1out = getSegmentVectorWorld(node1->outSegment());
	Vector3 seg2in  = getSegmentVectorWorld(node2->inSegment());
	Vector3 seg2out = getSegmentVectorWorld(node2->outSegment());
	double dot1in  = DotProduct(seg1in, seg1out) / Length(seg1in) / Length(seg1out);
	double dot1out = DotProduct(seg2in, seg2out) / Length(seg2in) / Length(seg2out);
	double dot2in  = DotProduct(seg1in, seg2out) / Length(seg1in) / Length(seg2out);
	double dot2out = DotProduct(seg2in, seg1out) / Length(seg2in) / Length(seg1out);

	// Do not flip if existing segments would annihilate.
	if(dot1in < -1.0 + CAFLOAT_EPSILON || dot1out < -1.0 + CAFLOAT_EPSILON)
		return false;

	// Flip junction if resulting segments would annihilate.
	if(dot2in < -1.0 + CAFLOAT_EPSILON || dot2out < -1.0 + CAFLOAT_EPSILON)
		return true;

	// Otherwise use curvature-based criterion.
	double dot1 = dot1in + dot1out;
	double dot2 = dot2in + dot2out;

	if(dot2 > dot1 + CAFLOAT_EPSILON)
		return true;
	return false;
}

/*******************************************************************************
* Checks whether the data structure is in good order. This is for debugging only.
*
* We verify the following conditions:
*
*  - PREV and NEXT pointers in doubly-linked node/segment list must be consistent.
*  - Kink dislocation segments must have a line vector that is a valid kink direction.
*  - Screw dislocation segments must be parallel to the Z axis.
*  - Line vectors of segments must be consistent with the absolute coordinates of
*    the two adjacent nodes.
*  - The total Burgers vector content must be preserved (one periodic screw dislocation).
********************************************************************************/
void DislocationNetwork::validate(bool strict)
{
	context().msgLogger(VERBOSITY_HIGH) << "Validating dislocation network." << endl;

	SIMULATION_ASSERT(dislocationLoops().empty() == false);
	int numPeriodicLoops = 0;
	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		NodalVector lineVectorSum(NULL_VECTOR);
		SegmentHandle segment = loop->segments();
		SIMULATION_ASSERT(segment->nextSegment() != segment);	// Loop must contain at least two segments.
		do {
			NodeHandle node1 = segment->node1();
			NodeHandle node2 = segment->node2();
			SIMULATION_ASSERT(node1->nextNode() == node2);
			SIMULATION_ASSERT(node2->prevNode() == node1);
			SIMULATION_ASSERT(segment->prevSegment()->nextSegment() == segment);
			if(segment->isKink()) {
				SIMULATION_ASSERT(node2->pos().X != node1->pos().X || node2->pos().Y != node1->pos().Y);
				bool isKinkDir = false;
				for(auto kd = params().kinkDirections.begin(); kd != params().kinkDirections.end(); ++kd)
					if(segment->lineVector().X == kd->h.X && segment->lineVector().Y == kd->h.Y && positionCompare(segment->lineVector().Z, kd->h.Z))
						isKinkDir = true;
				//SIMULATION_ASSERT(isKinkDir);
			}
			else {
				SIMULATION_ASSERT(node2->pos().X == node1->pos().X && node2->pos().Y == node1->pos().Y);
				SIMULATION_ASSERT(!strict || segment->lineVector().Z != 0);
				//SIMULATION_ASSERT(segment->getHSegmentLength() >= -CAFLOAT_EPSILON);
			}
			NodalVector delta = node2->pos() - node1->pos();
			SIMULATION_ASSERT(delta.X == segment->lineVector().X && delta.Y == segment->lineVector().Y);
			SIMULATION_ASSERT(positionCompare(delta.Z, segment->lineVector().Z));
			SIMULATION_ASSERT(!strict || segment->lineVector() != NULL_VECTOR);
			SIMULATION_ASSERT(segment->node1()->pos().X + segment->lineVector().X == segment->node2()->pos().X);
			SIMULATION_ASSERT(segment->node1()->pos().Y + segment->lineVector().Y == segment->node2()->pos().Y);
			SIMULATION_ASSERT(positionCompare(segment->node1()->pos().Z + segment->lineVector().Z, segment->node2()->pos().Z));
			lineVectorSum += segment->lineVector();
			segment = segment->nextSegment();
		}
		while(segment != loop->segments());

		if(lineVectorSum.X != 0 || lineVectorSum.Y != 0 || positionCompare(lineVectorSum.Z, 0) == false) {
			context().error("Loop with inconsistent line vector sum detected: %i %i %f", lineVectorSum.X, lineVectorSum.Y, lineVectorSum.Z);
		}
		numPeriodicLoops += (int)floor(lineVectorSum.Z / params().lineLength + 0.5);
	}
	SIMULATION_ASSERT(numPeriodicLoops == 1);
}
