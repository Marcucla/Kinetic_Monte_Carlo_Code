#include "Dislocations.h"
#include "../simulation/Simulation.h"
#include "../pointdefects/PointDefects.h"

using namespace std;

/******************************************************************************
* Constructor.
*****************************************************************************/
DislocationNetwork::DislocationNetwork(Context& _context,Simulation* sim, SimulationParameters& params) : UsingContext(_context), _sim(sim), _params(params)
{
}

/******************************************************************************
* Create a single straight screw dislocation.
*****************************************************************************/
void DislocationNetwork::initialize()
{
	// Create the initial straight dislocation line.
	DislocationNode* node1 = _nodePool.construct(ORIGIN, params().lineLength);
	node1->setNext(node1);
	node1->setPrev(node1);
	createLoop(node1);
}

/*******************************************************************************
* Inserts a new loop into the dislocation configuration.
********************************************************************************/
DislocationLoop& DislocationNetwork::createLoop(NodeHandle headNode)
{
	// Make sure loop consists of at least two segments and two nodes.
	eliminateSelfSegment(headNode);

	DislocationLoop line(headNode);
	_dislocations.push_back(line);
	return _dislocations.back();
}

/*******************************************************************************
* Makes sure that a node is not connected to itself with a segment.
* Inserts an additional node into the segment if necessary.
********************************************************************************/
void DislocationNetwork::eliminateSelfSegment(NodeHandle headNode)
{
	if(headNode->nextNode() == headNode) {
		SegmentHandle seg = headNode->outSegment();
		SIMULATION_ASSERT(seg->lineVector().X == 0 && seg->lineVector().Y == 0);
		SIMULATION_ASSERT(fabs(seg->lineVector().Z - params().lineLength) < CAFLOAT_EPSILON);
		// Insert a second node if necessary.
		seg->setScrewLineVector(params().lineLength / 2);
		DislocationNode* node2 = _nodePool.construct(headNode->pos() + seg->lineVector(), params().lineLength - seg->getHSegmentLength());
		context().msgLogger(VERBOSITY_HIGH) << "Created second node at " << node2->pos() << " linevec=" << node2->lineVector() << endl;
		SIMULATION_ASSERT(positionCompare(headNode->pos(), node2->pos()) == false);
		headNode->setNext(node2);
		headNode->setPrev(node2);
		node2->setNext(headNode);
		node2->setPrev(headNode);
	}
}

/*******************************************************************************
* Computes the world space vector connecting the two nodes of the given segment.
********************************************************************************/
Vector3 DislocationNetwork::getSegmentVectorWorld(SegmentHandle segment)  const
{
	// Transform from lattice to world space coordinates.
	return params().unitCell * segment->lineVector();
}

/*******************************************************************************
* Computes the world space center of mass of a segment.
********************************************************************************/
Point3 DislocationNetwork::getSegmentCenter(SegmentHandle segment) const
{
	Point3 p1 = params().unitCell * segment->node1()->pos();
	return p1 + 0.5 * (params().unitCell * segment->lineVector());
}

/*******************************************************************************
* Compares two lattice coordinates for equality.
* Takes into account periodicity of the simulation domain.
********************************************************************************/
bool DislocationNetwork::positionCompare(const NodalPosition& p1, const NodalPosition& p2) const
{
	// X and Y coordinates must be exactly equal.
	// The Z coordinates must be equal modulo the periodicity length.
	return (p1.X == p2.X) && (p1.Y == p2.Y) && positionCompare(p1.Z, p2.Z);
}

/*******************************************************************************
* Compares two lattice coordinates for equality.
* Takes into account periodicity of the simulation domain.
********************************************************************************/
bool DislocationNetwork::positionCompare(double z1, double z2) const
{
	// The Z coordinates must be equal modulo the periodicity length.
	double m = fabs(fmod(z2 - z1, params().lineLength));

	return (fabs(m) <= NODAL_POS_EPSILON || fabs(m - params().lineLength) <= NODAL_POS_EPSILON);
}
