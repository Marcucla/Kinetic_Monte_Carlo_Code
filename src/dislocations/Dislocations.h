#ifndef __DISLOCATIONS_H
#define __DISLOCATIONS_H

#include "../StandardIncludes.h"
#include "../Settings.h"
#include "../context/Context.h"

class Simulation;

struct SimulationParameters;


/// Stores the position of a node. X and Y components are always integer lattice positions, while Z can vary continuously.
typedef Point_3<int,int,double> NodalPosition;

/// Stores a vector between two nodes. X and Y components are always integer lattice positions, while Z can vary continuously.
typedef Vector_3<int,int,double> NodalVector;

/// Epsilon value that defines when two positions are considered equal.
#define NODAL_POS_EPSILON		1e-6

/// Epsilon value that defines when two nodal velocities are considered equal.
#define NODAL_VELOCITY_EPSILON	1e-6

/*
 * A nodal point of a discretized dislocation line.
 *
 * This class also represents the adjacent segment,
 * connecting the node with the following node in the discretized line.
 *
 * Nodes/segments form a circular doubly-linked list (closed dislocation loop).
 *
 * There are two types of segments: H-segments (pure screw) and V-segments (kinks).
 */
class DislocationNode
{
public:

	/// Screw segment constructor.
	DislocationNode(const NodalPosition& pos, double screwSegmentLength) :
		_pos(pos), _lineVector(0,0,screwSegmentLength), _kinkVelocity(0), _nodeVelocity(0), _kinkDirection(-1) {}

	/// Kink segment constructor.
	DislocationNode(const NodalPosition& pos, const NodalVector& kinkVector, int kinkDirection) :
		_pos(pos), _lineVector(kinkVector), _kinkVelocity(0), _nodeVelocity(0), _kinkDirection(kinkDirection) {
		SIMULATION_ASSERT(kinkDirection != -1);
	}

	/**************************** Node specific functions *****************************/

	/// Returns the next node in the dislocation loop.
	DislocationNode* nextNode() const { return _next; }

	/// Returns the previous node in the dislocation loop.
	DislocationNode* prevNode() const { return _prev; }

	/// Returns the node's position on the lattice grid.
	const NodalPosition& pos() const { return _pos; }

	/// Changes the node's position on the lattice grid.
	void setPos(const NodalPosition& p) { _pos = p; }

	/// Returns the ingoing segment adjacent to this node.
	DislocationNode* inSegment() const { return _prev; }

	/// Returns the outgoing segment adjacent to this node.
	DislocationNode* outSegment() const { return const_cast<DislocationNode*>(this); }

	/// Returns velocity of this node in units of Burgers vectors per second.
	double nodeVelocity() const { return _nodeVelocity; }

	/// Sets the velocity of this node in units of Burgers vectors per second.
	void setNodeVelocity(double v) { _nodeVelocity = v; }

	/**************************** Segment specific functions *****************************/

	/// Returns the lattice space vector of this segment.
	const NodalVector& lineVector() const { return _lineVector; }

	/// Changes the lattice space vector of this screw segment.
	void setScrewLineVector(double length) { _lineVector = NodalVector(0,0,length); }

	/// Changes the lattice space vector of this kink segment.
	void setKinkLineVector(int kinkDirection, const NodalVector& v) { _lineVector = v; _kinkDirection = kinkDirection; }

	/// Returns the next segment in the dislocation loop.
	DislocationNode* nextSegment() const { return _next; }

	/// Returns the previous segment in the dislocation loop.
	DislocationNode* prevSegment() const { return _prev; }

	/// True if this is a kink segment (V-segment).
	/// False if this is a pure screw segment (H-segment).
	bool isKink() const { return _lineVector.X != 0 || _lineVector.Y != 0; }

	/// True if this is a pure screw segment (H-segment).
	/// False if this is a kink segment (V-segment).
	bool isScrew() const { return _lineVector.X == 0 && _lineVector.Y == 0; }

	/// Returns a pointer to the first node of this segment.
	DislocationNode* node1() const { return const_cast<DislocationNode*>(this); }

	/// Returns a pointer to the second node of this segment.
	DislocationNode* node2() const { return _next; }

	/// Returns the length of this H-segment (can be positive or negative) [in units of b].
	double getHSegmentLength() const {
		SIMULATION_ASSERT(isKink() == false);
		SIMULATION_ASSERT(isScrew() == true);
		SIMULATION_ASSERT(node1()->pos().X == node2()->pos().X);
		SIMULATION_ASSERT(node1()->pos().Y == node2()->pos().Y);
		return _lineVector.Z;
	}

	/// Returns the direction index if this is a kink segment.
	int kinkDirection() const { SIMULATION_ASSERT(isKink()); return _kinkDirection; }

	/// Returns velocity of this V-segment in units of Burgers vectors per second.
	double kinkVelocity() const { return _kinkVelocity; }

	/// Sets the velocity of this V-segment in units of Burgers vectors per second.
	void setKinkVelocity(double v) { _kinkVelocity = v; }

	/// Returns true if this segment is part of an infinite straight screw segment without kinks.
	bool isStraightScrewLine() const {
		SIMULATION_ASSERT(isScrew());
		SIMULATION_ASSERT(nextSegment()->nextSegment() != this ||
				nextSegment()->isScrew());
		return nextSegment()->nextSegment() == this;
	}

	/**************************** General functions *****************************/

	/// Sets the next node/segment in the doubly-linked list.
	void setNext(DislocationNode* next) { _next = next; }

	/// Sets the previous node/segment in the doubly-linked list.
	void setPrev(DislocationNode* prev) { _prev = prev; }

private:

	/// The node's position in lattice-space.
	NodalPosition _pos;

	/// The lattice-space vector of this segment.
	NodalVector _lineVector;

	/// Pointer to the next node / segment in the dislocation loop.
	DislocationNode* _next;

	/// Pointer to the previous node / segment in the dislocation loop.
	DislocationNode* _prev;

	/// The direction index if this is a kink segment.
	int _kinkDirection;

	/// The velocity of the V-segment in units of Burgers vectors per second.
	double _kinkVelocity;

	/// The velocity of the node in units of Burgers vectors per second.
	double _nodeVelocity;

	/// The random diffusion direction of a kink segment (-1=negative direction, +1=positive direction).
	int _kinkDiffusionDir;

	/// The random diffusion direction of a node (-1=negative direction, +1=positive direction).
	int _nodeDiffusionDir;

	friend class DislocationNetwork;
};

/// Define the linear segment type. This is the same as the node type because they are stored
/// in the same structure.
typedef DislocationNode DislocationSegment;

/// Define handle types for nodes.
typedef DislocationNode* NodeHandle;

/// Define handle types for nodes.
typedef DislocationSegment* SegmentHandle;

struct DislocationBindingEventList
{
	double rate;
	DislocationNode* p;
	int t;
};

/*
 * A dislocation line.
 */
class DislocationLoop
{
public:

	/// Constructor
	DislocationLoop(DislocationNode* headNode) : _nodesAndSegments(headNode) {
		SIMULATION_ASSERT(headNode != NULL);
		SIMULATION_ASSERT(headNode->nextNode() != headNode);	// Loop must contain at least two nodes/segments.
	}

	/// Returns the linked list of nodes of this dislocation line.
	NodeHandle nodes() const { return _nodesAndSegments; }

	/// Returns the linked list of segments of this dislocation line.
	SegmentHandle segments() const { return _nodesAndSegments; }

	/// Sets the head node of the linked list.
	void setHead(NodeHandle node) {
		SIMULATION_ASSERT(node->nextNode() != node);	// Loop must contain at least two nodes/segments.
		_nodesAndSegments = node;
	}

private:

	/// Pointer to the head of the linked-list of nodes/segments.
	DislocationNode* _nodesAndSegments;
};

/**
 * Stores global and statistical information, characterizing the current dislocation configuration.
 */
struct DislocationStats
{
	double screwPositionX;	// Center of mass position of the screw dislocation along X-axis (units: Angstrom).
	double screwPositionY;	// Center of mass position of the screw dislocation along Y-axis (units: Angstrom).
	double lineWidthX;		// Width of dislocation line in X direction (units: Angstrom).
	double lineWidthY;		// Width of dislocation line in Y direction (units: Angstrom).
	int numLoops;			// Number of loops.
	int numSegments;		// Total number of dislocation segments.
	int numScrewSegments;	// Number of pure screw dislocation segments.
	int numKinkSegments;	// Number of kink dislocation segments.
};

/**
 * Stores the current dislocation configuration.
 */
class DislocationNetwork : public UsingContext
{
public:

	/// Constructor.
	DislocationNetwork(Context& context, Simulation* sim, SimulationParameters& params);

	/// Returns a reference to the list of dislocation loops.
	std::list<DislocationLoop>& dislocationLoops() { return _dislocations; }

	/// Returns a const reference to the list of dislocation loops.
	const std::list<DislocationLoop>& dislocationLoops() const { return _dislocations; }

	/// Returns a reference to the simulation's global parameter set.
	const SimulationParameters& params() const { return _params; }

	Simulation& simulation() const { return *_sim; }

	/// Checks whether the data structure is in good order. This is for debugging only.
	void validate(bool strict);

	/// Prints a text representation of the current dislocation structure to the console.
	void printNetwork();

	/// Dumps the current dislocation configuration to an output file.
	void dumpToVTK(std::ostream& stream, double simulationTime);

	/********************************** Query functions *********************************/

	/// Computes the world space vector connecting the two nodes of the given segment.
	Vector3 getSegmentVectorWorld(SegmentHandle segment) const;

	/// Computes the world space center of mass of a segment.
	Point3 getSegmentCenter(SegmentHandle segment) const;

	/// Determines the earliest time until two nodes will collide based on their current velocities.
	double freeMotionTime();

	/// Determines the earliest time until segement will bind with solute atoms based on their current velocities.
	double minTimeForBind(double t_mig);

	/// Determines the time the fastest kink segment needs to travel the given maximum distance (in Burgers vector units).
	double minTimeForDistance(double distance);

	/// Compares two lattice coordinates for equality. Takes into account periodicity of the simulation domain.
	bool positionCompare(const NodalPosition& p1, const NodalPosition& p2) const;

	/// Compares two lattice coordinates for equality. Takes into account periodicity of the simulation domain.
	bool positionCompare(double z1, double z2) const;

	/// Determines whether the given 4-junction should be flipped based on a curvature criterion.
	bool shouldFlipJunction(NodeHandle node1, NodeHandle node2);

	/// Calculates the current position of the screw dislocation in the XY plane. This is used to
	/// measure the dislocation velocity.
	DislocationStats computeDislocationPosition();

	/********************************** Elastic energy functions *********************************/

	/// Calculates the local stress tensor at the given point.
	SymmetricTensor2 calculateLocalStress(const Point3& p, bool includeExternalStress = true);

	/// Calculates the local resolved shear stress at the given point and along the given kink direction.
	double calculateResolvedShearStress(const Point3& p, int kinkDirection, bool includeExternalStress = true);

	/// Computes the elastic interaction energy of a (virtual) dislocation loop with the existing dislocations
	/// and its self energy. The rectangular loop is given by four points.
	double elasticLoopEnergy(const Point3& p1, const Point3& p2, const Point3& p3, const Point3& p4);

	/// Computes the elastic self-energy of a dislocation segment.
	double elasticSelfEnergy(const Vector3& lineVector);

	/// Calculates the elastic interaction energy of a (virtual) straight segment and the existing dislocations.
	double elasticSegmentNetworkEnergy(const Point3& x1, const Point3& x2);

	/// Computes the elastic interaction energy between two straight dislocation segments.
	/// Does not include effects of periodic images.
	double elasticSegmentSegmentEnergy(const Point3& x1, const Point3& x2, const Point3& x3, const Point3& x4, const Vector3& b);

	/// Computes the stress at a point due to the given dislocation segment.
	SymmetricTensor2 segmentStress(const Point3& p, const Point3& p1, const Point3& p2, const Vector3& b);

	/********************************** Modification functions *********************************/

	/// Create a single straight screw dislocation.
	void initialize();

	/// Inserts a new loop into the dislocation configuration.
	DislocationLoop& createLoop(NodeHandle headNode);

	/// Makes sure that a node is not connected to itself with a segment.
	/// Inserts an additional node into the segment if necessary.
	void eliminateSelfSegment(NodeHandle headNode);

	/// Deletes the given segment from the given loop.
	/// Returns the segment that follows the deleted segment in the loop.
	SegmentHandle deleteSegment(SegmentHandle segment, DislocationLoop& loop);

	/// Inserts a kink pair into an H-segment.
	bool insertKinkPair(SegmentHandle segment, double kinkPairPos, int kinkDirectionIndex, double kinkPairWidth);

	/// Removes redundant segments and nodes as well as degenerate loops.
	void cleanupNetwork();

	/// Finds attracting cross-kinks that form a locked configuration.
	bool detectLockedCrossKinks();

	/// Determines whether a cross-kink wants to stay locked or rather dissociate.
	bool areCrossKinksLocked(SegmentHandle screwSegment);

	/// Collapses a short screw segment such that the adjacent kinks form a locked cross-kink.
	SegmentHandle lockCrossKink(SegmentHandle screwSegment, DislocationLoop& loop);

	/// Calculates the velocity from the force acting on a kink segment.
	double calculateKinkVelocity(SegmentHandle kinkSegment);

	/// Assigns velocities to the kink segments.
	void calculateKinkVelocities(boost::mt19937& rng);

	/// Advances the simulation by the given amount and lets kinks move.
	void propagateKinks(double delta_t);

	/// Moves a kink segment along the screw direction by the given amount (in units of b).
	/// The caller of this method has to make sure that the displacement will not lead to
	/// an invalid dislocation configuration.
	void displaceKink(SegmentHandle kink, double deltaZ);

	bool isSoluteOnTheDislocation(const Point3& p, bool& iskinkbool);

	//bool executeBindingEvent(SegmentHandle& devent);
	bool executeBindingEvent(DislocationBindingEventList& devent);

	double calculateLocalStressKink(SegmentHandle kinkSegment);

private:

	/// The list of dislocation loops.
	std::list<DislocationLoop> _dislocations;

	/// The memory pool for allocating new nodes.
	boost::object_pool<DislocationNode> _nodePool;

	/// Reference to the simulation's global parameters.
	SimulationParameters& _params;

	Simulation* _sim;
};

/**
 * An iterator class that is used to iterate over all dislocation segments in
 * all dislocation loops of the network.
 */
class SegmentIterator
{
public:

	/// Constructor.
	SegmentIterator(const DislocationNetwork& network) : _loop(network.dislocationLoops().begin()), _loopEnd(network.dislocationLoops().end()) {
		if(_loop != _loopEnd) {
			_segment = _loop->segments();
			SIMULATION_ASSERT(_segment->nextSegment() != _segment);
		}
		else _segment = NULL;
	}

	/// Returns the current segment handle pointed to by the iterator.
	SegmentHandle operator*() const { SIMULATION_ASSERT(_segment != NULL); return _segment; }

	/// Returns the current segment handle pointed to by the iterator.
	SegmentHandle operator->() const { SIMULATION_ASSERT(_segment != NULL); return _segment; }

	/// Tests for equality of two iterators.
	bool operator==(const SegmentIterator& other) const { return _segment == other._segment; }

	/// Tests for inequality of two iterators.
	bool operator!=(const SegmentIterator& other) const { return _segment != other._segment; }

	/// Prefix increment operator.
	SegmentIterator& operator++() {
		SIMULATION_ASSERT(_loop != _loopEnd && _segment != NULL);
		_segment = _segment->nextSegment();
		if(_segment == _loop->segments()) {
			++_loop;
			if(_loop != _loopEnd) {
				_segment = _loop->segments();
				SIMULATION_ASSERT(_segment->nextSegment() != _segment);
			}
			else _segment = NULL;
		}
		return *this;
	}

	/// Returns whether this iterator has reached the end of the segment list.
	bool atEnd() const { return _segment == NULL; }

private:
	/// The current loop.
	std::list<DislocationLoop>::const_iterator _loop;
	/// end-of-loop-list pointer.
	std::list<DislocationLoop>::const_iterator _loopEnd;
	/// The current segment.
	SegmentHandle _segment;
};

#endif // __DISLOCATIONS_H
