#include "Dislocations.h"
#include "../simulation/Simulation.h"

using namespace std;

/*******************************************************************************
* Produces a text representation of the dislocation network for diagnostic purposes.
********************************************************************************/
void DislocationNetwork::printNetwork()
{
	context().msgLogger() << "Number of dislocation loops: " << dislocationLoops().size() << endl;
	context().msgLogger() << "Box size: " << params().lineLength << " [b] = " << params().pbcLength << " [A]" << endl;
	size_t loopIndex = 1;
	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop, ++loopIndex) {
		context().msgLogger() << "  Loop " << loopIndex << ":" << endl;
		size_t nodeIndex = 0;
		SegmentHandle segment = loop->segments();
		do {
			context().msgLogger() << "    Node " << nodeIndex << " at " << segment->node1()->pos() << (params().unitCell * segment->node1()->pos()) << (segment->isKink() ? " V-segment" : " H-segment") << "  delta=" << segment->lineVector() << endl;
			segment = segment->nextSegment();
			nodeIndex++;
		}
		while(segment != loop->segments());
	}
}

/*******************************************************************************
* Dumps the current dislocation configuration to an output file.
********************************************************************************/
void DislocationNetwork::dumpToVTK(ostream& stream, double simulationTime)
{
	vector<Point3> outputPoints;
	vector< pair< pair<int,int>, int> > outputSegments;
	vector< pair<int, int> > nodalPoints;

	double box_size = params().pbcLength;
	int loopIndex = 0;
	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop, ++loopIndex) {
		SegmentHandle segment = loop->segments();
		do {
			Point3 p1 = params().unitCell * segment->node1()->pos();
			while(p1.Z < 0) p1.Z += box_size;
			while(p1.Z > box_size) p1.Z -= box_size;

			if(outputPoints.empty() || outputPoints.back() != p1)
				outputPoints.push_back(p1);

			nodalPoints.push_back(make_pair(outputPoints.size()-1, loopIndex));

			Vector3 delta = params().unitCell * segment->lineVector();
			Point3 p2 = p1 + delta;
			if(p2.Z > box_size) {
				Point3 t = p1 + delta * ((box_size - p1.Z)/delta.Z);
				outputSegments.push_back(make_pair(make_pair(outputPoints.size()-1, outputPoints.size()), loopIndex));
				outputPoints.push_back(t);
				outputSegments.push_back(make_pair(make_pair(outputPoints.size(), outputPoints.size()+1), loopIndex));
				outputPoints.push_back(t + Vector3(0,0,-box_size));
				outputPoints.push_back(p2 + Vector3(0,0,-box_size));
			}
			else if(p2.Z < 0) {
				Point3 t = p1 + delta * (-p1.Z/delta.Z);
				outputSegments.push_back(make_pair(make_pair(outputPoints.size()-1, outputPoints.size()), loopIndex));
				outputPoints.push_back(t);
				outputSegments.push_back(make_pair(make_pair(outputPoints.size(), outputPoints.size()+1), loopIndex));
				outputPoints.push_back(t + Vector3(0,0,box_size));
				outputPoints.push_back(p2 + Vector3(0,0,box_size));
			}
			else {
				outputSegments.push_back(make_pair(make_pair(outputPoints.size()-1, outputPoints.size()), loopIndex));
				outputPoints.push_back(p2);
			}
			segment = segment->nextSegment();
		}
		while(segment != loop->segments());
	}

	stream << "# vtk DataFile Version 3.0" << endl;
	stream << "# Dislocation (t=" << simulationTime << ")" << endl;
	stream << "ASCII" << endl;
	stream << "DATASET UNSTRUCTURED_GRID" << endl;
	stream << "POINTS " << outputPoints.size() << " double" << endl;
	for(auto p = outputPoints.begin(); p != outputPoints.end(); ++p)
		stream << p->X << " " << p->Y << " " << (p->Z * params().outputScalingFactor) << endl;

	stream << endl << "CELLS " << (nodalPoints.size() + outputSegments.size()) << " " << (outputSegments.size() * 3 + nodalPoints.size() * 2) << endl;
	for(auto seg = outputSegments.begin(); seg != outputSegments.end(); ++seg)
		stream << "2 " << seg->first.first << " " << seg->first.second << endl;
	for(auto p = nodalPoints.begin(); p != nodalPoints.end(); ++p)
		stream << "1 " << p->first << endl;

	stream << endl << "CELL_TYPES " << (outputSegments.size() + nodalPoints.size()) << endl;
	for(size_t i = 0; i < outputSegments.size(); i++)
		stream << "3" << endl;	// Line
	for(size_t i = 0; i < nodalPoints.size(); i++)
		stream << "1" << endl;	// Vertex

	stream << endl << "CELL_DATA " << (outputSegments.size() + nodalPoints.size()) << endl;

	stream << endl << "SCALARS loop_index int 1" << endl;
	stream << endl << "LOOKUP_TABLE default" << endl;
	for(auto seg = outputSegments.begin(); seg != outputSegments.end(); ++seg)
		stream << seg->second << endl;
	for(auto p = nodalPoints.begin(); p != nodalPoints.end(); ++p)
		stream << p->second << endl;
}

/*******************************************************************************
* Calculates the current position of the screw dislocation in the XY plane.
* This is used to measure the dislocation velocity.
********************************************************************************/
DislocationStats DislocationNetwork::computeDislocationPosition()
{
	DislocationStats result;
	result.screwPositionX = 0;
	result.screwPositionY = 0;
	result.lineWidthX = 0;
	result.lineWidthY = 0;
	result.numLoops = dislocationLoops().size();
	result.numSegments = 0;
	result.numScrewSegments = 0;
	result.numKinkSegments = 0;

	// Find the loop which spans the simulation box and is periodic.

	// Note that there might be more than one infinite loop if dipoles occurred during the simulation.
	// In such a case, pick the line that is farthest from the origin.

	for(auto loop = dislocationLoops().begin(); loop != dislocationLoops().end(); ++loop) {
		NodalVector lineVectorSum(NULL_VECTOR);
		NodalPosition bboxMin(loop->nodes()->pos());
		NodalPosition bboxMax(loop->nodes()->pos());
		Point3 center(ORIGIN);
		double numHSegments = 0;
		SegmentHandle segment = loop->segments();
		do {
			lineVectorSum += segment->lineVector();
			result.numSegments++;
			if(segment->isScrew()) {
				numHSegments += fabs(segment->getHSegmentLength());
				result.numScrewSegments++;
				center.X += fabs(segment->getHSegmentLength()) * segment->node1()->pos().X;
				center.Y += fabs(segment->getHSegmentLength()) * segment->node1()->pos().Y;
				if(segment->node1()->pos().X > bboxMax.X) bboxMax.X = segment->node1()->pos().X;
				if(segment->node1()->pos().X < bboxMin.X) bboxMin.X = segment->node1()->pos().X;
				if(segment->node1()->pos().Y > bboxMax.Y) bboxMax.Y = segment->node1()->pos().Y;
				if(segment->node1()->pos().Y < bboxMin.Y) bboxMin.Y = segment->node1()->pos().Y;
			}
			else result.numKinkSegments++;
			segment = segment->nextSegment();
		}
		while(segment != loop->segments());
		SIMULATION_ASSERT(lineVectorSum.X == 0 && lineVectorSum.Y == 0);

		if(positionCompare(lineVectorSum.Z, 0)) {
			Point3 centerWorld = params().unitCell * Point3(center.X / numHSegments, center.Y / numHSegments, 0);
			if(square(centerWorld.X) + square(centerWorld.Y) >= square(result.screwPositionX) + square(result.screwPositionY)) {
				result.screwPositionX = centerWorld.X;
				result.screwPositionY = centerWorld.Y;
				Vector3 sizeWorld = params().unitCell * Vector3((bboxMax.X - bboxMin.X), (bboxMax.Y - bboxMin.Y), 0);
				result.lineWidthX = fabs(sizeWorld.X);
				result.lineWidthY = fabs(sizeWorld.Y);
			}
		}
	}
	return result;
}
