#include "PointRes.h"
#include "../simulation/Simulation.h"

using namespace std;

/*******************************************************************************
* Dumps the point defects to an output file.
********************************************************************************/
void PointRes::dumpToVTK(ostream& stream, double simulationTime)
{
	// Count the current number of defects.
	size_t reCount = 0;
	for(const auto& row : res()) {
		for(PointRe* pd = row.second; pd != nullptr; pd = pd->next)
			reCount++;
	}

	stream << "# vtk DataFile Version 3.0" << endl;
	stream << "# Point res (t=" << simulationTime << ")" << endl;
	stream << "ASCII" << endl;
	stream << "DATASET UNSTRUCTURED_GRID" << endl;
	stream << "POINTS " << reCount << " double" << endl;
	for(const auto& row : res()) {
		for(PointRe* pd = row.second; pd != nullptr; pd = pd->next) {
			Point3 wp = getWorldPosition(row.first.first, row.first.second, pd->position);
			stream << wp.X << " " << wp.Y << " " << (wp.Z * params().outputScalingFactor) << "\n";
		}
	}

	stream << endl << "CELLS " << reCount << " " << (reCount*2) << endl;
	for(size_t i = 0; i < reCount; i++)
		stream << "1 " << i << "\n";

	stream << endl << "CELL_TYPES " << reCount << "\n";
	for(size_t i = 0; i < reCount; i++)
		stream << "1" << "\n";	// Vertex

	stream << endl << "CELL_DATA " << reCount << "\n";

	stream << endl << "SCALARS re_type int 1" << "\n";
	stream << endl << "LOOKUP_TABLE default" << "\n";
	for(const auto& row : res())
		for(PointRe* pd = row.second; pd != nullptr; pd = pd->next)
			stream << pd->type << "\n";
}
