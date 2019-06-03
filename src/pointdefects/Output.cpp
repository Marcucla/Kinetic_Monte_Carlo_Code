#include "PointDefects.h"
#include "../simulation/Simulation.h"

using namespace std;

/*******************************************************************************
* Dumps the point defects to an output file.
********************************************************************************/
void PointDefects::dumpToVTK(ostream& stream, double simulationTime)
{
	// Count the current number of defects.
	size_t defectCount = 0;
	for(const auto& row : defects()) {
		for(PointDefect* pd = row.second; pd != nullptr; pd = pd->next)
			defectCount++;
	}

	stream << "# vtk DataFile Version 3.0" << endl;
	stream << "# Point defects (t=" << simulationTime << ")" << endl;
	stream << "ASCII" << endl;
	stream << "DATASET UNSTRUCTURED_GRID" << endl;
	stream << "POINTS " << defectCount << " double" << endl;
	for(const auto& row : defects()) {
		for(PointDefect* pd = row.second; pd != nullptr; pd = pd->next) {
			//Point3 wp = getWorldPosition(row.first.first, row.first.second, pd->position);
			Point3 wp = getWorldPosition(pd->x, pd->y, pd->position);
			stream << wp.X << " " << wp.Y << " " << (wp.Z * params().outputScalingFactor) << "\n";
		}
	}

	stream << endl << "CELLS " << defectCount << " " << (defectCount*2) << endl;
	for(size_t i = 0; i < defectCount; i++)
		stream << "1 " << i << "\n";

	stream << endl << "CELL_TYPES " << defectCount << "\n";
	for(size_t i = 0; i < defectCount; i++)
		stream << "1" << "\n";	// Vertex

	stream << endl << "CELL_DATA " << defectCount << "\n";

	stream << endl << "SCALARS defect_type int 1" << "\n";
	stream << endl << "LOOKUP_TABLE default" << "\n";
	for(const auto& row : defects())
		for(PointDefect* pd = row.second; pd != nullptr; pd = pd->next)
			stream << pd->type << "\n";
}
