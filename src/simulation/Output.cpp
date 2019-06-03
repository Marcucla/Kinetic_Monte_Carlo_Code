#include "Simulation.h"

using namespace std;

/******************************************************************************
* Prints a line to the log file.
*****************************************************************************/
void Simulation::printSimulationInfo()
{
	if(params().printEveryNTimesteps <= 0) return;
	if((simulationStep() % params().printEveryNTimesteps) != 0) return;

	DislocationStats dislocInfo = network().computeDislocationPosition();

	context().msgLogger(VERBOSITY_NORMAL) << simulationStep() << " " << simulationTime() << " " <<
			dislocInfo.screwPositionX << " " << dislocInfo.screwPositionY << " " <<
			dislocInfo.lineWidthX << " " << dislocInfo.lineWidthY << " " <<
			dislocInfo.numLoops << " " << dislocInfo.numSegments << " " << _simulationDiffusionRate << " " << _simulationNucleationRate << " " << _simulationDetachmentRate << " "<< _simulationPdefectMeanSquareMotion <<endl;

	_dislocPositionFile << simulationStep() << " " << simulationTime() << " " <<
			dislocInfo.screwPositionX << " " << dislocInfo.screwPositionY << " " <<
			dislocInfo.lineWidthX << " " << dislocInfo.lineWidthY << " " <<
			dislocInfo.numLoops << " " << dislocInfo.numSegments << " " << _simulationDiffusionRate << " " << _simulationNucleationRate << " " << _simulationDetachmentRate << " "<< _simulationPdefectMeanSquareMotion << endl << flush;
}

/******************************************************************************
* Dumps the current dislocation structure to a VTK file.
*****************************************************************************/
void Simulation::outputDislocations()
{
	if(params().outputEveryNTimesteps <= 0) return;
	if((simulationStep() % params().outputEveryNTimesteps) != 0) return;
	if(params().vtkOutputFilename.empty()) return;

	// Generate a filename for this snapshot.
	stringstream filename;
	size_t dot = params().vtkOutputFilename.rfind('.');
	if(dot != string::npos)
		filename << params().vtkOutputFilename.substr(0, dot+1) << simulationStep() << params().vtkOutputFilename.substr(dot);
	else
		filename << params().vtkOutputFilename << simulationStep();

	// Open output file for writing.
	ofstream stream;
	stream.open(filename.str().c_str());
	if(!stream)
		context().error("Failed to open output file for writing: %s", filename.str().c_str());

	// Export dislocation structure to file.
	network().dumpToVTK(stream, _simulationTime);
}

void Simulation::outputPointDefects()
{
	if(params().outputEveryNTimesteps <= 0) return;
	if((simulationStep() % params().outputEveryNTimesteps) != 0) return;
	if(params().pointDefectOutputFilename.empty()) return;

	// Generate a filename for this snapshot.
	stringstream filename;
	size_t dot = params().pointDefectOutputFilename.rfind('.');
	if(dot != string::npos)
		filename << params().pointDefectOutputFilename.substr(0, dot+1) << simulationStep() << params().pointDefectOutputFilename.substr(dot);
	else
		filename << params().pointDefectOutputFilename << simulationStep();

	// Open output file for writing.
	ofstream stream;
	stream.open(filename.str().c_str());
	if(!stream)
		context().error("Failed to open output file for writing: %s", filename.str().c_str());

	// Write point defects to file.
	pointDefects().dumpToVTK(stream, _simulationTime);
}

