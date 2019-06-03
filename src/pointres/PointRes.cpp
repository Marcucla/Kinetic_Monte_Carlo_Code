#include "PointRes.h"
#include "../simulation/Simulation.h"
#include "../dislocations/Dislocations.h"
#include <random>

using namespace std;

/******************************************************************************
* Constructor.
*****************************************************************************/
PointRes::PointRes(Context& _context, Simulation* sim, SimulationParameters& params) : UsingContext(_context), _sim(sim), _params(params)
{
}

/******************************************************************************
* Initializes the point defect cloud.
*****************************************************************************/
void PointRes::initialize()
{
	//_latticeOffsetX = 0;
	//_latticeOffsetY = 0;
	_latticeOffsetX = params().latticeParameter * sqrt(6.0) / 6;
	_latticeOffsetY = params().latticeParameter * sqrt(2.0) / 6;
}

Point3 PointRes::getWorldPosition(int x, int y, int z) const {
		// Transform from lattice to world space coordinates.
		Point3 p = params().unitCell * Point3(x,y,z);
		p.X -= _latticeOffsetX;
		p.Y -= _latticeOffsetY;
		return p;
}
/******************************************************************************
* Adjusts the window around the dislocation within which point defects are explicitly modeled.
*****************************************************************************/
void PointRes::updateWindow()
{
	// Do nothing if point defects are not enabled.
	if(params().pointReConcentration == 0 || params().pointReWindowRadius == 0) return;

	// Calculate the current position of the screw dislocation in the XY plane.
	DislocationStats stats = simulation().network().computeDislocationPosition();

	// Convert to lattice space.
	Point3 lsp = params().inverseUnitCell * Point3(stats.screwPositionX - _latticeOffsetX, stats.screwPositionY - _latticeOffsetY, 0);
	int newCenterX = (int)floor(lsp.X + 0.5);
	int newCenterY = (int)floor(lsp.Y + 0.5);
	int newRadius = (int)floor(params().pointReWindowRadius / 0.816496580);

	static const int sectorToLattice[6][2][2] = {
			{{1,0},{0,1}},
			{{0,-1},{1,1}},
			{{-1,-1},{1,0}},
			{{-1,0},{0,-1}},
			{{0,1},{-1,-1}},
			{{1,1},{-1,0}}
	};

	// Remove rows of defects which are no longer within the window.
	for(auto row = res().begin(); row != res().end(); ) {
		int x = row->first.first - newCenterX;
		int y = row->first.second - newCenterY;
		int dist;
		if(x >= 0) {
			if(y >= 0) dist = x+y;
			else if(x >= -y) dist = x;
			else dist = -y;
		}
		else {
			if(y <= 0) dist = -(x+y);
			else if(y <= -x) dist = -x;
			else dist = y;
		}
		if(dist <= newRadius)
			++row;
		else {
			PointRe* pd = row->second;
			while(pd != nullptr) {
				PointRe* nextpd = pd->next;
				_rePool.destroy(pd);
				pd = nextpd;
			}
			row = res().erase(row);
		}
	}

	// Average number of defects per [111] atomic row.
	double averageResPerRow = params().pointReConcentration * params().lineLength;
	std::poisson_distribution<int> resPerRowDistr(averageResPerRow);
	std::uniform_int_distribution<int> reDistr(0, params().lineLength - 1);
	std::random_device rd;
    std::mt19937 gen(rd());
	

	// This lambda function randomly positions new point defects on the given atomic row.

	// Create new rows of defects.
	for(int i = 1; i <= newRadius; i++) {
		for(int j = 0; j <= newRadius - i; j++) {
			for(int sector = 0; sector < 6; sector++) {
				int x = sectorToLattice[sector][0][0] * i + sectorToLattice[sector][0][1] * j + newCenterX;
				int y = sectorToLattice[sector][1][0] * i + sectorToLattice[sector][1][1] * j + newCenterY;
						auto entry = res().find(make_pair(x,y));
		if(entry != res().end()) continue;
		// Determine how many defects to distribute on this new atomic row.
		int nres = resPerRowDistr(gen);
		// Create point defects.
		PointRe* head = nullptr;
		for(int k = 0; k < nres; k++) {
			int z;
			for(;;) {
				z = reDistr(gen);
				// Check if lattice site is already occupied by another defect.
				bool isOccupied = false;
				for(PointRe* pd = head; pd != nullptr; pd = pd->next) {
					if(pd->position == z) {
						isOccupied = true;
						break;
					}
				}
				if(!isOccupied) break;
			}
			PointRe* re = _rePool.construct();
			re->position = z;
			re->type = 1;
			re->next = head;
			head = re;
		}
		res().insert(make_pair(make_pair(x,y),head));
			}
		}
	}
}
bool PointRes::migration(PointRe* p, int x, int y, int z, int x1, int y1, int z1){
	auto entry = res().find(make_pair(x1,y1));
	SIMULATION_ASSERT(p != nullptr);
	if(entry != res().end()){
		PointRe* pd = entry->second;
		bool isOccupied = false;
		for(PointRe* pt = pd; pt != nullptr; pt = pt->next) {
					if(pt->position == z1) {
						isOccupied = true;
						break;
					}
				}
		if(!isOccupied){
			PointRe* re = _rePool.construct();
			re->position = z1;
			re->type = 1;
			re->next = nullptr;
			PointRe* pl = pd;
			if(pd != nullptr){
			while(pl->next != nullptr){
				pl = pl->next;
			}
			pl->next = re;
			}
			else
			{
				pd = re;
			}
		}
	}//create the new pointdefect poisition after migration;

	auto entry2 = res().find(make_pair(x,y));
	if(entry2 != res().end()){
		PointRe* pw = entry2->second;
		SIMULATION_ASSERT(pw != nullptr);
		int count = 0;
		for(PointRe* ph = pw; ph != nullptr; ph = ph->next) {
			count++;
		}
		SIMULATION_ASSERT(count > 0);
		if(count == 1){
		_rePool.destroy(pw);
		entry2->second = nullptr;
		}
		else{
			if (pw->position == p->position){
				SIMULATION_ASSERT(pw!=nullptr);
				PointRe* pwnext = pw->next;
				_rePool.destroy(pw);
				entry2->second = pwnext;
			}
			else{
				PointRe* pybefore;
			for(PointRe* py = pw; py != nullptr; py = py->next) {
				SIMULATION_ASSERT(py->next != nullptr);
				if (py->next->position == p->position){
					pybefore = py;
					PointRe* pynext = py->next;
					pybefore->next = pynext->next;
					_rePool.destroy(pynext);
					break;
				}
			}
			}
		}
	}
	return true;

}

double PointRes::bindPointRes(const Point3& p0,const Point3& p1,double& deltaZ){
	double Z = deltaZ;
	double Z2 = deltaZ;
	double Z3 = deltaZ;
	for(auto row = res().begin(); row != res().end(); ++row) {
		int x = row->first.first;
		int y = row->first.second;
		PointRe* head = row->second;
		for(PointRe* pd = head; pd != nullptr; pd = pd->next) {
			int z = pd->position;
			Point3 p2 = getWorldPosition(x,y,z);
			Point3 p3 = params().inverseUnitCell * p2;
			if(deltaZ >0)
			{
				if((fabs(p0.X-p3.X)<=0.0001)&&(fabs(p0.Y-p3.Y)<=0.0001)&&(p3.Z >= p0.Z && (p0.Z + deltaZ) >= p3.Z))
				{
					Z = p3.Z - p0.Z;
					if(Z < Z3)
						Z3 = Z;
				}
				if((fabs(p0.X-p3.X)<=0.0001)&&(fabs(p0.Y-p3.Y)<=0.0001)&&(p3.Z >= (p0.Z + params().lineLength) && (p0.Z + params().lineLength + deltaZ) >= p3.Z))
				{
					Z = p3.Z - (p0.Z + params().lineLength);
					if(Z < Z3)
						Z3 = Z;
				}
				if((fabs(p0.X-p3.X)<=0.0001)&&(fabs(p0.Y-p3.Y)<=0.0001)&&(p3.Z >= (p0.Z - params().lineLength) && (p0.Z - params().lineLength + deltaZ) >= p3.Z))
				{
					Z = p3.Z - (p0.Z - params().lineLength);
					if(Z < Z3)
						Z3 = Z;
				}


				if((fabs(p1.X-p3.X)<=0.0001)&&(fabs(p1.Y-p3.Y)<=0.0001)&&(p3.Z >= p1.Z && (p1.Z + deltaZ) >= p3.Z))
				{
					Z2 = p3.Z - p1.Z;
					if(Z2 < Z3)
						Z3 = Z2;
				}
				if((fabs(p1.X-p3.X)<=0.0001)&&(fabs(p1.Y-p3.Y)<=0.0001)&&(p3.Z>= (p1.Z+params().lineLength) && (p1.Z + params().lineLength + deltaZ) >= p3.Z))
				{
					Z2 = p3.Z - (p1.Z+params().lineLength);
					if(Z2 < Z3)
						Z3 = Z2;
				}
				if((fabs(p1.X-p3.X)<=0.0001)&&(fabs(p1.Y-p3.Y)<=0.0001)&&(p3.Z >= (p1.Z - params().lineLength) && (p1.Z - params().lineLength+ deltaZ) >= p3.Z))
				{
					Z2 = p3.Z - (p1.Z-params().lineLength);
					if(Z2 < Z3)
						Z3 = Z2;
				}

			}
			else
			{

				if((fabs(p0.X-p3.X)<=0.0001)&&(fabs(p0.Y-p3.Y)<=0.0001)&&(p3.Z <= p0.Z && (p0.Z + deltaZ) <= p3.Z))
				{
					Z = p3.Z - p0.Z;
					if (Z > Z3)
						Z3 = Z;
				}
				if((fabs(p0.X-p3.X)<=0.0001)&&(fabs(p0.Y-p3.Y)<=0.0001)&&(p3.Z<= (p0.Z + params().lineLength) && (p0.Z + params().lineLength + deltaZ) <= p3.Z))
				{
					Z = p3.Z - (p0.Z + params().lineLength);
					if (Z > Z3)
						Z3 = Z;
				}
				if((fabs(p0.X-p3.X)<=0.0001)&&(fabs(p0.Y-p3.Y)<=0.0001)&&(p3.Z <= (p0.Z - params().lineLength) && (p0.Z - params().lineLength+ deltaZ) <= p3.Z))
				{
					Z = p3.Z - (p0.Z - params().lineLength);
					if (Z > Z3)
						Z3 = Z;
				}


				if((fabs(p1.X-p3.X)<=0.0001)&&(fabs(p1.Y-p3.Y)<=0.0001)&&(p3.Z <= p1.Z && (p1.Z + deltaZ) <= p3.Z))
				{
					Z2 = p3.Z - p1.Z;
					if (Z2 > Z3)
						Z3 = Z2;
				}
				if((fabs(p1.X-p3.X)<=0.0001)&&(fabs(p1.Y-p3.Y)<=0.0001)&&(p3.Z<= (p1.Z + params().lineLength) && (p1.Z + params().lineLength + deltaZ) <= p3.Z))
				{
					Z2 = p3.Z - (p1.Z + params().lineLength);
					if (Z2 > Z3)
						Z3 = Z2;
				}
				if((fabs(p1.X-p3.X)<=0.0001)&&(fabs(p1.Y-p3.Y)<=0.0001)&&(p3.Z <= (p1.Z - params().lineLength) && (p1.Z - params().lineLength+ deltaZ) <= p3.Z))
				{
					Z2 = p3.Z - (p1.Z - params().lineLength);
					if (Z2 > Z3)
						Z3 = Z2;
				}
			}
			//Vector3 l = p1 - p0;
			//SIMULATION_ASSERT(l.Z == 0);
			//Vector3 l1 = p3 - p0;
			/*if(sqrt((l.X * l.X +l.Y * l.Y) * (l1.X * l1.X +l1.Y * l1.Y )) == (l.X * l1.X + l.Y * l1.Y ))
			{
				if((l.X * l.X +l.Y * l.Y ) >= (l1.X * l1.X +l1.Y * l1.Y))
				{
					if(deltaZ >= 0)
					{
						if((l.X = 0)&&((p3.Z >= p0.Z)&&((p3.Z - p0.Z)<= deltaZ)))
							return (p3.Z - p0.Z);
						else
						{
							double templ = sqrt(l1.X * l1.X + l1.Y * l1.Y)/sqrt(l.X * l.X + l.Y * l.Y) * (p1.Z - p0.Z);
							if ((p3.Z >= p0.Z + templ) && (p3.Z <= (p0.Z + templ + deltaZ)))
							{
							if((p3.Z - p0.Z - templ)< Z)
							Z = p3.Z - p0.Z - templ;
							}
						}
					}
					else{
						if((l.X = 0)&&((p3.Z <= p0.Z)&&((p3.Z - p0.Z)>= deltaZ)))
							return (p3.Z - p0.Z);
						else{
							double templ = sqrt(l1.X * l1.X + l1.Y * l1.Y)/sqrt(l.X * l.X + l.Y * l.Y) * (p1.Z - p0.Z);
							if ((p3.Z <= p0.Z + templ) && (p3.Z >= (p0.Z + deltaZ + templ)))
							{
							if((p3.Z - p0.Z - templ)> Z)
							Z = p3.Z - p0.Z - templ;
							}
						}
					}
				}
			}*/

		}
	}
	return Z3;
}

bool PointRes::isSoluteOnTheDislocation(const Point3& p1, const Point3& p2)
{		
		bool Is = false;
		Vector3 l = p2 - p1;

		for(auto row = res().begin(); row != res().end(); ++row) {
		int x = row->first.first;
		int y = row->first.second;
		PointRe* head = row->second;
		for(PointRe* pd = head; pd != nullptr; pd = pd->next) {
			int z = pd->position;
			Point3 p3 = getWorldPosition(x,y,z);
			Point3 p4 = params().inverseUnitCell * p3;
			/*Vector3 l1 = p3 - p1;
			if(fabs(sqrt((l.X * l.X +l.Y * l.Y + l.Z * l.Z)*(l1.X * l1.X +l1.Y * l1.Y + l1.Z * l1.Z)) - (l.X * l1.X + l.Y * l1.Y + l.Z * l1.Z)) <= 0.01)
			{
				if((l.X * l.X + l.Y * l.Y + l.Z * l.Z) >= (l1.X * l1.X +l1.Y * l1.Y + l1.Z * l1.Z))
				{
					Is = true;
				}
			}*/
			if((fabs(p1.X-p4.X)<=0.000001)&&(fabs(p1.Y-p4.Y)<=0.000001)&&((fabs(p1.Z-p4.Z)<=0.000001)||(fabs(p1.Z-p4.Z+params().lineLength)<=0.000001)||(fabs(p1.Z-p4.Z-params().lineLength)<=0.000001)))
				Is = true;
			if((fabs(p2.X-p4.X)<=0.000001)&&(fabs(p2.Y-p4.Y)<=0.000001)&&((fabs(p2.Z-p4.Z)<=0.000001)||(fabs(p2.Z-p4.Z+params().lineLength)<=0.000001)||(fabs(p2.Z-p4.Z-params().lineLength)<=0.000001)))
				Is = true;
		}
		}
		return Is;
}

bool PointRes::isSoluteOnTheScrewDislocation(const Point3& p1, const Point3& p2)
{
		bool Is = false;
		SIMULATION_ASSERT(p1.X == p2.X && p1.Y == p2.Y);

		for(auto row = res().begin(); row != res().end(); ++row) {
		int x = row->first.first;
		int y = row->first.second;
		PointRe* head = row->second;
			for(PointRe* pd = head; pd != nullptr; pd = pd->next) {
			int z = pd->position;
			Point3 p3 = getWorldPosition(x,y,z);
			Point3 p4 = params().inverseUnitCell * p3;
			if((fabs(p1.X - p4.X)<=0.000001)&&(fabs(p1.Y - p4.Y)<=0.000001)&&(((p1.Z <= p4.Z)&&(p2.Z >= p4.Z))||((p1.Z <= (p4.Z + params().lineLength))&&(p2.Z >= (p4.Z + params().lineLength)))||((p1.Z <= (p4.Z - params().lineLength))&&(p2.Z >= (p4.Z - params().lineLength)))))
				Is = true;
			}
		}
		return Is;
}
