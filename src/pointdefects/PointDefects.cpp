#include "PointDefects.h"
#include "../simulation/Simulation.h"
#include "../dislocations/Dislocations.h"
#include <random>

using namespace std;

/******************************************************************************
* Constructor.
*****************************************************************************/
PointDefects::PointDefects(Context& _context, Simulation* sim, SimulationParameters& params) : UsingContext(_context), _sim(sim), _params(params)
{
}

/******************************************************************************
* Initializes the point defect cloud.
*****************************************************************************/
void PointDefects::initialize()
{
	_latticeOffsetX = 0;
	_latticeOffsetY = 0;
	//_latticeOffsetX2 = params().latticeParameter * sqrt(6.0) / 6;
	//_latticeOffsetY2 = params().latticeParameter * sqrt(2.0) / 6;
	//_latticeOffsetX = params().latticeParameter * sqrt(6.0) / 6;
	//_latticeOffsetY = params().latticeParameter * sqrt(2.0) / 6;
}

Point3 PointDefects::getWorldPosition(double x, double y, double z) const {
		// Transform from lattice to world space coordinates.
		Point3 p = params().unitCell * Point3(x,y,z);
		p.X -= _latticeOffsetX;
		p.Y -= _latticeOffsetY;
		return p;
}
/******************************************************************************
* Adjusts the window around the dislocation within which point defects are explicitly modeled.
*****************************************************************************/
void PointDefects::updateWindow()
{
	// Do nothing if point defects are not enabled.
	if(params().pointDefectConcentration == 0 || params().pointDefectWindowRadius == 0) return;

	// Calculate the current position of the screw dislocation in the XY plane.
	DislocationStats stats = simulation().network().computeDislocationPosition();

	// Convert to lattice space.
	Point3 lsp = params().inverseUnitCell * Point3(stats.screwPositionX - _latticeOffsetX, stats.screwPositionY - _latticeOffsetY, 0);
	int newCenterX = (int)floor(lsp.X + 0.5);
	int newCenterY = (int)floor(lsp.Y + 0.5);
	int newRadius = (int)floor(params().pointDefectWindowRadius / 0.816496580);

	static const int sectorToLattice[6][2][2] = {
			{{1,0},{0,1}},
			{{0,-1},{1,1}},
			{{-1,-1},{1,0}},
			{{-1,0},{0,-1}},
			{{0,1},{-1,-1}},
			{{1,1},{-1,0}}
	};

	// Remove rows of defects which are no longer within the window.

	for(auto row = defects().begin(); row != defects().end(); ) {
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
			PointDefect* pd = row->second;
			// period condition
			PointDefect* head = nullptr;
			double new_x;
			if(x>0){
			new_x=newCenterX-(x-1);
			}
			else if(x==0){
			new_x=newCenterX;
			}
			else{
			new_x=newCenterX-(x+1);
			}
			double new_y=y;
			if(y>0){
			new_y=newCenterY-(y-1);
			}
			else if(y==0){
			new_y=newCenterY;
			}
			else{
			new_y=newCenterY-(y+1);
			}
			auto entrynew = defects().find(make_pair(new_x,new_y));

			while(pd != nullptr) {
				PointDefect* nextpd = pd->next;
				//period condition
				if (entrynew != defects().end()){
				PointDefect* pnew = entrynew->second;
				bool isOccupied = false;
				for(PointDefect* pt = pnew; pt != nullptr; pt = pt->next) {
							if(pt->position == pd->position) {
								isOccupied = true;
								break;
							}
						}
				if(!isOccupied){
					PointDefect* defect = _defectPool.construct();
					defect->position = pd->position;
					defect->type = pd->type;
					defect->p_type = pd->p_type;
					defect->x = new_x;
					defect->y = new_y;
					defect->p0 = getWorldPosition(defect->x,defect->y,defect->position);//[modify for moving solute
					defect->p1 = getWorldPosition(defect->x,defect->y,defect->position);//
					defect->next = nullptr;

					PointDefect* pl = pnew;
					if(pnew != nullptr){
					while(pl->next != nullptr){
						pl = pl->next;
					}
					pl->next = defect;
					}
					else{
						entrynew->second = defect;//maybe wrong
					}
				}	
				}
				else{
				PointDefect* defect = _defectPool.construct();
				defect->position = pd->position;
				defect->type = pd->type;
				defect->p_type = pd->p_type;
				defect->x = new_x;
				defect->y = new_y;
				defect->p0 = getWorldPosition(defect->x,defect->y,defect->position);//[modify for moving solute
				defect->p1 = getWorldPosition(defect->x,defect->y,defect->position);//]
				defect->next = head;
				head=defect;
				}

				double msmr = simulation().simulationPdefectMeanSquareMotionReduced();
				msmr += (pd->p1.X - pd->p0.X)*(pd->p1.X - pd->p0.X) + (pd->p1.Y - pd->p0.Y)*(pd->p1.Y - pd->p0.Y) + (pd->p1.Z - pd->p0.Z)*(pd->p1.Z - pd->p0.Z);
				simulation().modifySimulationPdefectMeanSquareMotionReduced(msmr);

				_defectPool.destroy(pd);
				pd = nextpd;
			}
			row = defects().erase(row);
			//period condition
			if (entrynew != defects().end()) defects().insert(make_pair(make_pair(new_x,new_y),head));
		}
	}
	
	// Average number of defects per [111] atomic row.
	double averageDefectsPerRow = params().pointDefectConcentration * params().lineLength;
	std::poisson_distribution<int> defectsPerRowDistr(averageDefectsPerRow);
	std::uniform_int_distribution<int> defectDistr(0, 3 * (params().lineLength - 1));
	std::random_device rd;
    std::mt19937 gen(rd());
	

	// This lambda function randomly positions new point defects on the given atomic row.

	// Create new rows of defects.
	for(int i = 0; i <= newRadius; i++) {
		for(int j = 0; j <= newRadius - i; j++) {
			for(int sector = 0; sector < 6; sector++) {
				double x = sectorToLattice[sector][0][0] * i + sectorToLattice[sector][0][1] * j + newCenterX;
				double y = sectorToLattice[sector][1][0] * i + sectorToLattice[sector][1][1] * j + newCenterY;			
						auto entry = defects().find(make_pair(x,y));
		if(entry != defects().end()) continue;
		// Determine how many defects to distribute on this new atomic row.
		int ndefects = defectsPerRowDistr(gen);
		// Create point defects.
		PointDefect* head = nullptr;
		// Modify MSMR for new defects.
		double msmr = simulation().simulationPdefectMeanSquareMotionReduced();
		double msm = simulation().simulationPdefectMeanSquareMotion();
		msmr += msm * ndefects;
		simulation().modifySimulationPdefectMeanSquareMotionReduced(msmr);
		for(int k = 0; k < ndefects; k++) {
			int a;
			double z;
			int t;//the Position offset number.
			for(;;) {
				a = defectDistr(gen);
				t = a % 3;
				z = a/3.0;
				if(x>=y){
					if ((int)(x-y) % 3 != 0) z += 1.0-(int)(x-y) % 3 /3.0; //the first clockwise.
				}
				else{
					z += (int)(y-x) % 3 /3.0;//the second anti-clockwise.
				}
				// Check if lattice site is already occupied by another defect.
				bool isOccupied = false;
				for(PointDefect* pd = head; pd != nullptr; pd = pd->next) {
					if(pd->position == z) {
						isOccupied = true;
						break;
					}
				}
				if(!isOccupied) break;
			}
			PointDefect* defect = _defectPool.construct();
			defect->position = z;
			defect->type = 1;
			defect->p_type = t;
			double x_l;
			double y_l;
			if (t == 0){
				x_l= 2.0/12.0;
				y_l= -1.0/12.0;
			}
			if (t == 1){
				x_l= -1.0/12.0;
				y_l= -1.0/12.0;
			}
			if (t == 2){
				x_l= -1.0/12.0;
				y_l= 2.0/12.0;
			} 
			defect->x = x+x_l;
			defect->y = y+y_l;

			defect-> p0 = getWorldPosition(defect->x,defect->y,z);//[modify for moving solute
			defect-> p1 = getWorldPosition(defect->x,defect->y,z);//]

			defect->next = head;
			head = defect;
		}
		defects().insert(make_pair(make_pair(x,y),head));
		//context().msgLogger(VERBOSITY_NORMAL) << "x,y: " << x << y << endl;
			}
		}
	}

	// Create new rows of defects for second lattice.
	for(int i = 1; i <= newRadius; i++) {
		for(int j = 0; j <= newRadius - i; j++) {
			for(int sector = 0; sector < 6; sector++) {
				double x = sectorToLattice[sector][0][0] * i + sectorToLattice[sector][0][1] * j + newCenterX +1/3.0;
				double y = sectorToLattice[sector][1][0] * i + sectorToLattice[sector][1][1] * j + newCenterY +1/3.0;			
						auto entry = defects().find(make_pair(x,y));
		if(entry != defects().end()) continue;
		// Determine how many defects to distribute on this new atomic row.
		int ndefects = defectsPerRowDistr(gen);
		// Create point defects.
		PointDefect* head = nullptr;
		// Modify MSMR for new defects.
		double msmr = simulation().simulationPdefectMeanSquareMotionReduced();
		double msm = simulation().simulationPdefectMeanSquareMotion();
		msmr += msm * ndefects;
		simulation().modifySimulationPdefectMeanSquareMotionReduced(msmr);
		for(int k = 0; k < ndefects; k++) {
			int a;
			double z;
			int t;//the Position offset number.
			for(;;) {
				a = defectDistr(gen);
				t = a % 3 + 3;
				z = a/3.0;
				if(x>=y){
					if ((int)(x-y) % 3 != 0) z += 1.0-(int)(x-y) % 3 /3.0;
				}
				else{
					z += (int)(y-x) % 3 /3.0;
				}
				// Check if lattice site is already occupied by another defect.
				bool isOccupied = false;
				for(PointDefect* pd = head; pd != nullptr; pd = pd->next) {
					if(pd->position == z) {
						isOccupied = true;
						break;
					}
				}
				if(!isOccupied) break;
			}
			PointDefect* defect = _defectPool.construct();
			defect->position = z;
			defect->type = 1;
			defect->p_type = t;
			double x_l;
			double y_l;
			if (t == 5){
				x_l= -2.0/12.0;
				y_l= 1.0/12.0;
			}
			if (t == 4){
				x_l= 1.0/12.0;
				y_l= 1.0/12.0;
			}
			if (t == 3){
				x_l= 1.0/12.0;
				y_l= -2.0/12.0;
			} 
			defect->x = x+x_l;
			defect->y = y+y_l;
			defect-> p0 = getWorldPosition(defect->x,defect->y,z);//[modify for moving solute
			defect-> p1 = getWorldPosition(defect->x,defect->y,z);//]
			defect->next = head;
			head = defect;
		}
		defects().insert(make_pair(make_pair(x,y),head));
		//context().msgLogger(VERBOSITY_NORMAL) << "x1,y1: " << x << y << endl;
			}
		}
	}
}
bool PointDefects::migration(PointDefect* p, double x, double y, double z, double x_line, double y_line, double x1, double y1, double z1,double x1_line, double y1_line, int t){
	///
	auto entry = defects().find(make_pair(x1_line,y1_line));
	SIMULATION_ASSERT(p != nullptr);
	if(entry != defects().end()){
		//context().msgLogger(VERBOSITY_NORMAL) <<"good!!!!!!!"<< endl;
		PointDefect* pd = entry->second;
		bool isOccupied = false;
		for(PointDefect* pt = pd; pt != nullptr; pt = pt->next) {
					if(pt->position == z1) {
						isOccupied = true;
						break;
					}
				}
		//context().msgLogger(VERBOSITY_NORMAL) << "isOccupied: " << isOccupied << endl;
		if(!isOccupied){
			PointDefect* defect = _defectPool.construct();
			defect->position = z1;
			defect->x = x1;
			defect->y = y1;
			defect->p_type = t;
			defect->type = 1;
			defect->p0 = p->p0;
			defect->p1 = getWorldPosition(x1,y1,z1);
			
			double msm = simulation().simulationPdefectMeanSquareMotion();
			int msm_number = simulation().number();
			double msm_change = - (p->p1.X - p->p0.X)*(p->p1.X - p->p0.X) - (p->p1.Y - p->p0.Y)*(p->p1.Y - p->p0.Y) - (p->p1.Z - p->p0.Z)*(p->p1.Z - p->p0.Z) + (defect->p1.X - defect->p0.X)*(defect->p1.X - defect->p0.X) + (defect->p1.Y - defect->p0.Y)*(defect->p1.Y - defect->p0.Y) + (defect->p1.Z - defect->p0.Z)*(defect->p1.Z - defect->p0.Z);
			msm = msm + msm_change / msm_number;
			simulation().modifySimulationPdefectMeanSquareMotion(msm);

			defect->next = nullptr;
			PointDefect* pl = pd;
			if(pd != nullptr){
			while(pl->next != nullptr){
				pl = pl->next;
			}
			pl->next = defect;
			}
			else{
				entry->second = defect;//maybe wrong
			}
		}
		else{
			return false;
		}
	}
	else{
		PointDefect* head = nullptr;
		PointDefect* defect = _defectPool.construct();
		defect->position = z1;
		defect->x = x1;
		defect->y = y1;
		defect->p_type = t;
		defect->type = 1;
		defect->p0 = p->p0;
		defect->p1 = getWorldPosition(x1,y1,z1);
		defect->next = head;
		head = defect;

		double msm = simulation().simulationPdefectMeanSquareMotion();
		int msm_number = simulation().number();
		double msm_change = - (p->p1.X - p->p0.X)*(p->p1.X - p->p0.X) - (p->p1.Y - p->p0.Y)*(p->p1.Y - p->p0.Y) - (p->p1.Z - p->p0.Z)*(p->p1.Z - p->p0.Z) + (defect->p1.X - defect->p0.X)*(defect->p1.X - defect->p0.X) + (defect->p1.Y - defect->p0.Y)*(defect->p1.Y - defect->p0.Y) + (defect->p1.Z - defect->p0.Z)*(defect->p1.Z - defect->p0.Z);
		msm = msm + msm_change / msm_number;
		simulation().modifySimulationPdefectMeanSquareMotion(msm);
		
		defects().insert(make_pair(make_pair(x1_line,y1_line),head));	
	}//create the new pointdefect poisition after migration;

	auto entry2 = defects().find(make_pair(x_line,y_line));
	if(entry2 != defects().end()){
		PointDefect* pw = entry2->second;
		SIMULATION_ASSERT(pw != nullptr);
		int count = 0;
		for(PointDefect* ph = pw; ph != nullptr; ph = ph->next) {
			count++;
		}
		SIMULATION_ASSERT(count > 0);
		if(count == 1){
		_defectPool.destroy(pw);
		entry2->second = nullptr;
		}
		else{
			if (pw->position == p->position){
				SIMULATION_ASSERT(pw!=nullptr);
				PointDefect* pwnext = pw->next;
				_defectPool.destroy(pw);
				entry2->second = pwnext;
			}
			else{
				PointDefect* pybefore;
			for(PointDefect* py = pw; py != nullptr; py = py->next) {
				SIMULATION_ASSERT(py->next != nullptr);
				if (py->next->position == p->position){
					pybefore = py;
					PointDefect* pynext = py->next;
					pybefore->next = pynext->next;
					_defectPool.destroy(pynext);
					break;
				}
			}
			}
		}
	}
	return true;

}

double PointDefects::bindPointDefects(const Point3& p0,const Point3& p1,double& deltaZ){
	double Z = deltaZ;
	double Z2 = deltaZ;
	double Z3 = deltaZ;
	for(auto row = defects().begin(); row != defects().end(); ++row) {
		double x = row->first.first;
		double y = row->first.second;
		PointDefect* head = row->second;
		for(PointDefect* pd = head; pd != nullptr; pd = pd->next) {
			double z = pd->position;
			Point3 p2 = getWorldPosition(x,y,z);
			Point3 p3 = params().inverseUnitCell * p2;
			if(deltaZ >0)
			{
				if((fabs(p0.X-p3.X)<=0.001)&&(fabs(p0.Y-p3.Y)<=0.001)&&(p3.Z >= p0.Z && (p0.Z + deltaZ) >= p3.Z))
				{
					Z = p3.Z - p0.Z;
					if(Z < Z3)
						Z3 = Z;
				}
				if((fabs(p0.X-p3.X)<=0.001)&&(fabs(p0.Y-p3.Y)<=0.001)&&(p3.Z >= (p0.Z + params().lineLength) && (p0.Z + params().lineLength + deltaZ) >= p3.Z))
				{
					Z = p3.Z - (p0.Z + params().lineLength);
					if(Z < Z3)
						Z3 = Z;
				}
				if((fabs(p0.X-p3.X)<=0.001)&&(fabs(p0.Y-p3.Y)<=0.001)&&(p3.Z >= (p0.Z - params().lineLength) && (p0.Z - params().lineLength + deltaZ) >= p3.Z))
				{
					Z = p3.Z - (p0.Z - params().lineLength);
					if(Z < Z3)
						Z3 = Z;
				}


				if((fabs(p1.X-p3.X)<=0.001)&&(fabs(p1.Y-p3.Y)<=0.001)&&(p3.Z >= p1.Z && (p1.Z + deltaZ) >= p3.Z))
				{
					Z2 = p3.Z - p1.Z;
					if(Z2 < Z3)
						Z3 = Z2;
				}
				if((fabs(p1.X-p3.X)<=0.001)&&(fabs(p1.Y-p3.Y)<=0.001)&&(p3.Z>= (p1.Z+params().lineLength) && (p1.Z + params().lineLength + deltaZ) >= p3.Z))
				{
					Z2 = p3.Z - (p1.Z+params().lineLength);
					if(Z2 < Z3)
						Z3 = Z2;
				}
				if((fabs(p1.X-p3.X)<=0.001)&&(fabs(p1.Y-p3.Y)<=0.001)&&(p3.Z >= (p1.Z - params().lineLength) && (p1.Z - params().lineLength+ deltaZ) >= p3.Z))
				{
					Z2 = p3.Z - (p1.Z-params().lineLength);
					if(Z2 < Z3)
						Z3 = Z2;
				}

			}
			else
			{

				if((fabs(p0.X-p3.X)<=0.001)&&(fabs(p0.Y-p3.Y)<=0.001)&&(p3.Z <= p0.Z && (p0.Z + deltaZ) <= p3.Z))
				{
					Z = p3.Z - p0.Z;
					if (Z > Z3)
						Z3 = Z;
				}
				if((fabs(p0.X-p3.X)<=0.001)&&(fabs(p0.Y-p3.Y)<=0.001)&&(p3.Z<= (p0.Z + params().lineLength) && (p0.Z + params().lineLength + deltaZ) <= p3.Z))
				{
					Z = p3.Z - (p0.Z + params().lineLength);
					if (Z > Z3)
						Z3 = Z;
				}
				if((fabs(p0.X-p3.X)<=0.001)&&(fabs(p0.Y-p3.Y)<=0.001)&&(p3.Z <= (p0.Z - params().lineLength) && (p0.Z - params().lineLength+ deltaZ) <= p3.Z))
				{
					Z = p3.Z - (p0.Z - params().lineLength);
					if (Z > Z3)
						Z3 = Z;
				}


				if((fabs(p1.X-p3.X)<=0.001)&&(fabs(p1.Y-p3.Y)<=0.001)&&(p3.Z <= p1.Z && (p1.Z + deltaZ) <= p3.Z))
				{
					Z2 = p3.Z - p1.Z;
					if (Z2 > Z3)
						Z3 = Z2;
				}
				if((fabs(p1.X-p3.X)<=0.001)&&(fabs(p1.Y-p3.Y)<=0.001)&&(p3.Z<= (p1.Z + params().lineLength) && (p1.Z + params().lineLength + deltaZ) <= p3.Z))
				{
					Z2 = p3.Z - (p1.Z + params().lineLength);
					if (Z2 > Z3)
						Z3 = Z2;
				}
				if((fabs(p1.X-p3.X)<=0.001)&&(fabs(p1.Y-p3.Y)<=0.001)&&(p3.Z <= (p1.Z - params().lineLength) && (p1.Z - params().lineLength+ deltaZ) <= p3.Z))
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

bool PointDefects::isSoluteOnTheDislocation(const Point3& p1, const Point3& p2)
{		
		bool Is = false;
		Vector3 l = p2 - p1;
		//context().msgLogger(VERBOSITY_NORMAL) << "p1.x: " << p1.X <<"p1.y: " << p1.Y <<"p1.z: " << p1.Z << "p2.z: "<< p2.Z<< endl;
		for(auto row = defects().begin(); row != defects().end(); ++row) {
		double x = row->first.first;
		double y = row->first.second;
		PointDefect* head = row->second;
		for(PointDefect* pd = head; pd != nullptr; pd = pd->next) {
			double z = pd->position;
			Point3 p3 = getWorldPosition(x,y,z);
			Point3 p4 = params().inverseUnitCell * p3;
			//context().msgLogger(VERBOSITY_NORMAL) << "p4.x: " << p4.X <<"p4.y: " << p4.Y <<"p4.z: " << p4.Z << endl;
			//Point3 p4 = params().inverseUnitCell * pd->p1;
			/*Vector3 l1 = p3 - p1;
			if(fabs(sqrt((l.X * l.X +l.Y * l.Y + l.Z * l.Z)*(l1.X * l1.X +l1.Y * l1.Y + l1.Z * l1.Z)) - (l.X * l1.X + l.Y * l1.Y + l.Z * l1.Z)) <= 0.01)
			{
				if((l.X * l.X + l.Y * l.Y + l.Z * l.Z) >= (l1.X * l1.X +l1.Y * l1.Y + l1.Z * l1.Z))
				{
					Is = true;
				}
			}*/
			if((fabs(p1.X-p4.X)<=0.001)&&(fabs(p1.Y-p4.Y)<=0.001)&&((fabs(p1.Z-p4.Z)<=0.001)||(fabs(p1.Z-p4.Z+params().lineLength)<=0.001)||(fabs(p1.Z-p4.Z-params().lineLength)<=0.001)))
				Is = true;
			if((fabs(p2.X-p4.X)<=0.001)&&(fabs(p2.Y-p4.Y)<=0.001)&&((fabs(p2.Z-p4.Z)<=0.001)||(fabs(p2.Z-p4.Z+params().lineLength)<=0.001)||(fabs(p2.Z-p4.Z-params().lineLength)<=0.001)))
				Is = true;
		}
		}
		//context().msgLogger(VERBOSITY_NORMAL) << "isSoluteOnTheDislocation: " << Is << endl;
		return Is;
}

bool PointDefects::isSoluteOnTheScrewDislocation(const Point3& p1, const Point3& p2)
{
		bool Is = false;
		SIMULATION_ASSERT(p1.X == p2.X && p1.Y == p2.Y);

		for(auto row = defects().begin(); row != defects().end(); ++row) {
		double x = row->first.first;
		double y = row->first.second;
		PointDefect* head = row->second;
			for(PointDefect* pd = head; pd != nullptr; pd = pd->next) {
			double z = pd->position;
			Point3 p3 = getWorldPosition(x,y,z);
			Point3 p4 = params().inverseUnitCell * p3;
			//Point3 p4 = params().inverseUnitCell * pd->p1;
			if((fabs(p1.X - p4.X)<=0.001)&&(fabs(p1.Y - p4.Y)<=0.001)&&(((p1.Z <= p4.Z)&&(p2.Z >= p4.Z))||((p1.Z <= (p4.Z + params().lineLength))&&(p2.Z >= (p4.Z + params().lineLength)))||((p1.Z <= (p4.Z - params().lineLength))&&(p2.Z >= (p4.Z - params().lineLength)))))
				Is = true;
			}
		}
		return Is;
}
