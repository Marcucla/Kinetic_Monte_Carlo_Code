#include <fstream>
#include <cmath>

using namespace std;

const double MU = 160e9;
const double NU = 0.28;
const double a0 = 3.18;
const double b = a0*sqrt(3.0)/2.0;
const double b_spatial[3] = {0,0,b};
const double h = a0*sqrt(6.0)/3.0;
const double kinkWidth = 20.0;

void segmentStress(const double p[3], const double p1[3], const double p2[3], const double b[3], double stress[6], double a);

void calculateStressAlongSharpKink(const char* filename, double a)
{
	ofstream stream(filename);
	double A[3] = {h,0,-2000};
	double B[3] = {h,0,0};
	double C[3] = {0,0,0};
	double D[3] = {0,0,2000};
	for(int i = -1000; i<=1000; i++) {
		double z = (double)i / 50.0;

		double probePoint[3] = { 0, 0, b*z };
		if(z < 0) probePoint[0] = h;
		double totalStress[6] = {0};
		segmentStress(probePoint, A, B, b_spatial, totalStress, a);
		segmentStress(probePoint, B, C, b_spatial, totalStress, a);
		segmentStress(probePoint, C, D, b_spatial, totalStress, a);
		stream << z << " " << totalStress[4] << endl;
	}
}

void calculateStressAlongBeveledKink(const char* filename, double a, bool onLine)
{
	ofstream stream(filename);
	double A[3] = {h,0,-2000};
	double B[3] = {h,0,-kinkWidth*0.5*b};
	double C[3] = {0,0,+kinkWidth*0.5*b};
	double D[3] = {0,0,2000};
	for(int i = -1000; i<=1000; i++) {
		double z = (double)i / 50.0;

		double probePoint[3] = { 0, 0, b*z };
		if(z < 0) probePoint[0] = h;

		if(onLine && fabs(z) < kinkWidth/2) {
			double t = 0.5 - z / kinkWidth;
			probePoint[0] = h * t;
		}

		double totalStress[6] = {0};
		segmentStress(probePoint, A, B, b_spatial, totalStress, a);
		segmentStress(probePoint, B, C, b_spatial, totalStress, a);
		segmentStress(probePoint, C, D, b_spatial, totalStress, a);
		stream << z << " " << totalStress[4] << endl;
	}
}

int main(int argc, char* argv[])
{
	calculateStressAlongSharpKink("sharp_kink_0.1b.dat", 0.1 * b);
	calculateStressAlongSharpKink("sharp_kink_0.5b.dat", 0.5 * b);
	calculateStressAlongSharpKink("sharp_kink_1b.dat", 1 * b);
	calculateStressAlongBeveledKink("wide_kink_0.5b.dat", 0.5 * b, false);
	calculateStressAlongBeveledKink("wide_kink_on_0.5b.dat", 0.5 * b, true);

	return 0;
}

/******************************************************************************
* Computes the stress at a point due to the given dislocation segment.
*****************************************************************************/
void segmentStress(const double p[3], const double p1[3], const double p2[3], const double b[3], double stress[6], double a)
{
	double oneoverLp, common;
	double vec1x, vec1y, vec1z;
	double tpx, tpy, tpz;
	double Rx, Ry, Rz, Rdt;
	double ndx, ndy, ndz;
	double d2, s1, s2, a2, a2_d2, a2d2inv;
	double Ra, Rainv, Ra3inv, sRa3inv;
	double s_03a, s_13a, s_05a, s_15a, s_25a;
	double s_03b, s_13b, s_05b, s_15b, s_25b;
	double s_03, s_13, s_05, s_15, s_25;
	double m4p, m8p, m4pn, mn4pn, a2m8p;
	double txbx, txby, txbz;
	double dxbx, dxby, dxbz;
	double dxbdt, dmdxx, dmdyy, dmdzz, dmdxy, dmdyz, dmdxz;
	double tmtxx, tmtyy, tmtzz, tmtxy, tmtyz, tmtxz;
	double tmdxx, tmdyy, tmdzz, tmdxy, tmdyz, tmdxz;
	double tmtxbxx, tmtxbyy, tmtxbzz, tmtxbxy, tmtxbyz, tmtxbxz;
	double dmtxbxx, dmtxbyy, dmtxbzz, dmtxbxy, dmtxbyz, dmtxbxz;
	double tmdxbxx, tmdxbyy, tmdxbzz, tmdxbxy, tmdxbyz, tmdxbxz;
	double I_03xx, I_03yy, I_03zz, I_03xy, I_03yz, I_03xz;
	double I_13xx, I_13yy, I_13zz, I_13xy, I_13yz, I_13xz;
	double I_05xx, I_05yy, I_05zz, I_05xy, I_05yz, I_05xz;
	double I_15xx, I_15yy, I_15zz, I_15xy, I_15yz, I_15xz;
	double I_25xx, I_25yy, I_25zz, I_25xy, I_25yz, I_25xz;

	vec1x = p2[0] - p1[0];
	vec1y = p2[1] - p1[1];
	vec1z = p2[2] - p1[2];

	oneoverLp = 1 / sqrt(vec1x*vec1x + vec1y*vec1y + vec1z*vec1z);

	tpx = vec1x * oneoverLp;
	tpy = vec1y * oneoverLp;
	tpz = vec1z * oneoverLp;

	Rx = p[0] - p1[0];
	Ry = p[1] - p1[1];
	Rz = p[2] - p1[2];

	Rdt = Rx*tpx + Ry*tpy + Rz*tpz;

	ndx = Rx - Rdt*tpx;
	ndy = Ry - Rdt*tpy;
	ndz = Rz - Rdt*tpz;

	d2 = ndx*ndx + ndy*ndy + ndz*ndz;

	s1 = -Rdt;
	s2 = -((p[0]-p2[0])*tpx + (p[1]-p2[1])*tpy + (p[2]-p2[2])*tpz);
	a2 = a * a;
	a2_d2 = a2 + d2;
	a2d2inv = 1 / a2_d2;

	Ra = sqrt(a2_d2 + s1*s1);
	Rainv = 1 / Ra;
	Ra3inv = Rainv * Rainv * Rainv;
	sRa3inv = s1 * Ra3inv;

	s_03a = s1 * Rainv * a2d2inv;
	s_13a = -Rainv;
	s_05a = (2*s_03a + sRa3inv) * a2d2inv;
	s_15a = -Ra3inv;
	s_25a = s_03a - sRa3inv;

	Ra = sqrt(a2_d2 + s2*s2);
	Rainv = 1 / Ra;
	Ra3inv = Rainv * Rainv * Rainv;
	sRa3inv = s2 * Ra3inv;

	s_03b = s2 * Rainv * a2d2inv;
	s_13b = -Rainv;
	s_05b = (2*s_03b + sRa3inv) * a2d2inv;
	s_15b = -Ra3inv;
	s_25b = s_03b - sRa3inv;

	s_03 = s_03b - s_03a;
	s_13 = s_13b - s_13a;
	s_05 = s_05b - s_05a;
	s_15 = s_15b - s_15a;
	s_25 = s_25b - s_25a;

	m4p = 0.25 * MU / M_PI;
	m8p = 0.5 * m4p;
	m4pn = m4p / (1 - NU);
	mn4pn = m4pn * NU;
	a2m8p = a2 * m8p;


	txbx = tpy*b[2] - tpz*b[1];
	txby = tpz*b[0] - tpx*b[2];
	txbz = tpx*b[1] - tpy*b[0];

	dxbx = ndy*b[2] - ndz*b[1];
	dxby = ndz*b[0] - ndx*b[2];
	dxbz = ndx*b[1] - ndy*b[0];

	dxbdt = dxbx*tpx + dxby*tpy + dxbz*tpz;

	dmdxx = ndx * ndx;
	dmdyy = ndy * ndy;
	dmdzz = ndz * ndz;
	dmdxy = ndx * ndy;
	dmdyz = ndy * ndz;
	dmdxz = ndx * ndz;

	tmtxx = tpx * tpx;
	tmtyy = tpy * tpy;
	tmtzz = tpz * tpz;
	tmtxy = tpx * tpy;
	tmtyz = tpy * tpz;
	tmtxz = tpx * tpz;

	tmdxx = 2 * tpx * ndx;
	tmdyy = 2 * tpy * ndy;
	tmdzz = 2 * tpz * ndz;
	tmdxy = tpx*ndy + tpy*ndx;
	tmdyz = tpy*ndz + tpz*ndy;
	tmdxz = tpx*ndz + tpz*ndx;


	tmtxbxx = 2 * tpx * txbx;
	tmtxbyy = 2 * tpy * txby;
	tmtxbzz = 2 * tpz * txbz;
	tmtxbxy = tpx*txby + tpy*txbx;
	tmtxbyz = tpy*txbz + tpz*txby;
	tmtxbxz = tpx*txbz + tpz*txbx;

	dmtxbxx = 2 * ndx * txbx;
	dmtxbyy = 2 * ndy * txby;
	dmtxbzz = 2 * ndz * txbz;
	dmtxbxy = ndx*txby + ndy*txbx;
	dmtxbyz = ndy*txbz + ndz*txby;
	dmtxbxz = ndx*txbz + ndz*txbx;


	tmdxbxx = 2 * tpx * dxbx;
	tmdxbyy = 2 * tpy * dxby;
	tmdxbzz = 2 * tpz * dxbz;
	tmdxbxy = tpx*dxby + tpy*dxbx;
	tmdxbyz = tpy*dxbz + tpz*dxby;
	tmdxbxz = tpx*dxbz + tpz*dxbx;

	common = m4pn * dxbdt;

	I_03xx = common + m4pn*dmtxbxx - m4p*tmdxbxx;
	I_03yy = common + m4pn*dmtxbyy - m4p*tmdxbyy;
	I_03zz = common + m4pn*dmtxbzz - m4p*tmdxbzz;
	I_03xy = m4pn*dmtxbxy - m4p*tmdxbxy;
	I_03yz = m4pn*dmtxbyz - m4p*tmdxbyz;
	I_03xz = m4pn*dmtxbxz - m4p*tmdxbxz;

	I_13xx = -mn4pn * tmtxbxx;
	I_13yy = -mn4pn * tmtxbyy;
	I_13zz = -mn4pn * tmtxbzz;
	I_13xy = -mn4pn * tmtxbxy;
	I_13yz = -mn4pn * tmtxbyz;
	I_13xz = -mn4pn * tmtxbxz;

	I_05xx = common*(a2+dmdxx) - a2m8p*tmdxbxx;
	I_05yy = common*(a2+dmdyy) - a2m8p*tmdxbyy;
	I_05zz = common*(a2+dmdzz) - a2m8p*tmdxbzz;
	I_05xy = common*dmdxy - a2m8p*tmdxbxy;
	I_05yz = common*dmdyz - a2m8p*tmdxbyz;
	I_05xz = common*dmdxz - a2m8p*tmdxbxz;

	I_15xx = a2m8p*tmtxbxx - common*tmdxx;
	I_15yy = a2m8p*tmtxbyy - common*tmdyy;
	I_15zz = a2m8p*tmtxbzz - common*tmdzz;
	I_15xy = a2m8p*tmtxbxy - common*tmdxy;
	I_15yz = a2m8p*tmtxbyz - common*tmdyz;
	I_15xz = a2m8p*tmtxbxz - common*tmdxz;

	I_25xx = common * tmtxx;
	I_25yy = common * tmtyy;
	I_25zz = common * tmtzz;
	I_25xy = common * tmtxy;
	I_25yz = common * tmtyz;
	I_25xz = common * tmtxz;

	stress[0] += I_03xx*s_03 + I_13xx*s_13 + I_05xx*s_05 +
				I_15xx*s_15 + I_25xx*s_25;

	stress[1] += I_03yy*s_03 + I_13yy*s_13 + I_05yy*s_05 +
				I_15yy*s_15 + I_25yy*s_25;

	stress[2] += I_03zz*s_03 + I_13zz*s_13 + I_05zz*s_05 +
				I_15zz*s_15 + I_25zz*s_25;

	stress[3] += I_03xy*s_03 + I_13xy*s_13 + I_05xy*s_05 +
				I_15xy*s_15 + I_25xy*s_25;

	stress[4] += I_03yz*s_03 + I_13yz*s_13 + I_05yz*s_05 +
				I_15yz*s_15 + I_25yz*s_25;

	stress[5] += I_03xz*s_03 + I_13xz*s_13 + I_05xz*s_05 +
				I_15xz*s_15 + I_25xz*s_25;
}
