#include "Dislocations.h"
#include "../simulation/Simulation.h"

using namespace std;

/******************************************************************************
* Calculates the local stress tensor at the given point.
*****************************************************************************/
SymmetricTensor2 DislocationNetwork::calculateLocalStress(const Point3& p, bool includeExternalStress)
{
	// Verify simulation parameters.
	SIMULATION_ASSERT(params().coreWidth > 0);
	SIMULATION_ASSERT(params().shearModulus > 0);
	SIMULATION_ASSERT(params().poissonRatio > 0);

	// Burgers vector in spatial coordinate system.
	Vector3 b = params().b_spatial;

	SymmetricTensor2 stress(0);

	// Take into account external stress contribution, then add internal contributions.
	if(includeExternalStress)
		stress = params().externalStress;

	// Take into account local stress contributions.
	if(params().enableLocalStress) {
		// Loop over all dislocation segments.
		for(SegmentIterator segment(*this); !segment.atEnd(); ++segment) {
			Point3 p1 = params().unitCell * segment->node1()->pos();
			Point3 p2 = p1 + (params().unitCell * segment->lineVector());

			// Compute stress contributions of several periodic images to either side of the point.
			int d = (int)floor((p.Z - p1.Z) / params().pbcLength);
			int nimages = 5;
			for(int i = d - nimages; i <= d + nimages; i++) {
				Vector3 shift = Vector3(0, 0, params().pbcLength * i);
				stress += segmentStress(p - shift, p1, p2, b);
			}
		}
	}

	return stress;
}

/******************************************************************************
* Calculates the local resolved shear stress at the given point and
* along the given kink direction.
*****************************************************************************/
double DislocationNetwork::calculateResolvedShearStress(const Point3& p, int kinkDirection, bool includeExternalStress)
{
	SIMULATION_ASSERT(kinkDirection >= 0 && kinkDirection < params().kinkDirections.size());

	// Compute local stress tensor.
	SymmetricTensor2 stress = calculateLocalStress(p, includeExternalStress);

	// Kink plane normal vector.
	const Vector3& n = params().kinkDirections[kinkDirection].normal;

	// Compute resolved shear stress (double contraction of stress and Schmid tensor).
	// Only XZ and YZ components matter.
	SIMULATION_ASSERT(n.Z == 0);
	SIMULATION_CHECK_VALUE(stress(0,2));
	SIMULATION_CHECK_VALUE(stress(1,2));
	return stress(0,2) * n.X + stress(1,2) * n.Y;
}

/******************************************************************************
* Computes the elastic interaction energy of a (virtual) dislocation loop
* with the existing dislocations and its self energy. The rectangular loop
* is given by four corner points.
*****************************************************************************/
double DislocationNetwork::elasticLoopEnergy(const Point3& p1, const Point3& p2, const Point3& p3, const Point3& p4)
{
	double loopEnergy = 0;

	// Self-energy of the four segments in the loop.
	loopEnergy += elasticSelfEnergy(p2 - p1);
	loopEnergy += elasticSelfEnergy(p3 - p2);
	loopEnergy += elasticSelfEnergy(p4 - p3);
	loopEnergy += elasticSelfEnergy(p1 - p4);

	// Elastic interaction among the four segments within the loop.
	Vector3 b = params().unitCell * params().b;
	loopEnergy += elasticSegmentSegmentEnergy(p1, p2, p2, p3, b);
	loopEnergy += elasticSegmentSegmentEnergy(p1, p2, p4, p1, b);
	loopEnergy += elasticSegmentSegmentEnergy(p1, p2, p3, p4, b);
	loopEnergy += elasticSegmentSegmentEnergy(p2, p3, p3, p4, b);
	loopEnergy += elasticSegmentSegmentEnergy(p2, p3, p4, p1, b);
	loopEnergy += elasticSegmentSegmentEnergy(p3, p4, p4, p1, b);

	// Elastic interaction of the four segments with the existing dislocation segments.
	loopEnergy += elasticSegmentNetworkEnergy(p1, p2);
	loopEnergy += elasticSegmentNetworkEnergy(p2, p3);
	loopEnergy += elasticSegmentNetworkEnergy(p3, p4);
	loopEnergy += elasticSegmentNetworkEnergy(p4, p1);

	SIMULATION_CHECK_VALUE(loopEnergy);
	return loopEnergy;
}

/******************************************************************************
* Computes the elastic self-energy of a dislocation segment.
*****************************************************************************/
double DislocationNetwork::elasticSelfEnergy(const Vector3& lineVector)
{
	double mu = params().shearModulus;
	double nu = params().poissonRatio;
	double a = params().coreWidth * params().blength;

	// Calculate self-energy of straight dislocation segment using non-singular expression.
	Vector3 b = params().unitCell * params().b;
	double L = Length(lineVector);
	SIMULATION_ASSERT(L > 0.0);
	Vector3 t = lineVector / L;
	context().msgLogger() << "------------------------------" << endl;
	double L_alpha = sqrt(L * L + a * a);
	double b_dot_t = DotProduct(t, b);
#if 0
	// Formula taken from paper:
	double p = (square(params().blength) - nu * square(b_dot_t)) * L * log((L_alpha + L) / a) - (3.0 - nu) * square(b_dot_t) * (L_alpha - a) / 2.0;
	double selfEnergy = mu * p / (4.0 * M_PI * (1.0 - nu));
#else
	// Formula taken from Matlab code:
	double besq = LengthSquared(b - b_dot_t * t);
	double loga = 2.0 * L * log(a / (L_alpha - L));
	double p = square(b_dot_t) * (a - L_alpha + loga) + besq / (1.0 - nu) * (2.0 * a - 2.0 * L_alpha + loga);
	double selfEnergy = mu * p / (8.0 * M_PI);
#endif

	SIMULATION_CHECK_VALUE(selfEnergy);
	return selfEnergy * 6.241506e-12;		// Convert from Pa*A^3 to eV
}

/******************************************************************************
* Calculates the elastic interaction energy of a (virtual) straight segment
* and the existing dislocations.
*****************************************************************************/
double DislocationNetwork::elasticSegmentNetworkEnergy(const Point3& x1, const Point3& x2)
{
	// Verify simulation parameters.
	SIMULATION_ASSERT(params().coreWidth > 0);
	SIMULATION_ASSERT(params().shearModulus > 0);
	SIMULATION_ASSERT(params().poissonRatio > 0);

	// Burgers vector in spatial coordinates.
	Vector3 b = params().unitCell * params().b;

	double energy = 0;

	// Loop over all existing dislocation segments.
	for(SegmentIterator segment(*this); !segment.atEnd(); ++segment) {

		Point3 x3 = params().unitCell * segment->node1()->pos();
		Point3 x4 = x3 + params().unitCell * segment->lineVector();

		// Compute interaction energy with the two closest images of the segment.
		double d = floor((x1.Z - x3.Z) / params().pbcLength);
		Vector3 shift1 = Vector3(0, 0, params().pbcLength * d);
		Vector3 shift2 = Vector3(0, 0, params().pbcLength * (d+1));

		energy += elasticSegmentSegmentEnergy(x1, x2, x3 + shift1, x4 + shift1, b);
		energy += elasticSegmentSegmentEnergy(x1, x2, x3 + shift2, x4 + shift2, b);
	}

	SIMULATION_CHECK_VALUE(energy);
	return energy;
}

/******************************************************************************
* Computes the elastic interaction energy between two straight dislocation segments.
* Does not include effects of periodic images.
*****************************************************************************/
double DislocationNetwork::elasticSegmentSegmentEnergy(const Point3& x1, const Point3& x2, const Point3& x3, const Point3& x4, const Vector3& b)
{
	Vector3 x12 = x2 - x1;
	Vector3 x34 = x4 - x3;
	double L_m = Length(x12);
	double L_n = Length(x34);
	Vector3 t = x12/L_m;
	Vector3 tprime = x34/L_n;
	Vector3 u = CrossProduct(t, tprime);
	if(fabs(LengthSquared(u)) > CAFLOAT_EPSILON) {
		// Non-parallel segments:

		Vector3 v = CrossProduct(u, t);
		Vector3 vprime = CrossProduct(tprime, u);
		double A1 = (1.0 + params().poissonRatio) * DotProduct(b, t) * DotProduct(b, tprime);
		double A2 = (LengthSquared(b) + square(DotProduct(b, t))) * DotProduct(t, tprime);
		double A2prime = (LengthSquared(b) + square(DotProduct(b, tprime))) * DotProduct(t, tprime);
		double A3 = 2.0 * DotProduct(b, u) * DotProduct(b, v) * DotProduct(t, tprime) / LengthSquared(u);
		double A3prime = 2.0 * DotProduct(b, u) * DotProduct(b, vprime) * DotProduct(t, tprime) / LengthSquared(u);
		double A4 = (DotProduct(b,t) * DotProduct(b, v) + DotProduct(b, tprime) * DotProduct(b, vprime)) * DotProduct(t, tprime);
		double A5 = 2.0 * LengthSquared(CrossProduct(b, u)) * DotProduct(t, tprime) / LengthSquared(u);
		SIMULATION_ASSERT(LengthSquared(u) != 0);
		SIMULATION_ASSERT(LengthSquared(v) != 0);

		// Define kernel function W_np(x) used in non-singular expression.
		auto W_np = [&](const Vector3& x) -> double {
			double asq = square(params().coreWidth * params().blength);
			double R_a = sqrt(LengthSquared(x) + asq);
			double W =
				log(R_a + DotProduct(x, tprime)) * DotProduct(x, ((A1 - A2prime)*vprime + A3prime*u)) +
				log(R_a + DotProduct(x, t)) * DotProduct(x, ((A1 - A2)*v + A3*u)) +
				A4 * R_a +
				(A1 - A5)*(2.0 * square(DotProduct(x,u)) + LengthSquared(u)*asq) / sqrt(square(DotProduct(x,u)) + LengthSquared(u)*asq) *
				atan(((1.0 + DotProduct(t,tprime))*R_a + DotProduct(x, t+tprime)) / sqrt(square(DotProduct(x,u)) + LengthSquared(u)*asq));
			SIMULATION_CHECK_VALUE(W);
			return W;
		};

		// Factor 6.241506e-12 is to convert from Pa*A^3 to eV.
		double prefactor = 6.241506e-12 * params().shearModulus / (4.0 * M_PI * (1.0 - params().poissonRatio) * LengthSquared(u));
		SIMULATION_CHECK_VALUE(prefactor);

		return prefactor * (W_np(x4 - x2) + W_np(x2 - x1) - W_np(x4 - x1) - W_np(x3 - x2));
	}
	else {
		// Parallel segments:

		// Define kernel function W_np(x) used in non-singular expression.
		double bmag = Length(b);
		auto W_p = [&](const Vector3& x) -> double {
			double asq = square(params().coreWidth * params().blength);
			double R_a = sqrt(LengthSquared(x) + asq);
			double W =
				(2.0 * bmag * DotProduct(b,x) - square(bmag)*DotProduct(t,x)*(3.0 - params().poissonRatio)) *
				log(R_a + DotProduct(x, t)) + R_a * square(bmag) * (2.0 - params().poissonRatio) -
				0.5 * R_a * (square(DotProduct(b, x) - bmag * DotProduct(t, x)) - asq*square(bmag)*(params().poissonRatio - 1.0)) /
				(square(R_a) - square(DotProduct(t, x)));
			SIMULATION_CHECK_VALUE(W);
			return W;
		};

		// Factor 6.241506e-12 is to convert from Pa*A^3 to eV.
		double prefactor = 6.241506e-12 * params().shearModulus / (4.0 * M_PI * (1.0 - params().poissonRatio));
		SIMULATION_CHECK_VALUE(prefactor);

		return prefactor * (W_p(x4 - x2) + W_p(x2 - x1) - W_p(x4 - x1) - W_p(x3 - x2));
	}
}

/******************************************************************************
* Computes the stress at a point due to the given dislocation segment.
*****************************************************************************/
SymmetricTensor2 DislocationNetwork::segmentStress(const Point3& p, const Point3& p1, const Point3& p2, const Vector3& b)
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
	double a = params().coreWidth * params().blength;
	double MU = params().shearModulus;
	double NU = params().poissonRatio;

	vec1x = p2.X - p1.X;
	vec1y = p2.Y - p1.Y;
	vec1z = p2.Z - p1.Z;

	double length = sqrt(vec1x*vec1x + vec1y*vec1y + vec1z*vec1z);
	if(length == 0) return SymmetricTensor2(0);

	oneoverLp = 1.0 / length;

	tpx = vec1x * oneoverLp;
	tpy = vec1y * oneoverLp;
	tpz = vec1z * oneoverLp;

	Rx = p.X - p1.X;
	Ry = p.Y - p1.Y;
	Rz = p.Z - p1.Z;

	Rdt = Rx*tpx + Ry*tpy + Rz*tpz;

	ndx = Rx - Rdt*tpx;
	ndy = Ry - Rdt*tpy;
	ndz = Rz - Rdt*tpz;

	d2 = ndx*ndx + ndy*ndy + ndz*ndz;

	s1 = -Rdt;
	s2 = -((p.X-p2.X)*tpx + (p.Y-p2.Y)*tpy + (p.Z-p2.Z)*tpz);
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


	txbx = tpy*b.Z - tpz*b.Y;
	txby = tpz*b.X - tpx*b.Z;
	txbz = tpx*b.Y - tpy*b.X;

	dxbx = ndy*b.Z - ndz*b.Y;
	dxby = ndz*b.X - ndx*b.Z;
	dxbz = ndx*b.Y - ndy*b.X;

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

	SymmetricTensor2 stress;
	stress.voigtVector()[0] = I_03xx*s_03 + I_13xx*s_13 + I_05xx*s_05 +
				I_15xx*s_15 + I_25xx*s_25;

	stress.voigtVector()[1] = I_03yy*s_03 + I_13yy*s_13 + I_05yy*s_05 +
				I_15yy*s_15 + I_25yy*s_25;

	stress.voigtVector()[2] = I_03zz*s_03 + I_13zz*s_13 + I_05zz*s_05 +
				I_15zz*s_15 + I_25zz*s_25;

	stress.voigtVector()[3] = I_03xy*s_03 + I_13xy*s_13 + I_05xy*s_05 +
				I_15xy*s_15 + I_25xy*s_25;

	stress.voigtVector()[4] = I_03yz*s_03 + I_13yz*s_13 + I_05yz*s_05 +
				I_15yz*s_15 + I_25yz*s_25;

	stress.voigtVector()[5] = I_03xz*s_03 + I_13xz*s_13 + I_05xz*s_05 +
				I_15xz*s_15 + I_25xz*s_25;

	return stress;
}
