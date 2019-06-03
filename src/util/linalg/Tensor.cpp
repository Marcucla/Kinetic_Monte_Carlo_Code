#include "../../StandardIncludes.h"

using namespace std;

/******************************************************************************
* Computes the eigenvalues and optionally eigenvectors of this symmetric tensor.
******************************************************************************/
void SymmetricTensor2::eigenvalues(CAFloat lambda[3], Matrix3* eigenvectors) const
{
	if(eigenvectors) *eigenvectors = IDENTITY;

	SymmetricTensor2 T2(*this);
	if(T2.a[3]*T2.a[3]+T2.a[4]*T2.a[4]+T2.a[5]*T2.a[5] > 1e3 * std::numeric_limits<CAFloat>::min()) {

		CAFloat delta = 1.e-20 * DoubleContraction(T2,T2); 
		for(size_t iter = 0; iter < 5; iter++) {
			// For large diagonal terms and small but nonzero off-diagonals,
			// do at least one sweep.
			T2.Givens(1, 0, eigenvectors);
			T2.Givens(2, 0, eigenvectors);
			T2.Givens(2, 1, eigenvectors);

			if(T2.a[3]*T2.a[3]+T2.a[4]*T2.a[4]+T2.a[5]*T2.a[5] <= delta)
				break;
		}
	}

	lambda[0] = T2.a[0];
	lambda[1] = T2.a[1];
	lambda[2] = T2.a[2];
}

/******************************************************************************
* Finds the maximum absolute eigenvalue of the symmetric tensor.
******************************************************************************/
CAFloat SymmetricTensor2::maxEigenvalue() const
{
	using namespace std;

	CAFloat lambda[3];
	eigenvalues(lambda);

	// Find maximum eigenvalue.
	CAFloat dmax = 0.0;
	for(size_t i=0; i<3; i++)
		dmax = max(dmax, fabs(lambda[i]));
	return dmax;
}

/******************************************************************************
* Finds the minimum absolute eigenvalue of the symmetric tensor.
******************************************************************************/
CAFloat SymmetricTensor2::minEigenvalue() const
{
	using namespace std;

	CAFloat lambda[3];
	eigenvalues(lambda);

	// Find minimum eigenvalue.
	CAFloat dmin = CAFLOAT_MAX;
	for(size_t i=0; i<3; i++)
		dmin = min(dmin, fabs(lambda[i]));
	return dmin;
}

/******************************************************************************
* Perform Givens rotation of this symmetric tensor
******************************************************************************/
void SymmetricTensor2::Givens(size_t p, size_t q, Matrix3* eigenvectors)
{ 
	using namespace std;

	/* Givens transformation for row and column p and q
	See Matrix Computations, Golub & Van Loan, 2nd ed.,
	section 5.1.8 and section 8.5.2 */

	CAFloat Apq, tau, t, s, c, taup, tauq;

	SIMULATION_ASSERT(p<3);
	SIMULATION_ASSERT(q<3);
	SIMULATION_ASSERT(p!=q);

	if(p<q) swap(p,q);

	Apq = (*this)(p,q);
	if(fabs(Apq) > 1e3 * numeric_limits<CAFloat>::min()) {
		CAFloat App = a[p]; //A(p,p)
		CAFloat Aqq = a[q]; //A(q,q)

		tau = (Aqq - App)/(2*Apq);
		t = 1/(fabs(tau)+sqrt(1+tau*tau));
		if(tau < 0) t = -t;
		c = 1/sqrt(1+t*t);
		s = t*c;

		CAFloat Apq2cs = Apq*2*c*s;
		CAFloat c2 = c*c;
		CAFloat s2 = s*s;

		// Pre multiply by J.transpose() and post multiply by J.
		a[p] = App*c2 + Aqq*s2 - Apq2cs; // A(p,p) 
		a[q] = Aqq*c2 + App*s2 + Apq2cs; // A(q,q) 
		(*this)(p,q) = 0;

		if(q==0)
			if(p==1) {
				CAFloat Arp = a[4]; //A(3,2)
				CAFloat Arq = a[5]; //A(3,1)
				a[4] = Arp*c-Arq*s;
				a[5] = Arp*s+Arq*c;
			}
			else { // p==2
				CAFloat Arp = a[4]; //A(3,2)
				CAFloat Arq = a[3]; //A(2,1)
				a[4] = Arp*c-Arq*s;
				a[3] = Arp*s+Arq*c;
			}
		else { // p==2 && q==1
			CAFloat Arp = a[5]; //A(3,1)
			CAFloat Arq = a[3]; //A(2,1)
			a[5] = Arp*c-Arq*s;
			a[3] = Arp*s+Arq*c;
		}

		// Post multiply 'eigenvectors' matrix by J to update eigen vectors (as columns).
		if(eigenvectors) {
			taup = (*eigenvectors)(0,p); tauq = (*eigenvectors)(0,q);
			(*eigenvectors)(0,p) = c*taup-s*tauq; (*eigenvectors)(0,q) = s*taup+c*tauq;
			taup = (*eigenvectors)(1,p); tauq = (*eigenvectors)(1,q);
			(*eigenvectors)(1,p) = c*taup-s*tauq; (*eigenvectors)(1,q) = s*taup+c*tauq;
			taup = (*eigenvectors)(2,p); tauq = (*eigenvectors)(2,q);
			(*eigenvectors)(2,p) = c*taup-s*tauq; (*eigenvectors)(2,q) = s*taup+c*tauq;
		}
	}
}
