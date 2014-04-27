#ifndef hohe_InitParticleDistributor_H
#define hohe_InitParticleDistributor_H

#include <hohe2Common/geo/basicGeos.h>
#include <sph/FluidSolverSimpleSph.h>

namespace hohehohe2
{

///Initial particle placement.
class InitParticleDistributor
{

public:

	///Create particles and set initial values to them. pos and velocity will be setSize() and allocate()d.
	static void set(PointSet& pos, PointSet& velocity, float restLength, unsigned int id);

private:

	///Single particle.
	static void placement0_(PointSet& pos, PointSet& velocity, float restLength);

	///2 particles.
	static void placement1_(PointSet& pos, PointSet& velocity, float restLength);

	///Small n x n x n.
	static void placement2_(PointSet& pos, PointSet& velocity, float restLength);

	///1 x n x 1.
	static void placement3_(PointSet& pos, PointSet& velocity, float restLength);

	///Large n x n x n.
	static void placement4_(PointSet& pos, PointSet& velocity, float restLength);

};

}

#endif
