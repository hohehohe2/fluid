#ifndef hohe_InitParticleDistributor_H
#define hohe_InitParticleDistributor_H

#include <hohe2Common/geo/basicGeos.h>
#include <sph/FluidSolverSimpleSph.h>

namespace hohehohe2
{

struct ParticlesFluid;
struct ParticlesWall;

///Initial particle placement.
class InitParticleDistributor
{

public:

	///Create particles and set initial values to them. pos and velocity will be setSize() and allocate()d.
	static void set(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance, unsigned int id);

private:

	///Single particle.
	static void placement0_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance);

	///2 particles.
	static void placement1_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance);

	///Small n x n x n.
	static void placement2_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance);

	///1 x n x 1.
	static void placement3_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance);

	///Large n x n x n.
	static void placement4_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance);

	///Small n x n x n with ground particles.
	static void placement5_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance);

	///4 particles.
	static void placement6_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance);

	///Single particle with ground particles.
	static void placement7_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance);

};

}

#endif
