#ifndef hohe_FluidSolverSimpleSph_H
#define hohe_FluidSolverSimpleSph_H

#include "FluidSolver.h"
#include <geo/Particles.h>

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Simple SPH.
class FluidSolverSimpleSph : public FluidSolver
{

public:

	///Constructor.
	FluidSolverSimpleSph(){}

	///Destructor.
	virtual ~FluidSolverSimpleSph(){}

	///Initialize this solver.
	void initSolver();

	///Step the simulation.
	virtual void step(float deltaT);

private:

	///Particles to solve.
	Particles::SPtr m_particles;

	///SPH kernel radius.
	float m_kernelRadius;

private:

	void updateNeighbors_();
	void applyForce_();
	void calcDensity_();
	void project_();
	void integrate_(float deltaT);
	float densityToPressure_(float ro);
};

}

#endif
