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

	///Step the simulation.
	virtual void step(float deltaT);


private:

	void calcDensity_(Particles& particles);
	void applyForce_(Particles& particles);
	void project_(Particles& particles);
	void integrate_(Particles& particles);
};

}

#endif
