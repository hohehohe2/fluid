#ifndef hohe_FluidSolver_H
#define hohe_FluidSolver_H


namespace hohehohe2
{

class Particles;

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Base class of every fluid solver.
class FluidSolver
{

public:

	///Constructor.
	FluidSolver(){}

	///Destructor.
	virtual ~FluidSolver(){}

	///Step the simulation.
	virtual void step(Particles& particles, float deltaT) = 0;
};

}

#endif
