#ifndef hohe_FluidSolver_H
#define hohe_FluidSolver_H


namespace hohehohe2
{

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
	virtual void step(float deltaT) = 0;
};

}

#endif
