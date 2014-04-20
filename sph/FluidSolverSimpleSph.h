#ifndef hohe_FluidSolverSimpleSph_H
#define hohe_FluidSolverSimpleSph_H

#include "FluidSolver.h"
#include "SphKernel.h"
#include <geo/Particles.h>

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Simple SPH.
/**
Physical quantities are measured in MKS, not normalized here.
**/
class FluidSolverSimpleSph : public FluidSolver
{

public:

	///Constructor.
	/**
	@param particleMass in MKS.
	**/
	FluidSolverSimpleSph(float particleMass=1.0);

	///Destructor.
	virtual ~FluidSolverSimpleSph(){}

	///Setup a particles for this solver.
	/**
	After colling this, particle pos and velocity must be filled before the first step.
	**/
	void setupParticles(Particles& particles, unsigned int numPartices);

	///Step the simulation.
	virtual void step(Particles& particles, float deltaT);

private:

	///Mass per single particle.
	float m_particleMass;

	///SPH kernel calculator.
	SphKernel m_sphKernel;

private:

	float calcMaxVelocity_(Particles& particles);
	void updateNeighbors_(Particles& particles);
	void calcDensity_(Particles& particles);
	void calcAcceleration_(Particles& particles);
	void integrate_(Particles& particles, float deltaT);

	float densityToPressure_(float density);

	///Density->pressure stiffness coefficient.
	/**
	It is not in Constants.h because this is an artificially enough physical value.
	**/
	static const float K;

    ///So called Courant number (pet peeve for Prof. Bridson ;).
    static const float C;

};

}

#endif
