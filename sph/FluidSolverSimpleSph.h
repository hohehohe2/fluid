#ifndef hohe_FluidSolverSimpleSph_H
#define hohe_FluidSolverSimpleSph_H

#include <hohe2Common/cuda/Buffer.h>
#include <hohe2Common/geo/basicGeos.h>
#include <hohe2Common/container/CompactHash.h>
#include "SphKernel.h"
#include "PressureCalculator.h"


namespace hohehohe2
{

struct FluidParticles;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Simple SPH.
/**
Physical quantities are measured in MKS, not normalized here.
**/
class FluidSolverSimpleSph
{

public:

	///Constructor.
	/**
	@param particleMass in MKS.
	**/
	FluidSolverSimpleSph(float particleMass=1.0f);

	///Destructor.
	~FluidSolverSimpleSph(){}

	///Step the simulation.
	void step(FluidParticles& particles, float deltaT);

	///Get the distance between particles at rest density.
	float restLength() const {return m_restLength;}

private:

	///Distance between particles at rest density, just for initial distribution hint. (volume per particle) ^ (1/3).
	float m_restLength;

	///Mass per single particle.
	float m_particleMass;

	///SPH kernel calculator.
	SphKernel m_sphKernel;

	///Compact hash for neighbor search.
	CompactHash m_cHash;

	PressureCalculator m_pressureCalculator;

private:

	float calcMaxVelocity_(const FluidParticles& particles);
	void updateNeighbors_(FluidParticles& particles);
	void calcDensity_host_(FluidParticles& particles);
	void calcAcceleration_host_(FluidParticles& particles);
	void integrate_(FluidParticles& particles, float deltaT);

	float densityToPressure_(float density);

	///Density->pressure stiffness coefficient.
	/**
	It is not in Constants.h because this is an artificially enough physical value.
	**/
	static const float K;

	///Viscosity coefficient.
	static const float MU;

	///So called Courant number (pet peeve for Prof. Bridson ;).
	static const float PET_PEEVE_COURANT_NUMBER;

	//Compact hash parameters. Need adjustment.
	static const unsigned int COMPACT_HASH_NUM_HASH_ENTRIES = 2048;
	static const unsigned int COMPACT_HASH_NUM_ELEMENTS_IN_A_LIST = 256;
	static const unsigned int COMPACT_HASH_NUM_LISTS = 1024;

};

}

#endif
