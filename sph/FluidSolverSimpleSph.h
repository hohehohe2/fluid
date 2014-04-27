#ifndef hohe_FluidSolverSimpleSph_H
#define hohe_FluidSolverSimpleSph_H

#include <hohe2Common/cuda/Buffer.h>
#include <hohe2Common/geo/basicGeos.h>
#include <hohe2Common/container/CompactHash.h>
#include "SphKernel.h"
#include "PressureCalculator.h"
#include "ViscosityCalculator.h"
#include "DensityCalculator.h"


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

	///SPH kernel calculator.
	SphKernel m_sphKernel;

	///Compact hash for neighbor search.
	CompactHash m_cHash;

	PressureCalculator m_pressureCalculator;

	ViscosityCalculator m_viscosityCalculator;

	DensityCalculator m_densityCalculator;

private:

	void updateNeighbors_(FluidParticles& particles);
	void initAcceleration_host_(FluidParticles& particles);
	void integrate_(FluidParticles& particles, float deltaT);

	///So called Courant number (pet peeve for Prof. Bridson ;).
	static const float PET_PEEVE_COURANT_NUMBER;

	//Compact hash parameters. Need adjustment.
	static const unsigned int COMPACT_HASH_NUM_HASH_ENTRIES = 2048;
	static const unsigned int COMPACT_HASH_NUM_ELEMENTS_IN_A_LIST = 256;
	static const unsigned int COMPACT_HASH_NUM_LISTS = 1024;

};

}

#endif
