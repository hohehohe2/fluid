#ifndef hohe_FluidSolverSimpleSph_H
#define hohe_FluidSolverSimpleSph_H

#include <hohe2Common/cuda/Buffer.h>
#include <hohe2Common/geo/basicGeos.h>
#include <hohe2Common/container/CompactHash.h>
#include "SphKernel.h"
#include "CalculatorVolume.h"
#include "CalculatorDensity.h"
#include "CalculatorPressure.h"
#include "CalculatorViscosity.h"
#include "CalculatorSemiImplicitEulerIntegrate.h"
#include "CalculatorPressurePciSph.h"

namespace hohehohe2
{

struct ParticlesFluid;
struct ParticlesWall;
class CellCodeCalculator;

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
	void step(ParticlesFluid& particles, ParticlesWall& particlesWall, float deltaT);

	///Get the distance between particles at rest density.
	float equilibriumDistance() const {return m_equilibriumDistance;}

private:

	///Distance between particles at rest density. (volume per particle) ^ (1/3).
	float m_equilibriumDistance;

	///SPH kernel calculator.
	SphKernel m_sphKernel;

	///Compact hash for neighbor fluid particle search.
	CompactHash m_cHash;

	///Compact hash for neighbor wall particle search.
	CompactHash m_cHashWall;

	CalculatorVolume m_volumeCalculator;

	CalculatorDensity m_densityCalculator;

	//CalculatorPressure m_pressureCalculator;

	CalculatorPressurePciSph m_pressurePciSphCalculator;

	CalculatorViscosity m_viscosityCalculator;

	CalculatorSemiImplicitEulerIntegrate m_semiImplicitEulerIntegrateCalculator;

private:

	void updateNeighbors_(ParticlesFluid& particles, ParticlesWall& particlesWall, CellCodeCalculator& ccc);
	void initAcceleration_host_(ParticlesFluid& particles);

	///So called Courant number (pet peeve for Prof. Bridson ;).
	static const float PET_PEEVE_COURANT_NUMBER;

	//Compact hash parameters. Need adjustment.
	static const unsigned int COMPACT_HASH_NUM_HASH_ENTRIES = 2048;
	static const unsigned int COMPACT_HASH_NUM_ELEMENTS_IN_A_LIST = 256;
	static const unsigned int COMPACT_HASH_NUM_LISTS = 1024;

};

}

#endif
