#include "FluidSolverSimpleSph.h"
#include <hohe2Common/util/BufferUtil.h>
#include <hohe2Common/container/CellCodeCalculator.h>
#include "Constants.h"
#include "ParticlesFluid.h"
#include "ParticlesWall.h"


using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
const float FluidSolverSimpleSph::PET_PEEVE_COURANT_NUMBER = 0.5f;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
FluidSolverSimpleSph::FluidSolverSimpleSph(float particleMass)
	:
	m_cHash(COMPACT_HASH_NUM_HASH_ENTRIES, COMPACT_HASH_NUM_ELEMENTS_IN_A_LIST, COMPACT_HASH_NUM_LISTS, HOST),
	m_cHashWall(COMPACT_HASH_NUM_HASH_ENTRIES, COMPACT_HASH_NUM_ELEMENTS_IN_A_LIST, COMPACT_HASH_NUM_LISTS, HOST),
	//m_pressureCalculator(particleMass),
	//m_pressurePciSphCalculator(particleMass, m_sphKernel, 0.01f, 0.01f, 4,
	m_pressurePciSphCalculator(particleMass, m_sphKernel, 0.01f, 0.01f, 1), //tako. for fluid stability.
	m_viscosityCalculator(particleMass), m_densityCalculator(particleMass)
{
	const unsigned int KERNEL_RADIUS_PER_EQUILIBRIUM_DISTANCE = 4;

	//Adjust the kernel radius so that several dozens of neighbor particles are in the radius at rest density.
	float particleVolume = particleMass / Constants::RO0;
	m_equilibriumDistance = (float)pow(particleVolume, 1.0 / 3.0);
	m_sphKernel.setKernelRadius(m_equilibriumDistance * KERNEL_RADIUS_PER_EQUILIBRIUM_DISTANCE);

	m_pressurePciSphCalculator.precompute(m_equilibriumDistance, KERNEL_RADIUS_PER_EQUILIBRIUM_DISTANCE);
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::step(ParticlesFluid& particles, ParticlesWall& particlesWall, float deltaT)
{
	float remaining = deltaT;
	bool loop = true;
	do
	{
		float maxVelocity = particles.m_velocity->calcMaxLength(HOST);
		float dt = PET_PEEVE_COURANT_NUMBER * m_sphKernel.r() / maxVelocity;
		if (dt > remaining)
		{
			dt = remaining;
			loop = false;
		}
		remaining -= dt;

		std::cout << remaining << ": updateNeighbors - ";
		CellCodeCalculator ccc;
		updateNeighbors_(particles, particlesWall, ccc);
		std::cout << "volume - ";
		m_volumeCalculator.calculation(particlesWall, m_sphKernel, ccc, m_cHashWall, HOST);
		std::cout << "density - ";
		m_densityCalculator.calculation(particles, m_sphKernel, ccc, m_cHash, HOST, &particlesWall, &m_cHashWall);
		initAcceleration_host_(particles);
		std::cout << "viscosity - ";
		//m_pressureCalculator.calculation(particles, m_sphKernel, ccc, m_cHash, HOST, &particlesWall, &m_cHashWall);
		m_viscosityCalculator.calculation(particles, m_sphKernel, ccc, m_cHash, HOST);
		std::cout << "pressure - ";
		m_pressurePciSphCalculator.calculation(particles, ccc, m_cHash, HOST);
		std::cout << "integrate\n";
		m_semiImplicitEulerIntegrateCalculator.integrate(particles, dt, HOST);

	} while(loop);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::updateNeighbors_(ParticlesFluid& particles, ParticlesWall& particlesWall, CellCodeCalculator& ccc)
{

	//---- Define spatial grid cells and create CellCodeCalculator.
	const float kernelRaidus = m_sphKernel.r();

	particles.sync(HOST);

	BoundingBox bbox;
	particles.m_pos->calcBoundingBox(bbox, HOST);
	BoundingBox bboxWall;
	particlesWall.m_pos->calcBoundingBox(bboxWall, HOST);
	bbox.makeUnion(bboxWall);

	//Shift min values a little bit to make sure the particles which have min values are inside the cell.
	bbox.m_min -= Point(kernelRaidus, kernelRaidus, kernelRaidus) / 100.0f;

	ccc.reset(bbox.m_min, kernelRaidus);

	//---- Create hash for the fluid particles.
	//Get the unsorted code list.
	BufferUInt codeSet;
	ccc.getCode32(codeSet, *particles.m_pos, HOST);

	//Sort it. Once we create the hash, we don't have to keep using it so it's not a particles member variable.
	BufferUInt sortedCodeSet;
	BufferUtil::sortByKey(*particles.m_sortedIdMap, sortedCodeSet, codeSet, DEVICE);

	//Create compact hash.
	if ( ! m_cHash.build(sortedCodeSet, HOST))
	{
		std::cerr << "Fluid hash full\n";
	}

	//----Do the same for the wall particles.
	ccc.getCode32(codeSet, *particlesWall.m_pos, HOST);
	BufferUtil::sortByKey(*particlesWall.m_sortedIdMap, sortedCodeSet, codeSet, DEVICE);
	if ( ! m_cHashWall.build(sortedCodeSet, HOST))
	{
		std::cerr << "Wall hash full\n";
	}

}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::initAcceleration_host_(ParticlesFluid& particles)
{
	float* axs = particles.m_acceleration->xs(HOST);
	float* ays = particles.m_acceleration->ys(HOST);
	float* azs = particles.m_acceleration->zs(HOST);

	unsigned int size = particles.size();

	for (int idP = 0; idP < (int)size; ++idP)
	{
		//----Gravity.
		axs[idP] = 0.0f;
		ays[idP] = Constants::G;
		azs[idP] = 0.0f;
	}

	particles.setClean(HOST);
}
