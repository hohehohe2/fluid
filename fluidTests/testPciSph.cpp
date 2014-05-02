#include "gtest/gtest.h"

#include <hohe2Common/container/CompactHash.h>
#include <hohe2Common/container/CellCodeCalculator.h>
#include <hohe2Common/util/BufferUtil.h>
#include "sph/Constants.h"
#include "sph/CalculatorPressurePciSph.h"
#include "sph/ParticlesFluid.h"
#include "sph/SphKernel.h"

#include <iostream>

using namespace hohehohe2;


TEST(PciSph, create)
{
	const int kernelRadiusPerEquilibriumDistance = 4;
	float equilibriumDistance = 0.1f;
	SphKernel sphKernel;
	float deltaT = 0.1f;
	float particleMass = 1.0f;
	float kernelRadius = equilibriumDistance * kernelRadiusPerEquilibriumDistance;
	sphKernel.setKernelRadius(kernelRadius);
	float maxRelativeDensityError = 0.01f;
	unsigned int numMaxIterations = 6;

	//--------------Precompute.
	CalculatorPressurePciSph pcSph(particleMass, sphKernel, deltaT, maxRelativeDensityError, numMaxIterations);
	pcSph.precompute(equilibriumDistance, kernelRadiusPerEquilibriumDistance);

	//--------------Create particles.

	const unsigned int numLines = kernelRadiusPerEquilibriumDistance * 2 + 1;
	const unsigned int numParticles = numLines * numLines * numLines;
	const unsigned int centerPerticleId = numLines * numLines * kernelRadiusPerEquilibriumDistance + numLines * kernelRadiusPerEquilibriumDistance + kernelRadiusPerEquilibriumDistance;

	ParticlesFluid* particles = ParticlesFluid::createInstance(numParticles, HOST);
	float* pxs = particles->m_pos->xs(HOST);
	float* pys = particles->m_pos->ys(HOST);
	float* pzs = particles->m_pos->zs(HOST);

	particles->m_acceleration->memset(0, HOST);

	unsigned int id = 0;
	for (int i = -kernelRadiusPerEquilibriumDistance; i <= kernelRadiusPerEquilibriumDistance; ++i)
	{
		for (int j = -kernelRadiusPerEquilibriumDistance; j <= kernelRadiusPerEquilibriumDistance; ++j)
		{
			for (int k = -kernelRadiusPerEquilibriumDistance; k <= kernelRadiusPerEquilibriumDistance; ++k)
			{
				if (id == centerPerticleId)
				{
					//pxs[id] = 0.1f;
					//pys[id] = 0.2f;
					//pzs[id] = 0.3f;
					pxs[id] = 0.1f;
					pys[id] = 0.0f;
					pzs[id] = 0.0f;
				}
				else
				{
					pxs[id] = (float)i * equilibriumDistance;
					pys[id] = (float)j * equilibriumDistance;
					pzs[id] = (float)k * equilibriumDistance;
				}
				++id;
			}
		}
	}

	//--------------Create compact hash.

	CompactHash cHash(2048, 256, 1024, HOST);

	BoundingBox bbox;
	particles->m_pos->calcBoundingBox(bbox, HOST);

	//Shift min values a little bit to make sure the particles which have min values are inside the cell.
	bbox.m_min -= Point(kernelRadius, kernelRadius, kernelRadius) / 100.0f;

	CellCodeCalculator ccc;
	ccc.reset(bbox.m_min, kernelRadius);

	//Get the unsorted code list.
	BufferUInt codeSet;
	ccc.getCode32(codeSet, *particles->m_pos, HOST);

	//Sort it. Once we create the hash, we don't have to keep using it so it's not a particles member variable.
	BufferUInt sortedCodeSet;
	BufferUtil::sortByKey(*particles->m_sortedIdMap, sortedCodeSet, codeSet, DEVICE);

	if ( ! cHash.build(sortedCodeSet, HOST))
	{
		std::cerr << "Fluid hash full\n";
	}

	//--------------PCISPH.

	pcSph.calculation(*particles, sphKernel, ccc, cHash, HOST);

	//--------------Time integrate.

	float* axs = particles->m_acceleration->xs(HOST);
	float* ays = particles->m_acceleration->ys(HOST);
	float* azs = particles->m_acceleration->zs(HOST);

	float* vxs = particles->m_velocity->xs(HOST);
	float* vys = particles->m_velocity->ys(HOST);
	float* vzs = particles->m_velocity->zs(HOST);

	vxs[centerPerticleId] += axs[centerPerticleId] * deltaT;
	vys[centerPerticleId] += ays[centerPerticleId] * deltaT;
	vzs[centerPerticleId] += azs[centerPerticleId] * deltaT;
	pxs[centerPerticleId] += vxs[centerPerticleId] * deltaT;
	pys[centerPerticleId] += vys[centerPerticleId] * deltaT;
	pzs[centerPerticleId] += vzs[centerPerticleId] * deltaT;

	std::cout
		<< particles->m_pos->xs(HOST)[centerPerticleId] << " "
		<< particles->m_pos->ys(HOST)[centerPerticleId] << " "
		<< particles->m_pos->zs(HOST)[centerPerticleId] << " "
		<< std::endl;
	delete particles;
}
