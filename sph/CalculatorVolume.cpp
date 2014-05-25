#include "CalculatorVolume.h"

#include "Constants.h"
#include <hohe2Common/container/CellCodeCalculator.h>
#include <hohe2Common/container/CompactHash.h>
#include "ParticlesWall.h"
#include "FluidSolverSimpleSph.h"

using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void CalculatorVolume::calculation_host_(ParticlesWall& particles, const GlobalFluidParameters& globalParam, const CellCodeCalculator& ccc, const CompactHash& cHash)
{
	m_sphKernelPoly6.setKernelRadius(globalParam.m_kernelRadius);

	particles.m_pos->sync(HOST);
	particles.m_sortedIdMap->sync(HOST);
	const float* pxs = particles.m_pos->xs(HOST);
	const float* pys = particles.m_pos->ys(HOST);
	const float* pzs = particles.m_pos->zs(HOST);
	const unsigned int* sortedIdMaps = particles.m_sortedIdMap->get(HOST);
	float* vs = particles.m_volume->get(HOST);

	const unsigned int size = particles.size();

	#pragma omp parallel for
	for (int idP = 0; idP < (int)size; ++idP)
	{
		float sumW = 0.0f;

		for (unsigned int i= 0; i < 27; ++i)
		{
			bool isFilled;
			const unsigned int code = ccc.getNeighborCode32(isFilled, pxs[idP], pys[idP], pzs[idP], i);
			if ( ! isFilled)
			{
				continue;
			}
			unsigned int index;
			const unsigned int numObjects = cHash.lookup(index, code);
			for (unsigned int j = 0; j < numObjects; ++j)
			{
				const unsigned int idN = sortedIdMaps[index + j];
				const float distx = pxs[idN] - pxs[idP];
				const float disty = pys[idN] - pys[idP];
				const float distz = pzs[idN] - pzs[idP];
				const float dist2 = distx * distx + disty * disty + distz * distz;
				sumW += m_sphKernelPoly6.wPart(dist2);
			}
		}
		sumW *= m_sphKernelPoly6.getConstant();

		//To be implemented.
		//Need to adjust the coefficient so that a particle which distance from the wall is the same as the length between particles at the rest density
		//gets the force from the wall which is equivalent from the force where the wall particle is uniformly distributed, i.e.
		//
		//     (A)              (B)
		//
		// ...........      ...........
		// ...........      ...........
		// .....o.....      .....o.....
		// XXXXXXXXXXX      xxxxxxxxxxx
		//                  xxxxxxxxxxx
		//                  xxxxxxxxxxx
		//
		// o is the fluid partcle we are interested in
		// . are the neighbor fluid particles
		// X are the weighted wall particles
		// x are the un-weighted wall particles (i.e. same influence as the neighbor fluid particles)
		//
		//Here (A) is Akinci2012 and (B) is the traditional multi-layer wall particle model.
		//We need to adjust the scalingFactor so that
		//sum of the forces applied to o from 'X's in (A) == sum of the forces applied to o from 'x's in (B)
		//According to the paper it increases the volume of a wall particle scales by a factor of 1.4 but
		//in my condition it's about 3 (can be solved analitically but I just tested with poly6 and cubic spline kernels).
		//I don't know what makes the difference yet since a factor of 3 seems to be natural to me.

		const float scalingFactor = 1.267f;
		vs[idP] = scalingFactor / sumW;
	}

	particles.setClean(HOST);
}
