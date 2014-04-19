#include "FluidSolverSimpleSph.h"
#include "Constants.h"

using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::initSolver()
{
	//TODO H: Impl.

	//m_kernelRadius = 3.0f;

	///Density->pressure stiffness coefficient.
	//float m_k;

}



//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::step(float deltaT)
{
	updateNeighbors_();
	calcDensity_();
	applyForce_();
	project_();
	integrate_(deltaT);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::updateNeighbors_()
{
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::calcDensity_()
{
	float* xs = m_particles->m_pos->m_xs->get(true);
	float* ys = m_particles->m_pos->m_ys->get(true);
	float* zs = m_particles->m_pos->m_zs->get(true);
	//float* ds = m_particles->m_de->m_zs(true);

	unsigned int size = m_particles->size();
	for (unsigned int pid = 0; pid < size; ++pid)
	{
		for (unsigned int nid = pid; nid < size; ++nid)
		{
			float diffx = xs[nid] - xs[pid];
			float diffy = ys[nid] - ys[pid];
			float diffz = zs[nid] - zs[pid];
			float lengthSquared = diffx * diffx + diffy * diffy + diffz * diffz;
			if (lengthSquared > m_kernelRadius * m_kernelRadius)
			{
				continue;
			}

			
		}
	}
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::applyForce_()
{
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::project_()
{
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::integrate_(float deltaT)
{
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
float FluidSolverSimpleSph::densityToPressure_(float ro)
{
	return Constants::K * (ro * ro / (Constants::RO0 * Constants::RO0) - 1.0f);

}
