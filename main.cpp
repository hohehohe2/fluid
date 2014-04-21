#include "sph/FluidSolverSimpleSph.h"

using namespace hohehohe2;


int main()
{
	FluidSolverSimpleSph ssph;
	FluidSolverSimpleSph::Particles* particles = FluidSolverSimpleSph::Particles::createInstance(256);
	BufferSet::SPtr sptr = particles->getSelfSptr();
	ssph.step(*particles, 0.1f);
	return 0;
}
