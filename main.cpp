#include <cudaCommon/defines.h>

#include "sph/FluidSolverSimpleSph.h"
#include "geo/basicGeos.h"

using namespace hohehohe2;


int main()
{
	FluidSolverSimpleSph ssph;
	Lines* lines1 = new Lines;
	Lines* lines2 = new Lines;
	std::shared_ptr < BufferSet > sptr1 = lines1->getSelfSptr();
	std::shared_ptr < BufferSet > sptr2 = lines2->getSelfSptr();
	*lines2 = *lines1;
	sptr1.reset();

	return 0;
}
