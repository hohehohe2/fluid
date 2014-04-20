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
	std::cout << "000\n";
	*lines2 = *lines1;
	std::cout << "\n---\n";
	sptr1.reset();
	std::cout << "\n===\n";

	return 0;
}
