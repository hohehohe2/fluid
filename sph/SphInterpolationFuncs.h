#ifndef hohe_SphInterpolationFuncs_H
#define hohe_SphInterpolationFuncs_H

namespace hohehohe2
{

class SphKernel;
class CellCodeCalculator;
class CompactHash;


///A bunch of SPH interpolation functions. It assumes seaching inside the nearby 27 cells is enough for neighbor particle search.
struct SphInterpolationFuncs
{

	///Sigma W -> result.
	static void sumW_host_(
		float& result,
		unsigned int idP, const float* pxs, const float* pys, const float* pzs,
		const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash);

	///Sigma gradW -> result[3].
	static void sumGradW_host_(
		float* result,
		unsigned int idP, const float* pxs, const float* pys, const float* pzs,
		const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash);

	///Sigma (gradW dot gradW) -> result.
	static void sumGradWDotGradW_host_(
		float& result,
		unsigned int idP, const float* pxs, const float* pys, const float* pzs,
		const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash);

	///Calculate the SPH interpolation of a physical value. Sigma (value[N] * W) -> result.
	static void sumValW_host_(
		float& result,
		unsigned int idP, const float* pxs, const float* pys, const float* pzs, const float* neighborScalarVals,
		const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash);

	///Calculate the SPH interpolation of the gradient of a physical value. Sigma (value[N] * gradW) -> result[3].
	static void sumValGradW_host_(
		float* result,
		unsigned int idP, const float* pxs, const float* pys, const float* pzs, const float* neighborScalarVals,
		const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash);

	///Calculate the SPH interpolation of the laplacian of a physical value. Sigma (value[N] * laplacianW) -> result.
	static void sumValLaplaceW_host_(
		float& result,
		unsigned int idP, const float* pxs, const float* pys, const float* pzs, const float* neighborScalarVals,
		const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash);

};

}

#endif
