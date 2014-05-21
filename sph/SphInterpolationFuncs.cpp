//#include "SphInterpolationFuncs.h"
//
//#include "Constants.h"
//#include <hohe2Common/container/CellCodeCalculator.h>
//#include <hohe2Common/container/CompactHash.h>
//#include <hohe2Common/geo/basicGeos.h>
//#include "SphKernel.h"
//
//
//using namespace hohehohe2;
//
//
////Neighbor iteration.
//namespace
//{
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//template < typename Operator , typename ReturnType >
//void neighborIterCalculation_host_(
//	ReturnType* result,
//	unsigned int idP, const float* pxs, const float* pys, const float* pzs,
//	const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash)
//{
//	for (unsigned int i= 0; i < 27; ++i)
//	{
//		bool isFilled;
//		const unsigned int code = ccc.getNeighborCode32(isFilled, pxs[idP], pys[idP], pzs[idP], i);
//		if ( ! isFilled)
//		{
//			continue;
//		}
//		unsigned int index;
//		const unsigned int numObjects = cHash.lookup(index, code);
//		for (unsigned int j = 0; j < numObjects; ++j)
//		{
//			const unsigned int idN = sortedIdMapsPtr[index + j];
//			Operator::calculate(result, idP, idN, pxs, pys, pzs, sphKernel);
//		}
//	}
//}
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//struct OperatorSumW_
//{
//    static void calculate(float* result, unsigned int idP, unsigned int idN, const float* pxs, const float* pys, const float* pzs, const SphKernel& sphKernel)
//	{
//		const float distx = pxs[idN] - pxs[idP];
//		const float disty = pys[idN] - pys[idP];
//		const float distz = pzs[idN] - pzs[idP];
//		const float dist2 = distx * distx + disty * disty + distz * distz;
//		*result += sphKernel.w(dist2);
//	}
//};
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//struct OperatorSumGradW_
//{
//    static void calculate(float* result, unsigned int idP, unsigned int idN, const float* pxs, const float* pys, const float* pzs, const SphKernel& sphKernel)
//	{
//		float gradW[3];
//		sphKernel.gradW(gradW, pxs[idP], pys[idP], pzs[idP], pxs[idN], pys[idN], pzs[idN]);
//		result[0] += gradW[0];
//		result[1] += gradW[1];
//		result[2] += gradW[2];
//	}
//};
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//struct OperatorSumGradWDotGradW_
//{
//    static void calculate(float* result, unsigned int idP, unsigned int idN, const float* pxs, const float* pys, const float* pzs, const SphKernel& sphKernel)
//	{
//		float gradW[3];
//		sphKernel.gradW(gradW, pxs[idP], pys[idP], pzs[idP], pxs[idN], pys[idN], pzs[idN]);
//		*result = gradW[0] * gradW[0] + gradW[1] * gradW[1] + gradW[2] * gradW[2];
//	}
//};
//
//
//}
//
//
////Neighbor scalar interpolation.
//namespace 
//{
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//template < typename Operator , typename ReturnType >
//void sphScalarInterpolation_host_(
//	ReturnType* result,
//	unsigned int idP, const float* pxs, const float* pys, const float* pzs, const float* neighborScalarVals,
//	const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash)
//{
//	for (unsigned int i= 0; i < 27; ++i)
//	{
//		bool isFilled;
//		const unsigned int code = ccc.getNeighborCode32(isFilled, pxs[idP], pys[idP], pzs[idP], i);
//		if ( ! isFilled)
//		{
//			continue;
//		}
//		unsigned int index;
//		const unsigned int numObjects = cHash.lookup(index, code);
//		for (unsigned int j = 0; j < numObjects; ++j)
//		{
//			const unsigned int idN = sortedIdMapsPtr[index + j];
//			Operator::calculate(result, idP, idN, pxs, pys, pzs, neighborScalarVals, sphKernel);
//		}
//	}
//}
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//struct OperatorSumWTimesNSclar_
//{
//    static void calculate(float* result, unsigned int idP, unsigned int idN, const float* pxs, const float* pys, const float* pzs, const float* neighborScalarVals, const SphKernel& sphKernel)
//	{
//		const float distx = pxs[idN] - pxs[idP];
//		const float disty = pys[idN] - pys[idP];
//		const float distz = pzs[idN] - pzs[idP];
//		const float dist2 = distx * distx + disty * disty + distz * distz;
//		*result += neighborScalarVals[idN] * sphKernel.w(dist2);
//	}
//};
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//struct OperatorSumGradWTimesNSclar_
//{
//    static void calculate(float* result, unsigned int idP, unsigned int idN, const float* pxs, const float* pys, const float* pzs, const float* neighborScalarVals, const SphKernel& sphKernel)
//	{
//		float gradW[3];
//		sphKernel.gradW(gradW, pxs[idP], pys[idP], pzs[idP], pxs[idN], pys[idN], pzs[idN]);
//		result[0] += neighborScalarVals[idN] * gradW[0];
//		result[1] += neighborScalarVals[idN] * gradW[1];
//		result[2] += neighborScalarVals[idN] * gradW[2];
//	}
//};
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//struct OperatorSumLaplaceWTimesNSclar_
//{
//    static void calculate(float* result, unsigned int idP, unsigned int idN, const float* pxs, const float* pys, const float* pzs, const float* neighborScalarVals, const SphKernel& sphKernel)
//	{
//		const float distx = pxs[idN] - pxs[idP];
//		const float disty = pys[idN] - pys[idP];
//		const float distz = pzs[idN] - pzs[idP];
//		const float dist2 = distx * distx + disty * disty + distz * distz;
//		*result += neighborScalarVals[idN] * sphKernel.laplaceW(dist2);
//	}
//};
//
//}
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//void SphInterpolationFuncs::sumW_host_(
//	float& result,
//	unsigned int idP, const float* pxs, const float* pys, const float* pzs,
//	const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash)
//{
//	result = 0;
//	neighborIterCalculation_host_ < OperatorSumW_ , float > (&result, idP, pxs, pys, pzs, sortedIdMapsPtr, sphKernel, ccc, cHash);
//}
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//void SphInterpolationFuncs::sumGradW_host_(
//	float* result,
//	unsigned int idP, const float* pxs, const float* pys, const float* pzs,
//	const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash)
//{
//	result[0] = 0;
//	result[1] = 0;
//	result[2] = 0;
//	neighborIterCalculation_host_ < OperatorSumGradW_ , float > (result, idP, pxs, pys, pzs, sortedIdMapsPtr, sphKernel, ccc, cHash);
//}
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//void SphInterpolationFuncs::sumGradWDotGradW_host_(
//	float& result,
//	unsigned int idP, const float* pxs, const float* pys, const float* pzs,
//	const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash)
//{
//	result = 0;
//	neighborIterCalculation_host_ < OperatorSumGradWDotGradW_ , float > (&result, idP, pxs, pys, pzs, sortedIdMapsPtr, sphKernel, ccc, cHash);
//}
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//void SphInterpolationFuncs::sumValW_host_(
//	float& result,
//	unsigned int idP, const float* pxs, const float* pys, const float* pzs, const float* neighborScalarVals,
//	const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash)
//{
//	result = 0;
//	sphScalarInterpolation_host_ < OperatorSumWTimesNSclar_ , float > (&result, idP, pxs, pys, pzs, neighborScalarVals, sortedIdMapsPtr, sphKernel, ccc, cHash);
//}
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//void SphInterpolationFuncs::sumValGradW_host_(
//	float* result,
//	unsigned int idP, const float* pxs, const float* pys, const float* pzs, const float* neighborScalarVals,
//	const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash)
//{
//	result[0] = 0;
//	result[1] = 0;
//	result[2] = 0;
//	sphScalarInterpolation_host_ < OperatorSumGradWTimesNSclar_ , float > (result, idP, pxs, pys, pzs, neighborScalarVals, sortedIdMapsPtr, sphKernel, ccc, cHash);
//}
//
//
////-------------------------------------------------------------------
////-------------------------------------------------------------------
//void SphInterpolationFuncs::sumValLaplaceW_host_(
//	float& result,
//	unsigned int idP, const float* pxs, const float* pys, const float* pzs, const float* neighborScalarVals,
//	const unsigned int* sortedIdMapsPtr, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash)
//{
//	result = 0;
//	sphScalarInterpolation_host_ < OperatorSumLaplaceWTimesNSclar_ , float > (&result, idP, pxs, pys, pzs, neighborScalarVals, sortedIdMapsPtr, sphKernel, ccc, cHash);
//}
