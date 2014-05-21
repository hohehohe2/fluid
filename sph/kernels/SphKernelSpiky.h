#ifndef hohe_SphKernelSpiky_H
#define hohe_SphKernelSpiky_H

#define _USE_MATH_DEFINES
#include <math.h>

#include <hohe2Common/geo/basicGeos.h>
#include "SphKernel.h"

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///SPH kernel for pressure (Muller03).
class SphKernelSpiky : public SphKernel
{

public:

	///Gradient of wSpiky = gradWPart() * getConstant().
	inline float getConstant() const {return (float)(-45.0 / M_PI / m_r6);}

	///Gradient of wSpiky = gradWPart() * getConstant().
	inline void gradWPart(Point& result, const Point& pos, const Point& npos) const
	{
		result = pos - npos;
		float dist2 = result.squaredNorm();

		if (dist2 < m_r * m_r * 1.0E-4f || dist2 > m_r * m_r)
		{
			//Same particle or outside the kernel radius.
			result << 0.0f, 0.0f, 0.0f;
		}
		else
		{
			const float dist  = sqrt(dist2);
			const float l = m_r - dist;
			result *= l * l / dist;
		}
	}

};

}

#endif
