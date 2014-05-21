#ifndef hohe_SphKernelCubicSpline_H
#define hohe_SphKernelCubicSpline_H

#define _USE_MATH_DEFINES
#include <math.h>

#include "SphKernel.h"

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Traditional SPH kernel.
class SphKernelCubicSpline : public SphKernel
{

public:

	///wCubicSpline = wPart() * getConstant().
	inline float getConstant() const {return 8.0f / (float)(M_PI * m_r * m_r * m_r);}

	///wCubicSpline = wPart() * getConstant().
	inline float wPart(float dist2) const
	{
		if (dist2 < m_r * m_r)
		{
			const float q = sqrt(dist2) / m_r;
			if (q < 0.5f)
			{
				return 1.0f - 6.0f * q * q + 6 * q * q * q;
			}
			else
			{
				float oneMQ = 1.0f - q;
				return 2.0f * oneMQ * oneMQ * oneMQ;
			}
		}
		else
		{
			return 0.0f;
		}
	}
};

}

#endif
