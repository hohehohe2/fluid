#ifndef hohe_SphKernelPoly6_H
#define hohe_SphKernelPoly6_H

#define _USE_MATH_DEFINES
#include <math.h>

#include "SphKernel.h"

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///SPH kernel (Muller03).
class SphKernelPoly6 : public SphKernel
{

public:

	///wPoly6 = wPart() * getConstant().
	inline float getConstant() const {return (float)(315.0 / 64.0 / M_PI / m_r9);}

	///wPoly6 = wPart() * getConstant().
	inline float wPart(float dist2) const
	{
		if (dist2 < m_r * m_r)
		{
			float l = m_r * m_r - dist2;
			return (float)(l * l * l);
		}
		else
		{
			return 0.0f;
		}
	}
};

}

#endif
