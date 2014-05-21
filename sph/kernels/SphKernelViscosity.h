#ifndef hohe_SphKernelViscosity_H
#define hohe_SphKernelViscosity_H

#define _USE_MATH_DEFINES
#include <math.h>

#include "SphKernel.h"

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///SPH kernel for viscosity (Muller03).
class SphKernelViscosity : public SphKernel
{

public:

	///Laplacian of wViscosity = laplaceWPart() * getConstant().
	inline float getConstant() const {return (float)(45.0 / M_PI / m_r6);}

	///Laplacian of wViscosity = laplaceWPart() * getConstant().
	inline float laplaceWPart(const float dist2) const
	{
		if (dist2 < m_r * m_r * 1.0E-4f || dist2 > m_r * m_r)
		{
			//Same particle or outside the kernel radius.
			return 0.0f;
		}

		return m_r - dist2;
	}
};

}

#endif
