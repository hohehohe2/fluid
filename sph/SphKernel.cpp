#include "SphKernel.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace hohehohe2;


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void SphKernel::setKernelRadius(float radius)
{
	m_r = radius;
	m_r6 = m_r * m_r * m_r * m_r * m_r * m_r;
	m_r9 = m_r * m_r * m_r * m_r * m_r * m_r * m_r * m_r * m_r;
}

////---------------------------------------------------------------------------
////---------------------------------------------------------------------------
//float SphKernel::w(float dist2) const
//{
//	//wPoly6.
//
//	if (dist2 < m_r * m_r)
//	{
//		float l = m_r * m_r - dist2;
//		return (float)(315.0 / 64.0 / M_PI / m_r9 * l * l * l);
//	}
//	else
//	{
//		return 0.0f;
//	}
//}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
float SphKernel::w(float dist2) const
{
	//Cubic spline.

	if (dist2 < m_r * m_r)
	{
		const float q = sqrt(dist2) / m_r;
		if (q < 0.5f)
		{
			return 8.0f / ((float)M_PI * m_r * m_r * m_r) * (1.0f - 6.0f * q * q + 6 * q * q * q);
		}
		else
		{
			return 8.0f / ((float)M_PI * m_r * m_r * m_r) * 2 * ( 1.0f - q ) * ( 1.0f - q ) * ( 1.0f - q );
		}
	}
	else
	{
		return 0.0f;
	}
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void SphKernel::gradW(float result[3], float px, float py, float pz, float nx, float ny, float nz) const
{
	const float dx = px - nx;
	const float dy = py - ny;
	const float dz = pz - nz;

	const float dist2 = dx * dx + dy * dy + dz * dz;

	if (dist2 < m_r * m_r * 1.0E-4f || dist2 > m_r * m_r)
	{
		//Same particle or outside the kernel radius.
		result[0] = 0.0f;
		result[1] = 0.0f;
		result[2] = 0.0f;
		return;
	}

	const float dist  = sqrt(dist2);

	const float l = m_r - dist;
	const float gradWSpiky = (float)(-45.0 / M_PI / m_r6 * l * l);
	result[0] = gradWSpiky * dx / dist;
	result[1] = gradWSpiky * dy / dist;
	result[2] = gradWSpiky * dz / dist;
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
float SphKernel::laplaceW(const float dist2) const
{
	//Laplacian of wViscosity.

	if (dist2 < m_r * m_r * 1.0E-4f || dist2 > m_r * m_r)
	{
		//Same particle or outside the kernel radius.
		return 0.0f;
	}

	const float l = m_r - dist2;
	return (float)(45.0 / M_PI / m_r6 * l);
}
