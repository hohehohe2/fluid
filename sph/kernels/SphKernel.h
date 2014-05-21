#ifndef hohe_SphKernel_H
#define hohe_SphKernel_H


namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Base class of SPH kernels.
class SphKernel
{

public:

	///Set the kernel radius.
	void setKernelRadius(float radius)
	{
		m_r = radius;
		m_r6 = m_r * m_r * m_r * m_r * m_r * m_r;
		m_r9 = m_r * m_r * m_r * m_r * m_r * m_r * m_r * m_r * m_r;
	}

	///Get the kernel radius.
	float kernelRadius() const{return m_r;}

protected:

	///Kernel radius.
	float m_r;

	///Kernel radius ^ 6.
	float m_r6;

	///Kernel radius ^ 9.
	float m_r9;
};

}

#endif
