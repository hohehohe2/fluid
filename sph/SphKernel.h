#ifndef hohe_SphKernel_H
#define hohe_SphKernel_H


namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///SPH kernel.
/**
Technically we can fully define the parameters and can create a base interface class, but these methods
tend to be placed inside the inner most loops, so it will be good to leave the method non-virtual.
**/
class SphKernel
{

public:

	///Set the kernel radius.
	void setKernelRadius(float radius);

	///Calculate w.
	float w(float dist2) const;

	///Calculate grad of w.
	void gradW(float result[3], float px, float py, float pz, float nx, float ny, float nz) const;

	///Calculate laplacian of w.
	float laplaceW(const float dist2) const;

	///Get the kernel radius.
	float kernelRadius() const{return m_r;}

private:

	///Kernel radius.
	float m_r;

	///Kernel radius ^ 6.
	float m_r6;

	///Kernel radius ^ 9.
	float m_r9;
};

}

#endif
