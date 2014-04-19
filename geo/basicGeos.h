#ifndef hohe_basicGeos_H
#define hohe_basicGeos_H

#include <cudaCommon/Buffer.h>

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Points.
struct Points : public BufferSetSameSized
{

	typedef std::shared_ptr < Points > SPtr;

	//Constructor.
	Points(const std::string& name="") : BufferSetSameSized(name){}

	//Destructor.
	virtual ~Points(){}

	virtual unsigned int size() const {return (m_xs)? m_xs->size() : 0;}

	void init(unsigned int size, bool xs=true, bool ys=true, bool zs=true);

public:

	BufferFloat::SPtr m_xs;
	BufferFloat::SPtr m_ys;
	BufferFloat::SPtr m_zs;

};


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Linees.
struct Lines : public BufferSetSameSized
{

	typedef std::shared_ptr < Lines > SPtr;

	//Constructor.
	Lines(const std::string& name="") : BufferSetSameSized(name){}

	//Destructor.
	virtual ~Lines(){}

	virtual unsigned int size() const {return (m_starts)? m_starts->size() : 0;}

	void init(unsigned int size, bool starts=true, bool ends=true);

public:

	Points::SPtr m_starts;
	Points::SPtr m_ends;

};


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Triangles.
struct Triangles : public BufferSetSameSized
{

	typedef std::shared_ptr < Triangles > SPtr;

	//Constructor.
	Triangles(const std::string& name="") : BufferSetSameSized(name){}

	//Destructor.
	virtual ~Triangles(){}

	virtual unsigned int size() const{return (m_pos)? m_pos->size() : 0;}

	void init(unsigned int size, bool pos, bool v0s=true, bool v1s=true, bool v2s=true);

public:

	BufferUInt::SPtr m_v0s;
	BufferUInt::SPtr m_v1s;
	BufferUInt::SPtr m_v2s;
	Points::SPtr m_pos;
};


}

#endif
