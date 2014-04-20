#ifndef hohe_Geos_H
#define hohe_Geos_H

#include <cudaCommon/Buffer.h>

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Points.
struct Points : public BufferSetSameSizedChildren
{

	typedef std::shared_ptr < Points > SPtr;

	//Constructor.
	Points(const std::string& name="", unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST) : BufferSetSameSizedChildren(name)
	{
		addChild(m_xs = new BufferFloat("point xs"));
		addChild(m_ys = new BufferFloat("point ys"));
		addChild(m_zs = new BufferFloat("point zs"));
		setSize(size);
		allocate(allocMemoryType);
	}

	//Destructor.
	virtual ~Points(){}

	virtual unsigned int size() const {return m_xs->size();}

public:

	BufferFloat* m_xs;
	BufferFloat* m_ys;
	BufferFloat* m_zs;

};


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Linees.
struct Lines : public BufferSetSameSizedChildren
{

	typedef std::shared_ptr < Lines > SPtr;

	//Constructor.
	Lines(const std::string& name="", unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST) : BufferSetSameSizedChildren(name)
	{
		addChild(m_starts = new Points("line starts"));
		addChild(m_ends = new Points("line ends"));
		setSize(size);
		allocate(allocMemoryType);
	}

	//Destructor.
	virtual ~Lines(){}

	virtual unsigned int size() const {return m_starts->size();}

	void init(unsigned int size, bool starts=true, bool ends=true);

public:

	Points* m_starts;
	Points* m_ends;

};


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Triangles.
struct Triangles : public BufferSetSameSizedChildren
{

	typedef std::shared_ptr < Triangles > SPtr;

	//Constructor.
	Triangles(const std::string& name="", unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST) : BufferSetSameSizedChildren(name)
	{
		addChild(m_v0s = new BufferUInt("triangle v0s"));
		addChild(m_v1s = new BufferUInt("triangle v1s"));
		addChild(m_v2s = new BufferUInt("triangle v2s"));
		addChild(m_pos = new Points("triangle pos"));
		setSize(size);
		allocate(allocMemoryType);
	}

	//Destructor.
	virtual ~Triangles(){}

	virtual unsigned int size() const{return (m_pos)? m_pos->size() : 0;}

	void init(unsigned int size, bool pos, bool v0s=true, bool v1s=true, bool v2s=true);

public:

	BufferUInt* m_v0s;
	BufferUInt* m_v1s;
	BufferUInt* m_v2s;
	Points* m_pos;
};

}

#endif
