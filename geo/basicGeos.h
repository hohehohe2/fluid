#ifndef hohe_basicGeos_H
#define hohe_basicGeos_H

#include <cudaCommon/BufferSet.h>

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Points.
struct Points : public BufferSetSameSizedMember
{

	//Constructor.
	Points(unsigned int size=0)
		: m_xs("xs", size), m_ys("ys", size), m_zs("zs", size)
	{
		addMember_(m_xs);
		addMember_(m_ys);
		addMember_(m_zs);
	}

	//Destructor.
	virtual ~Points(){}

	virtual unsigned int size() const {return m_xs.size();}

public:

	Buffer < float > m_xs;
	Buffer < float > m_ys;
	Buffer < float > m_zs;

};


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Linees.
struct Lines : public BufferSetSameSizedMember
{

	//Constructor.
	Lines(unsigned int size=0) : m_start(size), m_end(size)
	{
		addMember_(m_start);
		addMember_(m_end);
	}

	//Destructor.
	virtual ~Lines(){}

	virtual unsigned int size() const {return m_start.size();}

public:

	Points m_start;
	Points m_end;

};


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Triangles.
struct Triangles : public BufferSetSameSizedMember
{

	//Constructor.
	Triangles(unsigned int size=0) : m_vtx0(size), m_vtx1(size), m_vtx2(size)
	{
		addMember_(m_vtx0);
		addMember_(m_vtx1);
		addMember_(m_vtx2);
	}

	//Destructor.
	virtual ~Triangles(){}

	virtual unsigned int size() const {return m_vtx0.size();}

public:

	Points m_vtx0;
	Points m_vtx1;
	Points m_vtx2;

};


}

#endif
