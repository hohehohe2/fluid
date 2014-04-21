#ifndef hohe_basicGeos_H
#define hohe_basicGeos_H

#include <cudaCommon/Buffer.h>

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Points.
struct Points : public BufferSetSized
{

	//Constructor.
	Points(const std::string& name="points", unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST) : BufferSetSized(name)
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

	//Make sure that every member has the same size.
	void setXs(BufferFloat* xs){removeChild(m_xs); addChild(m_xs = xs);}
	void setYs(BufferFloat* ys){removeChild(m_ys); addChild(m_ys = ys);}
	void setZs(BufferFloat* zs){removeChild(m_zs); addChild(m_zs = zs);}

public:

	BufferFloat* m_xs;
	BufferFloat* m_ys;
	BufferFloat* m_zs;

};


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Linees.
struct Lines : public BufferSetSized
{

	//Constructor.
	Lines(const std::string& name="lines", unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST) : BufferSetSized(name)
	{
		addChild(m_starts = new Points("line starts"));
		addChild(m_ends = new Points("line ends"));
		setSize(size);
		allocate(allocMemoryType);
	}

	//Destructor.
	virtual ~Lines(){}

	virtual unsigned int size() const {return m_starts->size();}

	//Make sure that both members have the same size.
	void setStarts(Points* starts){removeChild(m_starts); addChild(m_starts = starts);}
	void setEnds(Points* ends){removeChild(m_ends); addChild(m_ends = ends);}

public:

	Points* m_starts;
	Points* m_ends;

};


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Triangles.
struct Triangles : public BufferSet
{

	//Constructor.
	/**
	@param name Name of this object.
	@param size Data size of m_v0s, m_v1s, and m_v2s.
	@param allocMemoryType Set BufferSet::HOST, BufferSet::DEVICE, or BufferSet::HOST|BufferSet::DEVICE to allocate memory.
	@param vtxPos Position of the triangle vertices.
	**/
	Triangles(const std::string& name="triangles", unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST, Points* vtxPos=NULL) : BufferSet(name)
	{
		addChild(m_v0s = new BufferUInt("triangle v0s"));
		addChild(m_v1s = new BufferUInt("triangle v1s"));
		addChild(m_v2s = new BufferUInt("triangle v2s"));

		if ( ! vtxPos)
		{
			vtxPos = new Points("triangle pos");
		}
		m_vtxPos = vtxPos;
		addChild(vtxPos);

		m_v0s->setSize(size);
		m_v1s->setSize(size);
		m_v2s->setSize(size);
		allocate(allocMemoryType);
	}

	//Destructor.
	virtual ~Triangles(){}

	//Make sure that m_v0s, m_v1s, and m_v2s have the same size.
	void setV0s(BufferUInt* v0s){removeChild(m_v0s); addChild(m_v0s = v0s);}
	void setV1s(BufferUInt* v1s){removeChild(m_v1s); addChild(m_v1s = v1s);}
	void setV2s(BufferUInt* v2s){removeChild(m_v2s); addChild(m_v2s = v2s);}
	void setVtxPos(Points* vtxPos){removeChild(m_vtxPos); addChild(m_vtxPos = vtxPos);}

public:

	BufferUInt* m_v0s;
	BufferUInt* m_v1s;
	BufferUInt* m_v2s;
	Points* m_vtxPos;
};

}

#endif
