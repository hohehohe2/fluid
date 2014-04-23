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

	///Constructor.
	Points(const std::string& name="points", unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST) : BufferSetSized(name)
	{
		addChild(m_xyzs = new BufferFloat("point xyzs"));
		setSize(size);
		allocate(allocMemoryType);
	}

	///Destructor.
	virtual ~Points(){}

	///Overrided version of sizeSize().
	virtual void setSize(unsigned int size){m_xyzs->setSize(size * 3);}

	///Overrided version of size().
	virtual unsigned int size() const {return m_xyzs->size() / 3;}

	///Set the new buffer for this Points object. No memory allocation / resize will be made. New buffer xyzs cannot be NULL.
	void setXyzs(BufferFloat* xyzs)
	{
		assert(xyzs);
		removeChild(m_xyzs);
		addChild(m_xyzs = xyzs);
	}

	float* xs(MemoryType mType){return m_xyzs->get(mType);} ///Accessor.
	float* ys(MemoryType mType){return m_xyzs->get(mType) + size();} ///Accessor.
	float* zs(MemoryType mType){return m_xyzs->get(mType) + size() * 2;} ///Accessor.
	const float* xs(MemoryType mType) const{return m_xyzs->get(mType);} ///Accessor.
	const float* ys(MemoryType mType) const{return m_xyzs->get(mType) + size();} ///Accessor.
	const float* zs(MemoryType mType) const{return m_xyzs->get(mType) + size() * 2;} ///Accessor.
	float* xs(bool isHost){return m_xyzs->get(boolTo_(isHost));} ///Overloaded version of the accessor.
	float* ys(bool isHost){return m_xyzs->get(boolTo_(isHost)) + size();} ///Overloaded version of the accessor.
	float* zs(bool isHost){return m_xyzs->get(boolTo_(isHost)) + size() * 2;} ///Overloaded version of the accessor.
	const float* xs(bool isHost) const{return m_xyzs->get(boolTo_(isHost));} ///Overloaded version of the accessor.
	const float* ys(bool isHost) const{return m_xyzs->get(boolTo_(isHost)) + size();} ///Overloaded version of the accessor.
	const float* zs(bool isHost) const{return m_xyzs->get(boolTo_(isHost)) + size() * 2;} ///Overloaded version of the accessor.

public:

	///Buffer for the points. Use setXyzs() to assign a different buffer.
	BufferFloat* m_xyzs;

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

	//Make sure that m_v0s, m_v1s, m_v2s and m_vtxPos data are consistent.
	void setVtxPos(Points* vtxPos){removeChild(m_vtxPos); addChild(m_vtxPos = vtxPos);}

public:

	BufferUInt* m_v0s;
	BufferUInt* m_v1s;
	BufferUInt* m_v2s;
	Points* m_vtxPos;
};

}

#endif
