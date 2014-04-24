#ifndef hohe_FluidSolverSimpleSph_H
#define hohe_FluidSolverSimpleSph_H

#include <hohe2Common/cuda/Buffer.h>
#include "SphKernel.h"
#include "hohe2Common/geo/basicGeos.h"


namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Simple SPH.
/**
Physical quantities are measured in MKS, not normalized here.
**/
class FluidSolverSimpleSph
{

public:

	///Particles for this solver.
	struct Particles : public BufferSetSized
	{

		///Particles position. Always exists.
		PointSet* m_pos;

		///Particles velocity.
		PointSet* m_velocity;

		///Particles acceleration.
		PointSet* m_acceleration;

		///Particles density.
		BufferFloat* m_density;

		///Constructor.
		Particles(const std::string& name="Particles") : BufferSetSized(name), m_pos(NULL), m_velocity(NULL), m_acceleration(NULL), m_density(NULL){}

		///Destructor.
		virtual ~Particles(){}

		virtual unsigned int size() const {return (m_pos)? m_pos->size() : 0;}

		//Make sure every size is the same.
		void setPos(PointSet* pos){removeChild(m_pos); addChild(m_pos = pos);}
		void setVelocity(PointSet* velocity){removeChild(m_velocity); addChild(m_velocity = velocity);}
		void setAcceleration(PointSet* acceleration){removeChild(m_acceleration); addChild(m_acceleration = acceleration);}
		void setDensity(BufferFloat* density){removeChild(m_density); addChild(m_density = density);}

		///Create and setup for simulation. m_pos and m_velocity must be filled with initial values before the first step().
		static Particles* createInstance(unsigned int size=0, MemoryType allocMemoryType=HOST)
		{
			//Separated object creation and filling member variables so that we can reuse fillMembers_() in a subclass.
			Particles* obj = new Particles;
			fillMembers_(obj, size, allocMemoryType);
			return obj;
		};

	protected:

		static void fillMembers_(Particles* obj, unsigned int size, MemoryType allocMemoryType)
		{
			obj->addChild(obj->m_pos = new PointSet("particle pos"));
			obj->addChild(obj->m_velocity = new PointSet("particle velocity"));
			obj->addChild(obj->m_acceleration = new PointSet("particle acceleration"));
			obj->addChild(obj->m_density = new BufferFloat("particle density"));
			obj->setSize(size);
			obj->allocate(allocMemoryType);
			obj->m_pos->memset(0, HOST);
			obj->m_velocity->memset(0, HOST);
		}

	};

	///Constructor.
	/**
	@param particleMass in MKS.
	**/
	FluidSolverSimpleSph(float particleMass=1.0f);

	///Destructor.
	~FluidSolverSimpleSph(){}

	///Step the simulation.
	void step(Particles& particles, float deltaT);

private:

	///Mass per single particle.
	float m_particleMass;

	///SPH kernel calculator.
	SphKernel m_sphKernel;

private:

	float calcMaxVelocity_(Particles& particles);
	void updateNeighbors_(Particles& particles);
	void calcDensity_(Particles& particles);
	void calcAcceleration_(Particles& particles);
	void integrate_(Particles& particles, float deltaT);

	float densityToPressure_(float density);

	///Density->pressure stiffness coefficient.
	/**
	It is not in Constants.h because this is an artificially enough physical value.
	**/
	static const float K;

    ///So called Courant number (pet peeve for Prof. Bridson ;).
    static const float PET_PEEVE_COURANT_NUMBER;

};

}

#endif
