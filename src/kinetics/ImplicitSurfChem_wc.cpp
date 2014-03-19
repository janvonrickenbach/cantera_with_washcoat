/**
 *  @file ImplicitSurfChem.cpp
 * Definitions for the implicit integration of surface site density equations
 *  (see \ref  kineticsmgr and class
 *  \link Cantera::ImplicitSurfChem ImplicitSurfChem\endlink).
 */
// Copyright 2001  California Institute of Technology

#include "cantera/kinetics/ImplicitSurfChem_wc.h"
#include "cantera/numerics/Integrator.h"
#include "cantera/kinetics/solveSP.h"
#include "cantera/thermo/SurfPhase.h"
#include <iostream>
#include <fstream>


namespace Cantera{
// Constructor
ImplicitSurfChem_wc::ImplicitSurfChem_wc(InterfaceKinetics* k) :
    FuncEval(),
    m_integ(0),
    m_atol(1.e-25),
    m_rtol(1.e-7),
    m_maxstep(0.0)

{

	m_wc_thickness = 1E-5;

	createGrid();
	printGrid();

    m_integ = newIntegrator("CVODE");


    // use backward differencing, with a full Jacobian computed
    // numerically, and use a Newton linear iterator

    m_integ->setMethod(BDF_Method);
    m_integ->setProblemType(DENSE + NOJAC);
    m_integ->setIterator(Newton_Iter);
    m_state.resize(0);

    if (k->nPhases() > 2){
    	std::cout << "Error more than two phases" << std::endl;
    }

    surface_phase = (SurfPhase* )&k->thermo(k->surfacePhaseIndex());
    for (int phase_idx=0;phase_idx< k->nPhases();++phase_idx){
      if (phase_idx /= k->surfacePhaseIndex()){
        gas_phase    = &k->thermo(phase_idx);
      }
    }

    m_vol_sp  = gas_phase->nSpecies();
    m_surf_sp = k->nTotalSpecies() - m_vol_sp;


    for (int comp_idx=0;comp_idx < m_vol_sp;++comp_idx){
    	vars_enum temp = static_cast<vars_enum>(en_end+comp_idx);
    	en_vol_comp.push_back(temp);
    }

    for (int comp_idx=0;comp_idx < m_surf_sp;++comp_idx){
    	vars_enum temp = static_cast<vars_enum>(en_end+m_vol_sp+comp_idx);
    	en_surf_comp.push_back(temp);
    }

    m_bulk_massfraction.resize(m_vol_sp);
    m_fluxes.resize(m_nx+1);

    m_rho.resize(m_nx);
    m_diffusion_coeffs.resize(m_vol_sp);

    m_x.resize(m_nx);
    m_dx.resize(m_nx);
    m_x_co.resize(m_nx+1);

    comp_vector::iterator comp_iter;

    for (comp_iter =  m_diffusion_coeffs.begin();
    	 comp_iter != m_diffusion_coeffs.end();
    	 ++comp_iter){
    	comp_iter->resize(m_nx);
    }

    m_nvars = 1+m_surf_sp + m_vol_sp;
}


/*
 * Destructor. Deletes the integrator.
 */
ImplicitSurfChem_wc::~ImplicitSurfChem_wc()
{
    delete m_integ;
}


/*
 *  Must be called before calling method 'advance'
 */
void ImplicitSurfChem_wc::initialize(doublereal t0)
{
    m_integ->setTolerances(m_rtol, m_atol);
    m_integ->initialize(t0, *this);
}

/*
 *  Set masstransfer coefficient
 */
void ImplicitSurfChem_wc::set_wc_coefficient(doublereal h)
{
    m_wc_coefficient = h;
}

void ImplicitSurfChem_wc::set_transport(Transport* t)
{
    m_transport = t;
}

// Integrate from t0 to t1. The integrator is reinitialized first.
/*
 *   This routine does a time accurate solve from t = t0 to t = t1.
 *   of the surface problem.
 *
 *  @param t0  Initial Time -> this is an input
 *  @param t1  Final Time -> This is an input
 */
void ImplicitSurfChem_wc::integrate(doublereal t0, doublereal t1)
{
    m_integ->initialize(t0, *this);
    m_integ->setMaxStepSize(t1 - t0);
    m_integ->integrate(t1);
}

// Integrate from t0 to t1 without reinitializing the integrator.
/*
 *  Use when the coverages have not changed from
 *  their values on return from the last call to integrate or
 *  integrate0.
 *
 *  @param t0  Initial Time -> this is an input
 *  @param t1  Final Time -> This is an input
 */
void ImplicitSurfChem_wc::integrate0(doublereal t0, doublereal t1)
{
    m_integ->integrate(t1);
}


/*
 * Called by the integrator to evaluate ydot given y at time 'time'.
 */
void ImplicitSurfChem_wc::eval(doublereal time, doublereal* y,
                            doublereal* ydot, doublereal* p)
{
}



doublereal ImplicitSurfChem_wc::getStateVar(vars_enum var,int loc_idx){
	return m_state[loc_idx*m_nvars+var];

}
void       ImplicitSurfChem_wc::setStateVar(double value,vars_enum var,int loc_idx){
	m_state[loc_idx*m_nvars+var] = value;
}

void ImplicitSurfChem_wc::createGrid(){
	//compute dx in every cell
	grid_vector::iterator vec_it;
	grid_vector::iterator vec_it_co;
	grid_vector::iterator vec_it_dx;

	for(vec_it  = m_dx.begin();
		vec_it != m_dx.end();
		++vec_it){
	   *vec_it = m_wc_thickness / m_nx;
	}

	//compute interpolation coeffcient in every cell
	for(vec_it  = m_fx.begin();
		vec_it != m_fx.end();
		++vec_it){
	   *vec_it = 0.5;
	}

	// Get the corner values
	std::partial_sum(m_dx.begin(),m_dx.end(),m_x_co.begin());

	// Get the cell_center values
	for(vec_it    = m_x.begin()
	   ,vec_it_co = m_x_co.begin()
	   ,vec_it_dx = m_dx.begin();
		vec_it != m_x.end();
		++vec_it,++vec_it_co,++vec_it_dx){
	   *vec_it = *vec_it_co + 0.5 * (*vec_it_dx);
	}
}

void ImplicitSurfChem_wc::printGrid(){

	std::ofstream myfile;
	myfile.open ("grid.dat");

	grid_vector::iterator vec_it;
	grid_vector::iterator vec_it_co;
	grid_vector::iterator vec_it_dx;

	for(vec_it    = m_x.begin()
	   ,vec_it_co = m_x_co.begin()
	   ,vec_it_dx = m_dx.begin();
		vec_it != m_x.end();
		++vec_it,++vec_it_co,++vec_it_dx){
		myfile << *vec_it    << " "
			   << *vec_it_co << " "
			   << *vec_it_dx << std::endl;
	}

	myfile.close();

}


void ImplicitSurfChem_wc::getInitialConditions(doublereal t0, size_t lenc,
        doublereal* c)
{
}

}
