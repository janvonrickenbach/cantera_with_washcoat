/**
 *  @file ImplicitSurfChem.cpp
 * Definitions for the implicit integration of surface site density equations
 *  (see \ref  kineticsmgr and class
 *  \link Cantera::ImplicitSurfChem ImplicitSurfChem\endlink).
 */
// Copyright 2001  California Institute of Technology

/*

*/

#include "cantera/kinetics/ImplicitSurfChem_wc.h"
#include "cantera/kinetics/wcdata.h"
#include "cantera/numerics/Integrator.h"
#include "cantera/kinetics/solveSP.h"
#include "cantera/thermo/SurfPhase.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "cantera/kinetics/SingleWc.h"

namespace Cantera{
// Constructor

ImplicitSurfChem_wc::ImplicitSurfChem_wc(InterfaceKinetics* k
                                        ,Transport* t
                                        ,double h,double h_temp
                                        ,double wc_thickness, double area_to_volume
                                        ,double porosity, double tortuosity
                                        ,double d_p, double lambda_solid
                                        ,double structure_porosity, double lambda_solid_st  
                                        ,double atol, double rtol
                                        ,int nx, bool with_energy
                                        ,int cells_x,double L_r, double vel, double dt, double A_V,bool from_file, int maxsteps
                                        ,double mintemp,double maxtemp, double trate, double rhocp, double rhocp_st, bool inf_ext_mt
                                        ,double bulk_pressure, double heat_source, double cell_ratio):
    m_atol(atol),
    m_rtol(rtol),
    FuncEval(),
    m_integ(0),
    m_maxstep(0.0),
    m_cells_x(cells_x)

{

   //Initial integrator
    m_integ = newIntegrator("CVODE");
    // use backward differencing, with a full Jacobian computed
    // numerically, and use a Newton linear iterator
    m_integ->setMethod(BDF_Method);
 //   m_integ->setProblemType(DENSE+NOJAC);
    m_integ->setProblemType(GMRES);
    m_integ->setIterator(Newton_Iter);
    m_integ->setMaxSteps(maxsteps);
    m_integ->setMaxErrTestFails(100);
    int bandwidth =  k->nTotalSpecies();
    if (with_energy) {
     //  bandwidth +=  bandwidth;
     //  bandwidth = bandwidth * nx;
       bandwidth = -1;
    }
    m_integ->setBandwidth(bandwidth,bandwidth);
//    m_integ->setMaxOrder(1);
    for (int cell_idx=0;cell_idx < m_cells_x;++cell_idx){
       wc_list.push_back(new SingleWc(this,k,t,h,h_temp,wc_thickness,area_to_volume,porosity,tortuosity
               ,d_p,lambda_solid,structure_porosity,lambda_solid_st,nx,with_energy,cell_idx,cells_x,L_r,vel,A_V,from_file,mintemp,maxtemp,trate
               ,rhocp, rhocp_st, inf_ext_mt,bulk_pressure,heat_source,cell_ratio));
    }

}


void ImplicitSurfChem_wc::eval(doublereal time, doublereal* y,
                               doublereal* ydot, doublereal* p)
{
    for (int cell_idx=0;cell_idx < m_cells_x;++cell_idx){
       wc_list[cell_idx]->eval(time,y,
                               ydot,p);

    }
}

/*
 * Destructor. Deletes the integrator.
 */
ImplicitSurfChem_wc::~ImplicitSurfChem_wc()
{
   delete m_integ;
    for (int cell_idx=0;cell_idx < m_cells_x;++cell_idx){
       delete wc_list[cell_idx];
    }
}

void ImplicitSurfChem_wc::getInitialConditions(doublereal t0, size_t lenc,
        doublereal* y)
{
    for (int cell_idx=0;cell_idx < m_cells_x;++cell_idx){
       wc_list[cell_idx]->getInitialConditions(y);

    }
}

/*
 *  Must be called before calling method 'advance'
 */
void ImplicitSurfChem_wc::initialize(doublereal t0)
{
   m_integ->setTolerances(m_rtol, m_atol);
   m_integ->initialize(t0, *this);
}

void ImplicitSurfChem_wc::reinitialize(doublereal t0)
{
   m_integ->reinitialize(t0, *this);
}


// Integrate from t0 to t1. The integrator is reinitialized first.
/*
 *   This routine does a time accurate solve from t = t0 to t = t1.
 *   of the surface problem.
 *
 *  @param t0  Initial Time -> this is an input
 *  @param t1  Final Time -> This is an input
 */
void ImplicitSurfChem_wc::integrate(doublereal& t0, doublereal t1,int maxsteps)
{
    m_integ->setMaxStepSize(t1 - t0);
//    m_integ->setMaxStepSize(100*t1/maxsteps);
//    m_integ->setMaxSteps(maxsteps);
//    m_integ->initialize(t0, *this);
//    m_integ->setMaxOrder(1);
    m_integ->integrate(t1,t0);
    t0 = t1;
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
    m_integ->integrate(t1,t0);
}


/*
 * Called by the integrator to evaluate ydot given y at time 'time'.
 */

size_t ImplicitSurfChem_wc::neq() {return wc_list[0]->neq()*m_cells_x;}

void ImplicitSurfChem_wc::write_wc(int step){
    for (int cell_idx=0;cell_idx < m_cells_x;++cell_idx){
       wc_list[cell_idx]->write_wc(step);

    }

}

void ImplicitSurfChem_wc::set_bulk_temperature(double temperature){
    for (int cell_idx=0;cell_idx < m_cells_x;++cell_idx){
       wc_list[cell_idx]->set_bulk_temperature(temperature);

    }

}


} // namespace Cantera
