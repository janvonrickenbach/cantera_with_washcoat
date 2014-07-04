/**
 *  @file ImplicitSurfChem_wc.h
 * Declarations for the implicit integration of surface site density equations
 *  (see \ref  kineticsmgr and class
 *  \link Cantera::ImplicitSurfChem ImplicitSurfChem\endlink).
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_IMPSURFCHEM_wc_H
#define CT_IMPSURFCHEM_wc_H

#include "cantera/numerics/FuncEval.h"
#include "cantera/numerics/Integrator.h"
//#include "cantera/kinetics/InterfaceKinetics.h"
//#include "cantera/thermo/SurfPhase.h"
//#include "solveSP.h"
#include <vector>

namespace Cantera
{

class solveSP;
class InterfaceKinetics;
class SurfPhase;
class Transport;
class ThermoPhase;
class wcdata;
class SingleWc;


//! Advances the surface coverages of the associated set of SurfacePhase
//! objects in time
/*!
 *  This function advances a set of SurfacePhase objects, each
 *  associated with one InterfaceKinetics object, in time.
 *  The following equation is used for each surface phase, <I>i</I>.
 *
 *   \f[
 *        \dot \theta_k = \dot s_k (\sigma_k / s_0)
 *   \f]
 *
 *  In this equation,
 *    \f$ \theta_k \f$ is the site coverage for the kth species.
 *    \f$ \dot s_k \f$ is the source term for the kth species
 *   \f$ \sigma_k \f$ is the number of surface sites covered by
 *  each species k.
 *   \f$ s_0 \f$ is the total site density of the interfacial phase.
 *
 *  Additionally, the 0'th equation in the set is discarded. Instead the
 *  alternate equation is solved for
 *
 *   \f[
 *        \sum_{k=0}^{N-1}  \dot \theta_k = 0
 *   \f]
 *
 *  This last equation serves to ensure that sum of the \f$ \theta_k \f$
 *  values stays constant.
 *
 *  The object uses the CVODE software to advance the surface equations.
 *
 *  The solution vector used by this object is as follows.
 *   For each surface phase with \f$ N_s \f$ surface sites,
 *   it consists of the surface coverages
 *       \f$ \theta_k \f$ for \f$ k = 0, N_s - 1 \f$
 *
 * @ingroup  kineticsmgr
 *
 */
class ImplicitSurfChem_wc : public FuncEval
{

public:


    //! Constructor for multiple surfaces.
    /*!
     * @param k  Vector of pointers to InterfaceKinetics objects
     *           Each object consists of a surface or an edge containing
     *           internal degrees of freedom representing the concentration
     *           of surface adsorbates.
     */
    ImplicitSurfChem_wc(InterfaceKinetics* k
    		           ,Transport* t
    		           ,double h,double h_temp
    		           ,double wc_thickness, double area_to_volume
    		           ,double porosity, double tortuosity
    		           ,double d_p, double lambda_solid
    		           ,double atol, double rtol
		               ,int nx,bool with_energy, int cells_x
		              ,double L_r, double vel, double dt, double A_V,bool from_file, int maxsteps
		              ,double mintemp, double maxtemp, double trate);

    /**
     * Destructor. Deletes the integrator.
     */
    virtual ~ImplicitSurfChem_wc();

    /**
     * Overloads the virtual function
     * declared in FuncEval.
     */
    virtual void initialize(doublereal t0 = 0.0);

    /**
     * Overloads the virtual function
     * declared in FuncEval.
     */
    virtual void reinitialize(doublereal t0=0.0);


    //! Integrate from t0 to t1. The integrator is reinitialized first.
    /*!
     *   This routine does a time accurate solve from t = t0 to t = t1.
     *   of the surface problem.
     *
     *  @param t0  Initial Time -> this is an input
     *  @param t1  Final Time -> This is an input
     */
    void integrate(doublereal& t0, doublereal t1, int maxiter);


    //! Integrate from t0 to t1 without reinitializing the integrator.
    /*!
     *  Use when the coverages have not changed from
     *  their values on return from the last call to integrate or
     *  integrate0.
     *
     *  @param t0  Initial Time -> this is an input
     *  @param t1  Final Time -> This is an input
     */
    void integrate0(doublereal t0, doublereal t1);



    //! Evaluate the value of ydot[k] at the current conditions
    /*!
     *  @param t   Time (seconds)
     *  @param y   Vector containing the current solution vector
     *  @param ydot   Output vector containing the value of the
     *                derivative of the surface coverages.
     *  @param p   Unused parameter pass-through parameter vector
     */
    virtual void eval(doublereal t, doublereal* y, doublereal* ydot,
                      doublereal* p);

    void set_bulk_temperature(double temperature);

    // This is called by the CVODE integrator in the beginning
    // to initialize the state vector
    void getInitialConditions(doublereal t0, size_t lenc,doublereal* c);

    // Returns the length of the state vector. Called by CVODE
    size_t neq();
    void write_wc(int step);
    Integrator* m_integ;



protected:



    //! Pointer to the cvode integrator
    doublereal m_atol, m_rtol;   // tolerances
    doublereal m_maxstep;        // max step size

    int m_cells_x;
    std::vector<SingleWc*> wc_list;



};

}

#endif

