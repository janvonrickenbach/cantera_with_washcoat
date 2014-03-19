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
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "solveSP.h"

namespace Cantera
{

class solveSP;


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
    ImplicitSurfChem_wc(std::vector<InterfaceKinetics*> k);

    /**
     * Destructor. Deletes the integrator.
     */
    virtual ~ImplicitSurfChem_wc();

    /**
     * Overloads the virtual function
     * declared in FuncEval.
     */
    virtual void initialize(doublereal t0 = 0.0);

    virtual void set_wc_coefficient(doublereal h);
    virtual void set_transport(Transport* t);


    //! Integrate from t0 to t1. The integrator is reinitialized first.
    /*!
     *   This routine does a time accurate solve from t = t0 to t = t1.
     *   of the surface problem.
     *
     *  @param t0  Initial Time -> this is an input
     *  @param t1  Final Time -> This is an input
     */
    void integrate(doublereal t0, doublereal t1);


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


    // overloaded methods of class FuncEval

    //! Return the number of equations
    virtual size_t neq() {
        return m_nv;
    }

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

    //! Set the initial conditions for the solution vector
    /*!
     *  @param t0  Initial time
     *  @param leny  Length of the solution vector
     *  @param y   Value of the solution vector to be used.
     *            On output, this contains the initial value
     *           of the solution.
     */


protected:


    //! vector of pointers to surface phases.
    SurfPhase* m_surf;

    //! Vector of pointers to bulk phases
    ThermoPhase* m_bulk;

    std::vector<double*> m_bulk_massfraction;
    std::vector<double*> m_diff_coeffs;

    Transport* m_transport;


    //! vector of pointers to InterfaceKinetics objects
    InterfaceKinetics* m_kin;

    //! Vector of number of species in each Surface Phase
    std::vector<size_t>  m_nsp_surf;
    std::vector<size_t>  m_nsp_tot;


    std::vector<vector_int> pLocVec;
    //! Pointer to the cvode integrator
    Integrator* m_integ;
    doublereal m_atol, m_rtol;   // tolerances
    doublereal m_maxstep;        // max step size

    doublereal m_wc_coefficient;


    friend class solveSS;

    double getMassFrac(int loc_idx,mass_idx);
    double getTemp(int loc_idx);
    double getP();

private:

    vector_fp m_state;

    //! Controls the amount of printing from this routine
    //! and underlying routines.
    int m_ioFlag;
};
}

#endif

