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
		               ,int nx,bool with_energy);

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

    // This is called by the CVODE integrator in the beginning
    // to initialize the state vector
    void getInitialConditions(doublereal t0, size_t lenc,doublereal* c);

    // Returns the length of the state vector. Called by CVODE
    size_t neq() {return m_nvars*m_nx;}


    // Writes all the variables to file
    void printGrid();

    // Initialize the bulk values from the state of
    // m_gas_phase
    void set_bulk_from_state();

	// Reset the state of m_gas_phase to
	// the bulk state
	void set_state_from_bulk();

	//Getters
	int get_nx() const {
		return m_nx;
	}

	int get_vol_sp() const {
		return m_vol_sp;
	}

	int get_surf_sp() const {
		return m_surf_sp;
	}


	// Get the state from a wcdata object
	void get_state(wcdata& data) const;
	void get_state_object(wcdata& data);
	// Set the state from to wcdata object
	void set_state(double* state, const wcdata& data);

	void set_wcdata(wcdata* wcdata);

   void get_fluxes(double* y);


protected:



    // Type to store varialbes that exist in every cell
	// for a number of components
    typedef std::vector<std::vector<double> > comp_vector;

    // Type for scalar variables on the grid
    typedef std::vector<double> grid_vector;


    // Bulk Values
    std::vector<doublereal> m_bulk_massfraction;
    std::vector<doublereal> m_bulk_diff_coeffs;
    doublereal m_bulk_diff_temp;
    doublereal m_bulk_density;
    doublereal m_bulk_temperature;
    doublereal m_bulk_pressure;

    // Fluxes for the species
    comp_vector m_fluxes;

    // Fluxes for the energy equation
    std::vector<doublereal> m_fluxes_temp;

    // Thermal conductivity for temperature
    std::vector<doublereal> m_diff_temp;

    // Temporary arrays for species
    std::vector<doublereal> m_temp_vol_massfraction;
    std::vector<doublereal> m_temp_diff_coeffs;
    std::vector<doublereal> m_temp_surf_massfraction;

    // Temporary arrays for species production rates
    std::vector<doublereal> m_temp_production_rates;
    // Temporary reaction rates
    std::vector<doublereal> m_temp_rates_of_progress;
    // Temporary arrays for reaction enthalpy
    std::vector<doublereal> m_temp_delta_enthalpy;

    // Number of cell centers
    int m_nx;

    // Number of cell corners
    int m_nco;

    // Number of cell centers +2
    // Includes values at the boundary
    int m_nx_var;

    // Corner coordinates
    grid_vector m_x_co;

    // Cell center coordinates
    // includes the coordinates of the two boundaries
    grid_vector m_x;

    // Dx values for all the cells
    grid_vector m_dx;

    // Interpolation coefficients
    grid_vector m_fx;

    // Washcoat thickness
    doublereal m_wc_thickness;

    // Area to volume ratio
    doublereal m_area_to_volume;

    // Washcoat pore diameter
    doublereal m_dp;

    // Washcoat porostiy
    doublereal  m_porosity;

    // Washcoat tortuosity
    doublereal m_tortuosity;

    // Washcoat solid conductivity
    doublereal m_lambda_solid;

    // Bulk mass-transfer coefficient
    doublereal m_wc_coefficient;

    // Bulk heat-transfer coefficient
    doublereal m_wc_coefficient_temp;

    // Density at each grid point including boundaries
    grid_vector m_rho;

    // Diffusion Coefficients at each grid point including boundaries
    comp_vector m_diff_coeffs;

    // Pointer to transport manager
    Transport* m_transport;

    // Pointer to gas phase
    ThermoPhase* m_gas_phase;

    // Pointer to surface phase
    SurfPhase* m_surface_phase;

    // Pointer to interface kinetics object
    InterfaceKinetics* m_kin;

    // Pointer wcdata object
    wcdata* m_wc_data;

    //! Pointer to the cvode integrator
    Integrator* m_integ;
    doublereal m_atol, m_rtol;   // tolerances
    doublereal m_maxstep;        // max step size


    std::vector<doublereal> m_fluxes_return;

    // Number of volume species
    int m_vol_sp;

    // Number of surface species
    int m_surf_sp;

    bool m_with_energy;

    // Number of variables for the
    // ODE solver (gridpoints *variables)
    int m_nvars;


    // Type for the variable enumerations
    // en_start makes sure they start at zero
    enum vars_enum{en_start=0};

    std::vector<vars_enum> en_vol_comp;
    std::vector<vars_enum> en_surf_comp;
    vars_enum en_temperature;

    // Takes a state vector gets the value of var at grid point
    // loc_idx
    // loc_idx is the index in the state vector, which means boundary
    // points are not included
	doublereal getStateVar(double* state,vars_enum var,int loc_idx) const;

    // Takes a state vector sets value of var at grid point
    // loc_idx
    // loc_idx is the index in the state vector, which means boundary
	// points are not included
	void setStateVar(double* state,double value,vars_enum var,int loc_idx);

	// Creates the grid. This routine takes nx and wc_thickness
	// and buids a grid. It set the corner coordinates and the
	// interpolation cofficinets for the faces
	void createGrid();

    // Updated the fluxes for all the species and the energy equation
	void update_fluxes(double* state);

    // Uses the interpolation coefficients to interpolate the
	// input values th a cell phase
	double interpolate_values(double fxp,double val,double valw);

	// Writes a single variable to a file myfile
	// The vector should have length m_nx_var
	void write_var(std::ofstream& myfile,const std::vector<double>& vec);

	// Updates the material properties
	// Diffusion coefficient, conductivity, density
	// Currently sets it to the bulk state
   void update_material_properties(double* y);


};

}

#endif

