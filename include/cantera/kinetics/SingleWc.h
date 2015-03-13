/*
 * SingleWc.cpp
 *
 *  Created on: Jun 26, 2014
 *      Author: vonrickenbach
 */

#include <vector>
#include"cantera/base/config.h"
#ifndef SINGLEWC_H_
#define SINGLEWC_H_


namespace Cantera {

class InterfaceKinetics;
class Transport;
class wcdata;
class SurfPhase;
class ThermoPhase;
class ImplicitSurfChem_wc;

class SingleWc{


public:
 SingleWc(ImplicitSurfChem_wc* surf_chem,InterfaceKinetics* k, Transport* t
           ,double h,double h_temp
           ,double wc_thickness, double area_to_volume
           ,double porosity, double tortuosity
           ,double d_p, double lambda_solid
           ,double structure_porosity, double lambda_solid_st   
           ,int nx, bool with_energy, int x_idx,int nxcells, double L_r, double vel, double A_V, bool from_file
           ,double mintemp, double maxtemp, double trate,double rhocp, double rhocp_st, bool inf_ext_mt
           ,double bulk_pressure, double heat_source,double cell_ratio);


~SingleWc();

// Get the state of the CVode state vector and
// copy it into the data object. Uses
// the "solution" method of CVodeInt
void get_state(wcdata& data);

// Get the current state of the gas_phase
// and the surface objects and copy it into the
// data object.
void get_state_object(wcdata& data);

// Takes the state in data object and copies it into
// the CVode state vector
void set_state(double* state, const wcdata& data);


 // Set the masstransfer coefficient according to the
 // paper of Mladenov. This is called if the
 // masstransfer coefficient is smaller than 0
 void set_mt_coefficient_mladenov();

 void set_mt_coefficient_kelvin(double coeff);

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

// Set the bulk temperature
void set_bulk_temperature(double temperature);

// Evaluate the rhs of the ODE. This is called by
// CVode
void eval(doublereal time, doublereal* y,
                    doublereal* ydot, doublereal* p);


// Set the initial conditions. This is called by CVode
void getInitialConditions(doublereal* y);

// Calls the write_wc routine of the wash-coat data class.
void write_wc(int step);

// Number of unknowns in the ODE. Called by CVode
std::size_t neq();

    // Type for the variable enumerations
    // en_start makes sure they start at zero
    enum vars_enum{en_start=0};

inline
void   setStateVar(double* state,double value,vars_enum var,int loc_idx,int m_x_idx){
   state[m_x_idx*(m_nvars*(m_nx+1))+ loc_idx*m_nvars+var] = value;
}

inline
doublereal getStateVar(double* state,vars_enum var,int loc_idx, int m_x_idx) const {
   return state[m_x_idx*(m_nvars*(m_nx+1))+loc_idx*m_nvars+var];
}

void setCanterafromState(double* state, int g_idx);

protected:
    // Type to store varialbes that exist in every cell
    // for a number of components
    typedef std::vector<std::vector<double> > comp_vector;

    // Type for scalar variables on the grid
    typedef std::vector<double> grid_vector;

    // Pointer to interface kinetics object
    InterfaceKinetics* m_kin;

    // Pointer to transport manager
    Transport* m_transport;

    // Bulk heat-transfer coefficient
    doublereal m_wc_coefficient_temp;

    // Bulk mass-transfer coefficient
    std::vector<doublereal> m_wc_coefficient;
    doublereal m_wc_coefficient_in;

    // Washcoat thickness
    doublereal m_wc_thickness;

    // Number of cell centers
    int m_nx;

    // Area to volume ratio
    doublereal m_area_to_volume;

    // Decide if energy equation is solved
    bool m_with_energy;

    // Washcoat porosity
    doublereal  m_porosity;

    // Washcoat solid conductivity
    doublereal m_lambda_solid;

    // Washcoat tortuosity
    doublereal m_tortuosity;

    // Washcoat pore diameter
    doublereal m_dp;
 
    // Structure porosity
    doublereal  m_structure_porosity;

    // strut solid conductivity 
    doublereal m_lambda_solid_st; 

    int m_x_idx;

    ImplicitSurfChem_wc* m_surf_chem;

    double m_L_r;
    double m_vel;

    int m_nxcells;

    double m_A_V;

    bool m_from_file;

    double m_mintemp;
    double m_maxtemp;
    double m_trate;

    double m_rhocp;
    double m_rhocp_st;

    bool m_inf_ext_mt;

    doublereal m_bulk_pressure;
    double m_heat_source;

    double m_cell_ratio;

    // Pointer to gas phase
    ThermoPhase* m_gas_phase;

    // Pointer to surface phase
    SurfPhase* m_surface_phase;

    // Bulk Values
    std::vector<doublereal> m_bulk_massfraction;
    std::vector<doublereal> m_bulk_diff_coeffs;
    doublereal m_bulk_diff_temp;
    doublereal m_bulk_density;
    doublereal m_bulk_temperature;

    // Fluxes for the species
    comp_vector m_fluxes;
    std::vector<double> m_inflow_comp;
    double m_inflow_temp;

    // Fluxes for the energy equation
    std::vector<doublereal> m_fluxes_temp;

    // Thermal conductivity for temperature
    std::vector<doublereal> m_diff_temp;

    // Temporary arrays for species
    std::vector<doublereal> m_temp_vol_massfraction;
    std::vector<doublereal> m_temp_diff_coeffs;
    std::vector<doublereal> m_temp_surf_coverages;

    // Temporary arrays for species production rates
    std::vector<doublereal> m_temp_production_rates;
    // Temporary reaction rates
    std::vector<doublereal> m_temp_rates_of_progress;
    // Temporary arrays for reaction enthalpy
    std::vector<doublereal> m_temp_delta_enthalpy;
    std::vector<doublereal> m_temp_species_enthalpy;

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

    std::vector<double> m_prefactor;

    doublereal m_update_temp;

    int m_update_counter;


    // Density at each grid point including boundaries
    grid_vector m_rho;

    // Diffusion Coefficients at each grid point including boundaries
    comp_vector m_diff_coeffs;

    // Pointer wcdata object
    wcdata* m_wcdata;

    std::vector<doublereal> m_fluxes_return;

    // Number of volume species
    int m_vol_sp;

    // Number of surface species
    int m_surf_sp;

    // Number of variables for the
    // ODE solver (gridpoints *variables)
    int m_nvars;

    // Small value to avoid division by zero
    double m_small;

    // Enumerations
    std::vector<vars_enum> en_vol_comp;
    std::vector<vars_enum> en_surf_comp;
    vars_enum en_temperature;



   // Creates the grid. This routine takes nx and wc_thickness
   // and buids a grid. It set the corner coordinates and the
   // interpolation cofficinets for the faces
   void createGrid();

      // Updated the fluxes for all the species and the energy equation
   void update_fluxes(double* state);

      // Uses the interpolation coefficients to interpolate the
   // input values th a cell phase
   double interpolate_values(double fxp,double val,double valw);


   // Updates the material properties
   // Diffusion coefficient, conductivity, density
   // Currently sets it to the bulk state
   void update_material_properties(double* y);

    // Takes a state vector gets the value of var at grid point
    // loc_idx
    // loc_idx is the index in the state vector, which means boundary
    // points are not included

    // Takes a state vector sets value of var at grid point
    // loc_idx
    // loc_idx is the index in the state vector, which means boundary
    // points are not included

   double bulk_massfraction(double* y,int nc);
   double bulk_temperature(double* y);

};
} /* namespace Cantera */

#endif /* SINGLEWC_H_ */
