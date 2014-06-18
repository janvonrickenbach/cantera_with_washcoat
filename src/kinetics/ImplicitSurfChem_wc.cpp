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

namespace Cantera{
// Constructor
ImplicitSurfChem_wc::ImplicitSurfChem_wc(InterfaceKinetics* k
                                        ,Transport* t
                                        ,double h,double h_temp
                                        ,double wc_thickness, double area_to_volume
                                        ,double porosity, double tortuosity
                                        ,double d_p, double lambda_solid
                                        ,double atol, double rtol
                                        ,int nx, bool with_energy):
    m_kin(k),
    m_transport(t),
    m_wc_coefficient(h),
    m_wc_coefficient_temp(h_temp),
    m_wc_thickness(wc_thickness),
    m_nx(nx),
    m_area_to_volume(area_to_volume),
    m_with_energy(with_energy),
    m_atol(atol),
    m_rtol(rtol),
    m_porosity(porosity),
    m_lambda_solid(lambda_solid),
    m_tortuosity(tortuosity),
    m_dp(d_p),
    FuncEval(),
    m_integ(0),
    m_maxstep(0.0)


{

   // Set number of corners
   m_nco    = m_nx+1;
   // Set points where the variables are stored
   m_nx_var = m_nx+2;

   // Initialize the grid
   createGrid();

   //Initial integrator
    m_integ = newIntegrator("CVODE");
    // use backward differencing, with a full Jacobian computed
    // numerically, and use a Newton linear iterator
    m_integ->setMethod(BDF_Method);
    m_integ->setProblemType(GMRES);
    m_integ->setIterator(Newton_Iter);

    if (k->nPhases() > 2){
       std::cout << "Error more than two phases" << std::endl;
    }

    // Identify surface phase and store pointer
    m_surface_phase = (SurfPhase* )&k->thermo(k->surfacePhaseIndex());

    // Identify gas phase and store pointer
    for (int phase_idx=0;phase_idx< k->nPhases();++phase_idx){
      if (phase_idx != k->surfacePhaseIndex()){
        m_gas_phase    = &k->thermo(phase_idx);
      }
    }

    // Get number of species
    m_vol_sp  = m_gas_phase->nSpecies();
    m_surf_sp = k->nTotalSpecies() - m_vol_sp;


    // Allocate temporary arrays
    m_temp_vol_massfraction.resize(m_vol_sp,0);
    m_temp_diff_coeffs.resize(m_vol_sp,0);
    m_temp_surf_massfraction.resize(m_surf_sp,0);
    m_temp_production_rates.resize(m_surf_sp+m_vol_sp,0);

    // Initialize bulk arrays
    m_bulk_massfraction.resize(m_vol_sp,0);
    m_bulk_diff_coeffs.resize(m_vol_sp,0);


    // Initialize arrays related to energy
    if (m_with_energy){
       m_temp_rates_of_progress.resize(m_kin->nReactions());
       m_temp_delta_enthalpy.resize(m_kin->nReactions());
    }


    // Initialize enumerations
    for (int comp_idx=0;comp_idx < m_vol_sp;++comp_idx){
       vars_enum temp = static_cast<vars_enum>(en_start+comp_idx);
       en_vol_comp.push_back(temp);
    }

    for (int comp_idx=0;comp_idx < m_surf_sp;++comp_idx){
       vars_enum temp = static_cast<vars_enum>(en_start+m_vol_sp+comp_idx);
       en_surf_comp.push_back(temp);
    }

    if (m_with_energy) {
       en_temperature = static_cast<vars_enum>(en_start+m_vol_sp+m_surf_sp);
    }

    // Set number of variables
    m_nvars = m_surf_sp + m_vol_sp;

    if (with_energy) m_nvars +=1;


    // Initialize fluxes arrays
    m_fluxes.resize(m_vol_sp);

    m_fluxes_return.resize(m_vol_sp,0);
    if (with_energy) m_fluxes_return.push_back(0);

    if (with_energy) m_fluxes_temp.resize(m_nco,0.0);

    comp_vector::iterator comp_iter;
    for (comp_iter =  m_fluxes.begin();
         comp_iter != m_fluxes.end();
         ++comp_iter){
       comp_iter->resize(m_nco,0);
    }

    // Initialize arrays for density and diffusion
    // coefficients
    m_rho.resize(m_nx_var);
    m_diff_coeffs.resize(m_nx_var);
    if (with_energy) m_diff_temp.resize(m_nx_var);

    for (comp_iter =  m_diff_coeffs.begin();
         comp_iter != m_diff_coeffs.end();
         ++comp_iter){
       comp_iter->resize(m_vol_sp);
    }



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
void ImplicitSurfChem_wc::integrate(doublereal t0, doublereal t1,int maxsteps)
{
    m_integ->setMaxStepSize(t1 - t0);
    m_integ->setMaxStepSize(100*t1/maxsteps);
    m_integ->setMaxSteps(maxsteps);
    m_integ->initialize(t0, *this);
//    m_integ->setMaxOrder(1);
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
   double source = 0.0;
   double temperature;
   update_fluxes(y);
   update_material_properties(y);
   double dx;

//   for (int idx_full =0;idx_full<neq();++idx_full)
//   {
//      if (y[idx_full] < 0.0) return -1;
//   }

    for (int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
      dx = m_x_co[g_idx_p+1] - m_x_co[g_idx_p];

      for (int nc =0;nc<m_vol_sp;++nc){
           m_temp_vol_massfraction[nc]  = getStateVar(y,en_vol_comp[nc],g_idx_p);
        }
        for (int nc =0;nc<m_surf_sp;++nc){
         m_temp_surf_massfraction[nc] = getStateVar(y,en_surf_comp[nc],g_idx_p);
        }

        if (m_with_energy) {
           temperature = std::max(getStateVar(y,en_temperature,g_idx_p),100.0);
        } else{
           temperature = m_bulk_temperature;
        }
        m_gas_phase->setState_TPY(temperature,m_bulk_pressure,&m_temp_vol_massfraction.front());
        m_surface_phase->setMassFractions(&m_temp_surf_massfraction.front());
        m_surface_phase->setTemperature(temperature);
        m_kin->getNetProductionRates(&m_temp_production_rates.front());

        for (int nc =0;nc<m_vol_sp;++nc){
           source = 0.0;
           source = (m_fluxes[nc][g_idx_p] - m_fluxes[nc][g_idx_p+1])/(dx+m_small);
           if (source != source){
              source = 0.0;
           }
           source = source + m_temp_production_rates[nc]
                           * m_gas_phase->molecularWeight(nc)
                           * m_area_to_volume;
           if(source != source){
              source = 0.0;
           }
           setStateVar(ydot,source,en_vol_comp[nc],g_idx_p);
        }
        for (int nc =0;nc<m_surf_sp;++nc){
           source = 0.0;
           source = source + m_temp_production_rates[nc+m_vol_sp]
                            / m_surface_phase->siteDensity()
                            * m_area_to_volume* m_wc_thickness;

           if(source != source){
              source = 0.0;
           }
           setStateVar(ydot,source,en_surf_comp[nc],g_idx_p);
        }

      }
    if (m_with_energy){
      for (int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
          dx = m_x_co[g_idx_p+1] - m_x_co[g_idx_p];
        source = 0.0;
        source = m_fluxes_temp[g_idx_p] - m_fluxes_temp[g_idx_p+1];

        for (int nc =0;nc<m_vol_sp;++nc){
           m_temp_vol_massfraction[nc]  = getStateVar(y,en_vol_comp[nc],g_idx_p);
        }
        for (int nc =0;nc<m_surf_sp;++nc){
           m_temp_surf_massfraction[nc] = getStateVar(y,en_surf_comp[nc],g_idx_p);
        }

        temperature = getStateVar(y,en_temperature,g_idx_p);
        m_gas_phase->setState_TPY(temperature,m_bulk_pressure,&m_temp_vol_massfraction.front());
        m_surface_phase->setMassFractions(&m_temp_surf_massfraction.front());
        m_surface_phase->setTemperature(temperature);

        m_kin->getNetRatesOfProgress(&m_temp_rates_of_progress.front());
        m_kin->getDeltaEnthalpy(&m_temp_delta_enthalpy.front());
        for (int nr=0;nr < m_kin->nReactions();++nr){
           source -= m_area_to_volume *
                     m_temp_rates_of_progress[nr] * m_temp_delta_enthalpy[nr]*dx;
        }
        // Taking into account the thermal inertia of the solid
        source = source/(1E6*dx);
        setStateVar(ydot,source,en_temperature,g_idx_p);
//        double res_sum = 0.0;
//        for (int state_idx=0;state_idx < m_nvars; state_idx++){
//           res_sum += std::abs(ydot[state_idx]/y[state_idx]);
//        }
//        cout << res_sum;
      }
   }


}



void       ImplicitSurfChem_wc::setStateVar(double* state,double value,vars_enum var,int loc_idx){
   state[loc_idx*m_nvars+var] = value;
}doublereal ImplicitSurfChem_wc::getStateVar(double* state,vars_enum var,int loc_idx) const {
   return state[loc_idx*m_nvars+var];

}

void ImplicitSurfChem_wc::createGrid(){
   //compute dx in every cell
   m_x.resize(m_nx);
   m_dx.resize(m_nx);
   m_x_co.resize(m_nco-1);
   m_fx.resize(m_nco);

   grid_vector::iterator vec_it;
   grid_vector::iterator vec_it_co;
   grid_vector::iterator vec_it_dx;
   grid_vector::iterator vec_it_x;

   double ratio = 1.2;
   double current_dx = m_wc_thickness *(1-ratio)/(1-pow(ratio,m_nx));
   for(vec_it  = m_dx.begin();
      vec_it != m_dx.end();
      ++vec_it){
      *vec_it = current_dx;
      current_dx *= ratio;
   }


   // Get the corner values
   std::partial_sum(m_dx.begin(),m_dx.end(),m_x_co.begin());
   m_x_co.insert(m_x_co.begin(),0.0);

   // Get the cell_center values
   for(vec_it    = m_x.begin()
      ,vec_it_co = m_x_co.begin()
      ,vec_it_dx = m_dx.begin();
      vec_it != m_x.end();
      ++vec_it,++vec_it_co,++vec_it_dx){
      *vec_it = *vec_it_co + 0.5 * (*vec_it_dx);
   }
   m_x.insert(m_x.begin(),0.0);
   m_x.push_back(m_x_co.back());

   //compute interpolation coeffcient in every cell

   for(int idx=0;idx<m_nco;++idx){
      m_fx[idx] = (m_x[idx+1] - m_x_co[idx])/(m_x[idx+1]-m_x[idx]);
   }
}

void ImplicitSurfChem_wc::printGrid(){

   std::ofstream myfile;
   myfile.open("grid.dat");

   std::vector<double> temp_vec;

   temp_vec.clear();
   temp_vec.push_back(m_x[0]);
   for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
      temp_vec.push_back(m_x[g_idx_p+1]);
   }

   write_var(myfile,temp_vec);

   for (int nc=0;nc<m_vol_sp;++nc){
      temp_vec.clear();
      temp_vec.push_back(m_bulk_massfraction[nc]);
      for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
         temp_vec.push_back(getStateVar(m_integ->solution(),en_vol_comp[nc],g_idx_p));
      }
      write_var(myfile,temp_vec);
   }

   if (m_with_energy){
      temp_vec.clear();
      temp_vec.push_back(0);
      for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
         temp_vec.push_back(getStateVar(m_integ->solution(),en_temperature,g_idx_p)-m_bulk_temperature);
      }
      write_var(myfile,temp_vec);
   }

   myfile.close();

}

void ImplicitSurfChem_wc::write_var(std::ofstream &myfile,const std::vector<double> & vec){

   std::vector<double>::const_iterator vec_it;

   for(vec_it    = vec.begin();
      vec_it   != vec.end();
      ++vec_it){
      myfile << *vec_it << " ";
   }
    myfile << std::endl;



}

double ImplicitSurfChem_wc::interpolate_values(double fxp,double valp,double valw){
   return fxp * valp + (1-fxp) * valw;
}

void ImplicitSurfChem_wc::update_fluxes(double* state){

   int var_idx_p, var_idx_w;
   int state_idx_p, state_idx_w;
   int g_idx_p, g_idx_w;

   double Y_1, T_1;
   double prefactor_inf, prefactor_bulk, prefactor_mix;
   double dx_inf;

   for(int nc=0;nc<m_vol_sp;++nc){
      std::vector<double> prefactor;
      prefactor.resize(m_nx_var);

      // Compute rho * diff coeff for all the cells
      for(g_idx_p=0; g_idx_p < m_nx_var; ++g_idx_p){
           prefactor[g_idx_p] =  m_rho[g_idx_p]
                          * m_diff_coeffs[g_idx_p][nc];
      }

      // Compute west fluxes for all the cells
      for(g_idx_w=1; g_idx_w < m_nx; ++g_idx_w){
         state_idx_p = g_idx_w;
         state_idx_w = g_idx_w-1;
         var_idx_p = g_idx_w+1;
         var_idx_w = g_idx_w;

         //prefactor
         m_fluxes[nc][g_idx_w] = interpolate_values(m_fx[var_idx_w]
                                    ,prefactor[var_idx_w],prefactor[var_idx_p]);

         // Multiply with component gradient
         m_fluxes[nc][g_idx_w] =  m_fluxes[nc][g_idx_w] *
                        (getStateVar(state,en_vol_comp[nc],state_idx_w) -
                         getStateVar(state,en_vol_comp[nc],state_idx_p)) /
                         (m_x[var_idx_p] - m_x[var_idx_w] + m_small);
      }


      //D_b*D_i*h*rho_b*rho_i*(-Y_1 + Y_b)/(D_b*dx*h*rho_b + D_i*rho_i)
      // Apply the boundary condition for the first cell
      if (m_nx == 1){
         dx_inf = 0.0;
      }else{
         dx_inf = m_x[1] - m_x[0];
      }
      prefactor_inf = interpolate_values(m_fx[0],prefactor[0],prefactor[1]);
      prefactor_bulk = m_bulk_diff_coeffs[nc] * m_bulk_density
                       * m_wc_coefficient * dx_inf;
      prefactor_mix = prefactor_inf * m_bulk_density *m_wc_coefficient
                       *m_bulk_diff_coeffs[nc];

      Y_1 = getStateVar(state,en_vol_comp[nc],0);
      m_fluxes[nc][0] =(prefactor_mix *
                (m_bulk_massfraction[nc]- Y_1 ))/
               (prefactor_bulk + prefactor_inf + m_small);

      }

       if (m_with_energy){
        for(g_idx_w=1; g_idx_w < m_nx; ++g_idx_w){
           state_idx_p = g_idx_w;
           state_idx_w = g_idx_w-1;
           var_idx_p = g_idx_w+1;
           var_idx_w = g_idx_w;

           m_fluxes_temp[g_idx_w] = interpolate_values(m_fx[var_idx_w]
                                                      ,m_diff_temp[var_idx_w],m_diff_temp[var_idx_p]);
           m_fluxes_temp[g_idx_w] = m_fluxes_temp[g_idx_w] *
                        (getStateVar(state,en_temperature,state_idx_w) -
                         getStateVar(state,en_temperature,state_idx_p)) /
                         (m_x[var_idx_p] - m_x[var_idx_w] + m_small);

        }

        T_1 = getStateVar(state,en_temperature,0);
        prefactor_inf    = interpolate_values(m_fx[0],m_diff_temp[0],m_diff_temp[1]);
        prefactor_bulk   = m_bulk_diff_temp * m_wc_coefficient_temp * dx_inf;
        prefactor_mix    = prefactor_inf * m_wc_coefficient_temp * m_bulk_diff_temp;
        m_fluxes_temp[0] = prefactor_mix * (m_bulk_temperature - T_1)/
                         (prefactor_bulk + prefactor_inf + m_small);
       }
}


void ImplicitSurfChem_wc::getInitialConditions(doublereal t0, size_t lenc,
        doublereal* y)
{
   set_state(y,*m_wc_data);

}

void ImplicitSurfChem_wc::update_material_properties(double* y){
   double temperature;
   double knudsen_diff_coeff;

   temperature = m_bulk_temperature;
   for (int loc_idx=0;loc_idx<m_nx_var;++loc_idx){

     // Compute the surface massfraction based on the flux condtion
     if (loc_idx == 0){
       for (int nc =0;nc<m_vol_sp;++nc){
          m_temp_vol_massfraction[nc] = -m_fluxes[nc][0]/
                                     (m_bulk_diff_coeffs[nc]*
                                      m_wc_coefficient * m_bulk_density +m_small)
                                      + m_bulk_massfraction[nc];
       }
       if (m_with_energy){
       temperature = -m_fluxes_temp[0]/
                   (m_bulk_diff_temp *
                    m_wc_coefficient_temp + m_small)
                  + m_bulk_temperature;
       }

     }

     // Use the same value at the boundary since the flux is zero
     else if (loc_idx == m_nx_var-1){
       for (int nc =0;nc<m_vol_sp;++nc){
          m_temp_vol_massfraction[nc] = getStateVar(y,en_vol_comp[nc],loc_idx-2);
     }
       if (m_with_energy){
          temperature = getStateVar(y,en_temperature,loc_idx-2);
       }
     }

     else{
       for (int nc =0;nc<m_vol_sp;++nc){
       m_temp_vol_massfraction[nc]  = getStateVar(y,en_vol_comp[nc],loc_idx-1);
      }
      if (m_with_energy){
         temperature = getStateVar(y,en_temperature,loc_idx-1);
      }
     }

    m_gas_phase->setState_TPY(temperature,m_bulk_pressure
                           ,&m_temp_vol_massfraction.front());

     m_rho[loc_idx]       = m_gas_phase->density();
     if (m_with_energy){
        m_diff_temp[loc_idx] = m_porosity * m_transport->thermalConductivity()
                            + (1.0-m_porosity) * m_lambda_solid;
     }

     m_transport->getMixDiffCoeffs(&m_temp_diff_coeffs.front());
     for (int nc=0;nc<m_vol_sp;++nc){
      knudsen_diff_coeff = m_dp/3.0 * sqrt(abs(8.0*GasConstant*temperature/
                                        (Pi * m_gas_phase->molecularWeight(nc))+m_small));
      m_diff_coeffs[loc_idx][nc] = m_porosity
                                  /(m_tortuosity + m_small) *
                                  (1.0/(1.0/(knudsen_diff_coeff + m_small) +1.0/(m_temp_diff_coeffs[nc] + m_small) + m_small));
    }
   }

}


void ImplicitSurfChem_wc::set_bulk_from_state(){
   m_gas_phase->getMassFractions(&m_bulk_massfraction.front());
   m_bulk_pressure    = m_gas_phase->pressure();
   m_bulk_density     = m_gas_phase->density();
   m_bulk_temperature = m_gas_phase->temperature();
   m_transport->getMixDiffCoeffs(&m_bulk_diff_coeffs.front());
   m_bulk_diff_temp = m_transport->thermalConductivity();

}

void ImplicitSurfChem_wc::set_state_from_bulk(){
   m_gas_phase->setState_TPY(m_bulk_temperature
                          ,m_bulk_pressure
                          ,&m_bulk_massfraction.front());

}

void ImplicitSurfChem_wc::get_state(wcdata& data) const{

   for (int nc=0;nc<m_vol_sp;++nc){
      data.get_vol_massfractions().at((m_nx+1)*nc)   = m_bulk_massfraction[nc];
      for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
         data.get_vol_massfractions().at((m_nx+1)*nc+g_idx_p+1) = getStateVar(m_integ->solution(),en_vol_comp[nc],g_idx_p);
      }
   }

   for (int nc=0;nc<m_surf_sp;++nc){
      for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
         data.get_surf_massfractions().at(m_nx*nc+g_idx_p) = getStateVar(m_integ->solution(),en_surf_comp[nc],g_idx_p);
      }
   }

   data.get_temperature().at(0)   = m_bulk_temperature;
   if (m_with_energy){
      for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
         data.get_temperature().at(g_idx_p+1) = getStateVar(m_integ->solution(),en_temperature,g_idx_p);
      }
   }
}

void ImplicitSurfChem_wc::get_state_object(wcdata& data){

   m_gas_phase->getMassFractions(&m_temp_vol_massfraction.front());
   for (int nc=0;nc<m_vol_sp;++nc){
      data.get_vol_massfractions().at((m_nx+1)*nc) = m_temp_vol_massfraction[nc];
      for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
         data.get_vol_massfractions().at((m_nx+1)*nc+g_idx_p+1) = m_temp_vol_massfraction[nc];
      }
   }

   m_surface_phase->getMassFractions(&m_temp_surf_massfraction.front());
   for (int nc=0;nc<m_surf_sp;++nc){
      for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
         data.get_surf_massfractions().at(m_nx*nc+g_idx_p) = m_temp_surf_massfraction[nc];
      }
   }

   double temperature = m_gas_phase->temperature();
   for(int g_idx_p=0;g_idx_p<m_nx+1;++g_idx_p){
      data.get_temperature().at(g_idx_p) = temperature;
   }
}


void ImplicitSurfChem_wc::set_state(double* state, const wcdata& data){
   double temp_value;
   double sum_val = 0.0;
   for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
      sum_val = 0.0;
      for (int nc=0;nc<m_vol_sp;++nc){
         temp_value = std::max(0.0,std::min(1.0,data.get_vol_massfractions().at((m_nx+1)*nc+g_idx_p+1)));
         sum_val += temp_value;
      }
      for (int nc=0;nc<m_vol_sp;++nc){
         temp_value = std::max(0.0,std::min(1.0,data.get_vol_massfractions().at((m_nx+1)*nc+g_idx_p+1)))/sum_val;
         setStateVar(state,temp_value,en_vol_comp[nc],g_idx_p);
      }
   }


   for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
      sum_val = 0.0;
      for (int nc=0;nc<m_surf_sp;++nc){
         temp_value = std::max(0.0,std::min(1.0,data.get_surf_massfractions().at(m_nx*nc+g_idx_p)));
         sum_val += temp_value;
      }
      for (int nc=0;nc<m_surf_sp;++nc){
         temp_value = std::max(0.0,std::min(1.0,data.get_surf_massfractions().at(m_nx*nc+g_idx_p)))/sum_val;
         setStateVar(state,temp_value,en_surf_comp[nc],g_idx_p);
      }
   }


   if (m_with_energy){
      for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
            temp_value = data.get_temperature().at(g_idx_p+1);
            setStateVar(state,temp_value,en_temperature,g_idx_p);
      }
   }

}

void ImplicitSurfChem_wc::set_wcdata(wcdata* wc_data_obj){
   m_wc_data = wc_data_obj;
}

void ImplicitSurfChem_wc::get_fluxes(double* y){
   for (int nc=0;nc<m_vol_sp;++nc){
      y[nc] = m_fluxes[nc][0];
   }
   if (m_with_energy){
      y[m_vol_sp] = m_fluxes_temp[0];
   }
}

} // namespace Cantera
