#include "cantera/kinetics/SingleWc.h"
#include "cantera/kinetics/wcdata.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/ImplicitSurfChem_wc.h"
#include "cantera/thermo/SurfPhase.h"
#include<iostream>
#include<fstream>

namespace Cantera {

SingleWc::SingleWc(ImplicitSurfChem_wc* surf_chem,InterfaceKinetics* k, Transport* t
           ,double h,double h_temp
           ,double wc_thickness, double area_to_volume
           ,double porosity, double tortuosity
           ,double d_p, double lambda_solid
           ,int nx, bool with_energy,int x_idx, int nxcells, double L_r, double vel, double A_V,bool from_file
           ,double mintemp, double maxtemp, double trate, double rhocp, double rhocp_st, bool inf_ext_mt
           ,double bulk_pressure, double heat_source, double cell_ratio):
      m_kin(k),
      m_transport(t),
      m_wc_coefficient_temp(h_temp),
      m_wc_coefficient_in(h),
      m_wc_thickness(wc_thickness),
      m_nx(nx),
      m_area_to_volume(area_to_volume),
      m_with_energy(with_energy),
      m_porosity(porosity),
      m_lambda_solid(lambda_solid),
      m_tortuosity(tortuosity),
      m_dp(d_p),
      m_x_idx(x_idx),
      m_surf_chem(surf_chem),
      m_L_r(L_r),
      m_vel(vel),
      m_nxcells(nxcells),
      m_A_V(A_V),
      m_from_file(from_file),
      m_mintemp(mintemp),
      m_maxtemp(maxtemp),
      m_trate(trate),
      m_rhocp(rhocp),
      m_rhocp_st(rhocp_st),
      m_inf_ext_mt(inf_ext_mt),
      m_bulk_pressure(bulk_pressure),
      m_heat_source(heat_source),
      m_cell_ratio(cell_ratio)

{
      // Set number of corners
      m_nco    = m_nx+1;
      // Set points where the variables are stored
      m_nx_var = m_nx+2;
      m_small = 1E-50;

      // Initialize the grid
      createGrid();

      // Not more than tow phases allowed
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


      // Initialize the inflow with the composition
      // of the gas object
      m_inflow_comp.resize(m_vol_sp,0.0);
      m_gas_phase->getMassFractions(&m_inflow_comp.front());

      if (m_with_energy) m_inflow_temp = m_gas_phase->temperature();

      // Allocate temporary arrays
      m_prefactor.resize(m_nx_var);
      m_temp_vol_massfraction.resize(m_vol_sp,0);
      m_temp_diff_coeffs.resize(m_vol_sp,0);
      m_temp_surf_coverages.resize(m_surf_sp,0);
      m_temp_production_rates.resize(m_surf_sp+m_vol_sp,0);

      m_bulk_diff_coeffs.resize(m_vol_sp,0);


      // Initialize arrays related to energy equation
      if (m_with_energy){
         m_temp_rates_of_progress.resize(m_kin->nReactions());
         m_temp_delta_enthalpy.resize(m_kin->nReactions());
         m_temp_species_enthalpy.resize(m_vol_sp,0.0);
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
      if (with_energy) m_fluxes_temp.resize(m_nco,0.0);

      comp_vector::iterator comp_iter;
      for (comp_iter =  m_fluxes.begin();
           comp_iter != m_fluxes.end();
           ++comp_iter){
         comp_iter->resize(m_nco,0.0);
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

      // Initialize from file or
      // From object
      if (m_from_file){
        m_wcdata = new wcdata(this,m_x_idx,0);
      }else{
        m_wcdata = new wcdata(this);
      }

      // Set masstransfer coefficient
      m_wc_coefficient.resize(m_vol_sp,h);

      m_update_temp = -1;
      m_update_counter = 0;

   }

SingleWc::~SingleWc(){

}


void SingleWc::eval(doublereal time, doublereal* y,
                    doublereal* ydot, doublereal* p)
{
   double source = 0.0;
   double dx = 0.0;
   double dx_temp;

   update_material_properties(y);
   update_fluxes(y);

   for (int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
      dx = m_x_co[g_idx_p+1] - m_x_co[g_idx_p];

      setCanterafromState(y,g_idx_p);
      m_kin->getNetProductionRates(&m_temp_production_rates.front());
      for (int nc =0;nc<m_vol_sp;++nc){
         source = 0.0;
         source = (m_fluxes[nc][g_idx_p] - m_fluxes[nc][g_idx_p+1])/(dx+m_small);

         if (source != source) source = 0.0;

         source = source + m_temp_production_rates[nc]
                         * m_gas_phase->molecularWeight(nc)
                         * m_area_to_volume;

         if(source != source)source = 0.0;

         setStateVar(ydot,source,en_vol_comp[nc],g_idx_p+1,m_x_idx);

      }

      for (int nc =0;nc<m_surf_sp;++nc){
         source = 0.0;
         source = source + m_temp_production_rates[nc+m_vol_sp]
                          / m_surface_phase->siteDensity()
                          * m_area_to_volume* m_wc_thickness;

         if(source != source) source = 0.0;
         setStateVar(ydot,source,en_surf_comp[nc],g_idx_p,m_x_idx);
      }

      if (m_with_energy){
         if (g_idx_p+1 == m_nx and m_rhocp_st > m_small){
            dx_temp = m_rhocp_st;
         }else{
            dx_temp = dx;
         }

         source = 0.0;
         source = (m_fluxes_temp[g_idx_p] - m_fluxes_temp[g_idx_p+1])/ dx_temp;

         m_kin->getNetRatesOfProgress(&m_temp_rates_of_progress.front());
         m_kin->getDeltaEnthalpy(&m_temp_delta_enthalpy.front());
         m_gas_phase->getEnthalpy_RT(&m_temp_species_enthalpy.front());
         for (int nr=0;nr < m_kin->nReactions();++nr){
            source -= m_area_to_volume *
                      m_temp_rates_of_progress[nr] * m_temp_delta_enthalpy[nr];
         }
         //for (int nc=0;nc < m_vol_sp;++nc){
         //   source -= m_area_to_volume *
         //             m_temp_production_rates[nc] * m_temp_species_enthalpy[nc];
         //}
         source += m_heat_source / (m_A_V * m_wc_thickness + m_small);
         // Taking into account the thermal inertia of the solid
         if (source != source) source = 0.0;
         source = source/(m_rhocp);
         setStateVar(ydot,source,en_temperature,g_idx_p+1,m_x_idx);
      }

      // Set ydot to zero value in extra cell
      for (int nc =0;nc<m_surf_sp;++nc){
         source = 0.0;
         setStateVar(ydot,source,en_surf_comp[nc],m_nx,m_x_idx);
      }
   }


   double bulk_y_dot = 0;
   double state_left = 0.0;
   double cp;
   cp = m_gas_phase->cp_mass();

   for (int nc =0;nc<m_vol_sp;++nc){
      if (!m_inf_ext_mt){
         if (m_x_idx == 0){
            state_left = m_inflow_comp[nc];
         }
         else{
            state_left = getStateVar(y,en_vol_comp[nc],0,m_x_idx-1);
         }
         bulk_y_dot = m_vel* (state_left-getStateVar(y,en_vol_comp[nc],0,m_x_idx))/(m_L_r/m_nxcells)- m_A_V* m_fluxes[nc][0]/m_bulk_density;

      }else{
         bulk_y_dot = 0.0;
      }

      if (bulk_y_dot != bulk_y_dot){
         bulk_y_dot = 0.0;
      }
      setStateVar(ydot,bulk_y_dot,en_vol_comp[nc],0,m_x_idx);
   }

   if (m_with_energy){
      if (!m_inf_ext_mt){
         if (m_x_idx == 0){
            state_left = m_inflow_temp;
         }else{
            state_left = getStateVar(y,en_temperature,0,m_x_idx-1);
         }

         bulk_y_dot = m_vel* (state_left-getStateVar(y,en_temperature,0,m_x_idx))/(m_L_r/m_nxcells)- m_A_V* m_fluxes_temp[0]/(m_bulk_density*cp);

         setStateVar(ydot,bulk_y_dot,en_temperature,0,m_x_idx);
      }else{
         setStateVar(ydot,0.0,en_temperature,0,m_x_idx);
      }
   }

}

void SingleWc::setCanterafromState(doublereal* y,int g_idx_p){
   double temperature;
   for (int nc =0;nc<m_vol_sp;++nc){
      m_temp_vol_massfraction[nc]  = getStateVar(y,en_vol_comp[nc],g_idx_p+1,m_x_idx);
   }

   for (int nc =0;nc<m_surf_sp;++nc){
      m_temp_surf_coverages[nc] = getStateVar(y,en_surf_comp[nc],g_idx_p,m_x_idx);
   }

   if (m_with_energy) {
      temperature = std::max(getStateVar(y,en_temperature,g_idx_p+1,m_x_idx),100.0);
   }else{
      temperature = bulk_temperature(y);
   }

   m_gas_phase->setState_TPY(temperature,m_bulk_pressure,&m_temp_vol_massfraction.front());
   m_surface_phase->setMoleFractions(&m_temp_surf_coverages.front());
   m_surface_phase->setTemperature(temperature);
}

void SingleWc::createGrid(){
   //compute dx in every cell
   m_x.resize(m_nx);
   m_dx.resize(m_nx);
   m_x_co.resize(m_nco-1);
   m_fx.resize(m_nco);

   grid_vector::iterator vec_it;
   grid_vector::iterator vec_it_co;
   grid_vector::iterator vec_it_dx;
   grid_vector::iterator vec_it_x;


   double current_dx = m_wc_thickness *(1-m_cell_ratio)/(1-pow(m_cell_ratio,m_nx));
   for(vec_it  = m_dx.begin();
      vec_it != m_dx.end();
      ++vec_it){
      *vec_it = current_dx;
      current_dx *= m_cell_ratio;
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

double SingleWc::interpolate_values(double fxp,double valp,double valw){
   return fxp * valp + (1-fxp) * valw;
}

void SingleWc::set_bulk_temperature(double temperature){
   m_bulk_temperature = temperature;
}

void SingleWc::update_fluxes(double* state){

   double Y_1, T_1;
   double prefactor_inf, prefactor_bulk, prefactor_mix;
   double dx_inf;

   for(int nc=0;nc<m_vol_sp;++nc){

      // Compute rho * diff coeff for all the cells
      for(int g_idx_p=0; g_idx_p < m_nx_var; ++g_idx_p){
           m_prefactor[g_idx_p] =  m_rho[g_idx_p]
                          * m_diff_coeffs[g_idx_p][nc];
      }

      // Compute west fluxes for all the cells
      for(int g_idx_w=1; g_idx_w < m_nx; ++g_idx_w){

         //prefactor
         m_fluxes[nc][g_idx_w] = interpolate_values(m_fx[g_idx_w]
                                    ,m_prefactor[g_idx_w],m_prefactor[g_idx_w+1]);

         // Multiply with component gradient
         m_fluxes[nc][g_idx_w] =  m_fluxes[nc][g_idx_w] *
                        (getStateVar(state,en_vol_comp[nc],g_idx_w,m_x_idx) -
                         getStateVar(state,en_vol_comp[nc],g_idx_w+1,m_x_idx)) /
                         (m_x[g_idx_w+1] - m_x[g_idx_w] + m_small);
      }


      //D_b*D_i*h*rho_b*rho_i*(-Y_1 + Y_b)/(D_b*dx*h*rho_b + D_i*rho_i)
      // Apply the boundary condition for the first cell

      // Special case for no wash coat
      if (m_nx == 1){
         dx_inf = 0.0;
      }else{
         dx_inf = m_x[1] - m_x[0];
      }
      Y_1 = getStateVar(state,en_vol_comp[nc],1,m_x_idx);
      prefactor_inf  = interpolate_values(m_fx[0],m_prefactor[0],m_prefactor[1]);
      if (m_inf_ext_mt){
         m_fluxes[nc][0] =(prefactor_inf *
                   (bulk_massfraction(state,nc)- Y_1 )) / (dx_inf+ m_small);

      }else{
         prefactor_bulk = m_bulk_diff_coeffs[nc] * m_bulk_density
                          * m_wc_coefficient[nc] * dx_inf;
         prefactor_mix  = prefactor_inf * m_bulk_density *m_wc_coefficient[nc]
                           *m_bulk_diff_coeffs[nc];

         m_fluxes[nc][0] =(prefactor_mix *
                   (bulk_massfraction(state,nc)- Y_1 ))/
                  (prefactor_bulk + prefactor_inf + m_small);

      }
   }

   // Fluxes for the energy equation
   if (m_with_energy){
      for(int g_idx_w=1; g_idx_w < m_nx; ++g_idx_w){

         m_fluxes_temp[g_idx_w] = interpolate_values(m_fx[g_idx_w]
                                                    ,m_diff_temp[g_idx_w],m_diff_temp[g_idx_w+1]);
         m_fluxes_temp[g_idx_w] = m_fluxes_temp[g_idx_w] *
                      (getStateVar(state,en_temperature,g_idx_w,m_x_idx) -
                       getStateVar(state,en_temperature,g_idx_w+1,m_x_idx)) /
                       (m_x[g_idx_w+1] - m_x[g_idx_w] + m_small);

      }

      // Boundary condition
      T_1 = getStateVar(state,en_temperature,1,m_x_idx);
      prefactor_inf    = interpolate_values(m_fx[0],m_diff_temp[0],m_diff_temp[1]);
      if (m_inf_ext_mt){
         m_fluxes_temp[0] = prefactor_inf * (bulk_temperature(state) - T_1) / (dx_inf+ m_small);
      }else{
         prefactor_bulk   = m_bulk_diff_temp * m_wc_coefficient_temp * dx_inf;
         prefactor_mix    = prefactor_inf * m_wc_coefficient_temp * m_bulk_diff_temp;
         m_fluxes_temp[0] = prefactor_mix * (bulk_temperature(state) - T_1)/
                          (prefactor_bulk + prefactor_inf + m_small);
      }
   }
}

void SingleWc::getInitialConditions(doublereal* y)
{
   set_state(y,*m_wcdata);
}

void SingleWc::update_material_properties(double* y){
   double temperature;
   double knudsen_diff_coeff;


   temperature = bulk_temperature(y);

   if (std::abs(temperature -m_update_temp) < 1.0 and m_update_counter < 100){
      ++m_update_counter;
      return;
   }
   m_update_counter = 0;
   m_update_temp = temperature;

   // Update Bulk properties
   for (int nc =0;nc<m_vol_sp;++nc){
      m_temp_vol_massfraction[nc] = bulk_massfraction(y,nc);
   }
   m_gas_phase->setState_TPY(temperature,m_bulk_pressure,&m_temp_vol_massfraction.front());
   m_transport->getMixDiffCoeffs(&m_bulk_diff_coeffs.front());
   m_bulk_diff_temp = m_transport->thermalConductivity();
   m_bulk_density   = m_gas_phase->density();

   if (m_wc_coefficient_in < 0){
      //set_mt_coefficient_mladenov();
      set_mt_coefficient_kelvin();
   }

   // Loop over all the cells in the was coat
   for (int loc_idx=0;loc_idx<m_nx_var;++loc_idx){

      // Compute the surface massfraction based on the flux condition for the first cell
      if (loc_idx == 0){
         for (int nc =0;nc<m_vol_sp;++nc){
            m_temp_vol_massfraction[nc] = -m_fluxes[nc][0]/
                                       (m_bulk_diff_coeffs[nc]*
                                        m_wc_coefficient[nc] * m_bulk_density +m_small)
                                        + bulk_massfraction(y,nc);
         }
         if (m_with_energy){
         temperature = -m_fluxes_temp[0]/
                     (m_bulk_diff_temp *
                      m_wc_coefficient_temp + m_small)
                    + bulk_temperature(y);
         }

      }

      // Copy from inside at the boundary since the flux is zero
      else if (loc_idx == m_nx_var-1){
         for (int nc =0;nc<m_vol_sp;++nc){
            m_temp_vol_massfraction[nc] = getStateVar(y,en_vol_comp[nc],loc_idx-1,m_x_idx);
         }
         if (m_with_energy){
            temperature = getStateVar(y,en_temperature,loc_idx-1,m_x_idx);
         }
      }

      else{
         for (int nc =0;nc<m_vol_sp;++nc){
            m_temp_vol_massfraction[nc]  = getStateVar(y,en_vol_comp[nc],loc_idx,m_x_idx);
         }
         if (m_with_energy){
            temperature = getStateVar(y,en_temperature,loc_idx,m_x_idx);
         }
      }


      m_gas_phase->setState_TPY(temperature,m_bulk_pressure
                              ,&m_temp_vol_massfraction.front());

      m_rho[loc_idx]  = m_gas_phase->density();


      // Bulk diffusion coefficients
      m_transport->getMixDiffCoeffs(&m_temp_diff_coeffs.front());
      for (int nc=0;nc<m_vol_sp;++nc){
         // Knudsen diffusion coefficients
         knudsen_diff_coeff = m_dp/3.0 * sqrt(abs(8.0*GasConstant*temperature/
                                           (Pi * m_gas_phase->molecularWeight(nc))+m_small));

         // Effective diffusion coefficient
         m_diff_coeffs[loc_idx][nc] = m_porosity
                                      /(m_tortuosity + m_small) *
                                      (1.0/(1.0/(knudsen_diff_coeff + m_small) +1.0/(m_temp_diff_coeffs[nc] + m_small) + m_small));
      }

      if (m_with_energy){
         m_diff_temp[loc_idx] = m_porosity * m_transport->thermalConductivity()
                               + (1.0-m_porosity) * m_lambda_solid;
      }
   }

}

void SingleWc::set_mt_coefficient_mladenov(){
   double dx = m_L_r/m_nxcells;
   double gz_prefac, loc_gz;
   gz_prefac = (m_x_idx * dx + dx/2) / (pow(0.001,2) * m_vel);
   for (int nc=0;nc <m_vol_sp;++nc){
      loc_gz = gz_prefac * m_bulk_diff_coeffs[nc];
      m_wc_coefficient[nc] = (3.675 + 8.827*std::pow((1000.*loc_gz),-0.545)
                                    *std::exp(-48.2*loc_gz))*1000.0;
   }
}

void SingleWc::set_mt_coefficient_kelvin(){
   double dx = m_L_r/m_nxcells;
   double L_Kelvin = 0.0023;
   double x = (m_x_idx * dx + dx/2);
   for (int nc=0;nc <m_vol_sp;++nc){
      if (x <= L_Kelvin){
        m_wc_coefficient[nc] = 18.07/L_Kelvin;
      }else{
        m_wc_coefficient[nc] = 14.3/L_Kelvin;
      }
   }
}

void SingleWc::get_state(wcdata& data) {

   double* solution = m_surf_chem->m_integ->solution();
   update_material_properties(solution);
   update_fluxes(solution);

   for (int nc=0;nc<m_vol_sp;++nc){
      for(int g_idx_p=0;g_idx_p<m_nx+1;++g_idx_p){
         data.get_vol_massfractions().at((m_nx+1)*nc+g_idx_p) = getStateVar(solution,en_vol_comp[nc],g_idx_p,m_x_idx);
      }
   }

   for (int nc=0;nc<m_surf_sp;++nc){
      for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
         data.get_surf_coverages().at(m_nx*nc+g_idx_p) = getStateVar(solution,en_surf_comp[nc],g_idx_p,m_x_idx);
      }
   }

   for(int g_idx_p=0;g_idx_p<m_nx+1;++g_idx_p){
      if (m_with_energy){
          data.get_temperature().at(g_idx_p) = getStateVar(solution,en_temperature,g_idx_p,m_x_idx);
      }else{
          data.get_temperature().at(g_idx_p)   = bulk_temperature(solution);
      }
   }

   for (int nc=0;nc<m_vol_sp;++nc){
      for(int g_idx_p=0;g_idx_p<m_nx+1;++g_idx_p){
         data.get_fluxes().at((m_nx+1)*nc+g_idx_p) = m_fluxes[nc][g_idx_p];
      }
   }
}

void SingleWc::get_state_object(wcdata& data)  {

   m_gas_phase->getMassFractions(&m_temp_vol_massfraction.front());
   for (int nc=0;nc<m_vol_sp;++nc){
      data.get_vol_massfractions().at((m_nx+1)*nc) = m_temp_vol_massfraction[nc];
      for(int g_idx_p=0;g_idx_p<m_nx+1;++g_idx_p){
         data.get_vol_massfractions().at((m_nx+1)*nc+g_idx_p) = m_temp_vol_massfraction[nc];
      }
   }

   m_surface_phase->getMoleFractions(&m_temp_surf_coverages.front());
   for (int nc=0;nc<m_surf_sp;++nc){
      for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
         data.get_surf_coverages().at(m_nx*nc+g_idx_p) = m_temp_surf_coverages[nc];
      }
   }

   double temperature = m_gas_phase->temperature();
   for(int g_idx_p=0;g_idx_p<m_nx+1;++g_idx_p){
      data.get_temperature().at(g_idx_p) = temperature;
   }
}

void SingleWc::set_state(double* state, const wcdata& data){
   double temp_value;
   double sum_val = 0.0;
   for(int g_idx_p=0;g_idx_p<m_nx+1;++g_idx_p){
      sum_val = 0.0;
      for (int nc=0;nc<m_vol_sp;++nc){
         temp_value = std::max(0.0,std::min(1.0,data.get_vol_massfractions().at((m_nx+1)*nc+g_idx_p)));
         sum_val += temp_value;
      }
      for (int nc=0;nc<m_vol_sp;++nc){
         temp_value = std::max(0.0,std::min(1.0,data.get_vol_massfractions().at((m_nx+1)*nc+g_idx_p)))/sum_val;
         setStateVar(state,temp_value,en_vol_comp[nc],g_idx_p,m_x_idx);
      }
   }

   if (m_inf_ext_mt){
      for (int nc=0;nc<m_vol_sp;++nc){
         setStateVar(state,m_inflow_comp[nc],en_vol_comp[nc],0,m_x_idx);
      }
   }



   for(int g_idx_p=0;g_idx_p<m_nx;++g_idx_p){
      sum_val = 0.0;
      for (int nc=0;nc<m_surf_sp;++nc){
         temp_value = std::max(0.0,std::min(1.0,data.get_surf_coverages().at(m_nx*nc+g_idx_p)));
         sum_val += temp_value;
      }
      for (int nc=0;nc<m_surf_sp;++nc){
         temp_value = std::max(0.0,std::min(1.0,data.get_surf_coverages().at(m_nx*nc+g_idx_p)))/sum_val;
         setStateVar(state,temp_value,en_surf_comp[nc],g_idx_p,m_x_idx);
      }
   }

   for (int nc=0;nc<m_surf_sp;++nc){
      setStateVar(state,1.0,en_surf_comp[nc],m_nx,m_x_idx);
   }


   if (m_with_energy){
      for(int g_idx_p=0;g_idx_p<m_nx+1;++g_idx_p){
            temp_value = data.get_temperature().at(g_idx_p);
            setStateVar(state,temp_value,en_temperature,g_idx_p,m_x_idx);
      }
   }else{
      m_bulk_temperature = data.get_temperature().at(0);
   }


}

size_t SingleWc::neq(){return m_nvars*(m_nx+1);}

void SingleWc::write_wc(int step){
   get_state(*m_wcdata);
   m_wcdata->write_data(m_x_idx,step);
}

double SingleWc::bulk_massfraction(double* y,int nc){
   if (m_inf_ext_mt){
      return m_inflow_comp[nc];
   }
   return getStateVar(y,en_vol_comp[nc],0,m_x_idx);

}

double SingleWc::bulk_temperature(double* y){
   if (m_inf_ext_mt and m_with_energy){
      return m_inflow_temp;
   }
   if (m_with_energy){
      m_bulk_temperature = getStateVar(y,en_temperature,0,m_x_idx);
 //  }else if (m_mintemp > 0.0){
 //     double deltat = (m_maxtemp-m_mintemp)/m_trate;
 //     if (time < deltat){
 //        m_bulk_temperature = m_mintemp+m_trate*time;
 //     }else{
 //        m_bulk_temperature = m_maxtemp-m_trate*(time-deltat);
 //     }
   }
   // Temperature is fixed
   return m_bulk_temperature;
}


} /* namespace Cantera */
