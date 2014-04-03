/*
 * wcdata.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: vonrickenbach
 */

#include "cantera/kinetics/wcdata.h"
#include "cantera/Cantera.h"
#include "cantera/kinetics/ImplicitSurfChem_wc.h"
#include <fstream>
#include <sstream>


namespace Cantera{
wcdata::wcdata(Cantera::ImplicitSurfChem_wc* wc_obj,wcdata* data):
  m_wc_obj(wc_obj)

{
   m_vol_sp  = wc_obj->get_vol_sp();
   m_surf_sp = wc_obj->get_surf_sp();
   m_nx      = wc_obj->get_nx();

   m_temperature        = data->m_temperature;
   m_surf_massfractions = data->m_surf_massfractions;
   m_vol_massfractions  = data->m_vol_massfractions;

}

wcdata::wcdata(Cantera::ImplicitSurfChem_wc* wc_obj):
   m_wc_obj(wc_obj)

{
   m_vol_sp  = wc_obj->get_vol_sp();
   m_surf_sp = wc_obj->get_surf_sp();
   m_nx      = wc_obj->get_nx();

   m_vol_massfractions.resize(m_nx*m_vol_sp);
   m_surf_massfractions.resize(m_nx*m_surf_sp);
   m_temperature.resize(m_nx);

	m_wc_obj->get_state_object(*this);

}

wcdata::grid_vec& wcdata::get_temperature()  {
   return m_temperature;
}

wcdata::grid_vec& wcdata::get_surf_massfractions()  {
   return m_surf_massfractions;
}

wcdata::grid_vec& wcdata::get_vol_massfractions()  {
   return m_vol_massfractions;
}

const wcdata::grid_vec& wcdata::get_temperature() const {
   return m_temperature;
}

const wcdata::grid_vec& wcdata::get_surf_massfractions() const {
   return m_surf_massfractions;
}

const wcdata::grid_vec& wcdata::get_vol_massfractions() const {
   return m_vol_massfractions;
}

wcdata::~wcdata() {

	// TODO Auto-generated destructor stub
}

void wcdata::write_data() const {
	std::ofstream myfile;
	std::stringstream ss;
	std::string temp_string;

	ss << "grid_" << m_nx << ".dat";

	temp_string = ss.str();
	myfile.open(temp_string.c_str());

    grid_vec::const_iterator it;
    int idx;

    for (idx=0,it = m_vol_massfractions.begin();
       it!=m_vol_massfractions.end();++it,++idx){
        myfile << *it << " ";
       if (!((idx+1) % m_nx)) myfile << std::endl;
    }


    for (idx=0,it = m_surf_massfractions.begin();
       it!=m_surf_massfractions.end();++it,++idx){
       myfile << *it << " ";
       if (!((idx+1) % m_nx)) myfile << std::endl;
    }

    for (idx=0,it = m_temperature.begin();
       it!=m_temperature.end();++it,++idx){
       myfile << *it-1000 << " ";
       if (!((idx+1) % m_nx)) myfile << std::endl;
    }

    myfile.close();
}

} // Namespace cantera

