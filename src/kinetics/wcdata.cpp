/*
 * wcdata.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: vonrickenbach
 */

#include "cantera/kinetics/wcdata.h"
#include "cantera/Cantera.h"
#include "cantera/kinetics/SingleWc.h"
#include <fstream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <iostream>
#include <cstdlib>


namespace Cantera{
wcdata::wcdata(Cantera::SingleWc* wc_obj,wcdata* data):
  m_wc_obj(wc_obj)

{
   m_vol_sp  = wc_obj->get_vol_sp();
   m_surf_sp = wc_obj->get_surf_sp();
   m_nx      = wc_obj->get_nx();

   m_temperature        = data->m_temperature;
   m_surf_coverages     = data->m_surf_coverages;
   m_vol_massfractions  = data->m_vol_massfractions;

}

wcdata::wcdata(Cantera::SingleWc* wc_obj):
   m_wc_obj(wc_obj)
{
   initialize_arrays();
	m_wc_obj->get_state_object(*this);

}

wcdata::wcdata(Cantera::SingleWc* wc_obj,int x_idx, int step):
   m_wc_obj(wc_obj)
{
   initialize_arrays();
   read_data(x_idx,step);
}

void wcdata::initialize_arrays(){
   m_vol_sp  = m_wc_obj->get_vol_sp();
   m_surf_sp = m_wc_obj->get_surf_sp();
   m_nx      = m_wc_obj->get_nx();

   m_vol_massfractions.resize((m_nx+1)*m_vol_sp);
   m_surf_coverages.resize(m_nx*m_surf_sp);
   m_temperature.resize(m_nx+1);

}

wcdata::grid_vec& wcdata::get_temperature()  {
   return m_temperature;
}

wcdata::grid_vec& wcdata::get_surf_coverages()  {
   return m_surf_coverages;
}

wcdata::grid_vec& wcdata::get_vol_massfractions()  {
   return m_vol_massfractions;
}

const wcdata::grid_vec& wcdata::get_temperature() const {
   return m_temperature;
}

const wcdata::grid_vec& wcdata::get_surf_coverages() const {
   return m_surf_coverages;
}

const wcdata::grid_vec& wcdata::get_vol_massfractions() const {
   return m_vol_massfractions;
}

wcdata::~wcdata() {

	// TODO Auto-generated destructor stub
}

void wcdata::write_data(int x_idx,int step) const {
	std::ofstream myfile;
	std::stringstream ss;
	std::string temp_string;

	ss << "grid_" << x_idx <<"_" << step << ".dat";

	std::cout.precision(15);

	temp_string = ss.str();
	myfile.open(temp_string.c_str());

    grid_vec::const_iterator it;
    int idx;

    for (idx=0,it = m_vol_massfractions.begin();
       it!=m_vol_massfractions.end();++it,++idx){
       myfile << std::scientific
              << std::setprecision(std::numeric_limits<double>::digits10/4)
              << *it << " ";
       if (!((idx+1) % (m_nx+1))) myfile << std::endl;
    }


    for (idx=0,it = m_surf_coverages.begin();
       it!=m_surf_coverages.end();++it,++idx){
       myfile << std::scientific
              << std::setprecision(std::numeric_limits<double>::digits10/4)
              << *it << " ";
       if (!((idx+1) % (m_nx))) myfile << 0  << std::endl;
    }

    for (idx=0,it = m_temperature.begin();
       it!=m_temperature.end();++it,++idx){
       myfile << std::scientific
              << std::setprecision(std::numeric_limits<double>::digits10/4)
              << *it << " ";
       if (!((idx+1) % (m_nx+1))) myfile << std::endl;
    }

    myfile.close();
}

void wcdata::read_data(int x_idx, int step) {
	std::ifstream myfile;
	std::stringstream ss;
	std::string temp_string;
	std::string header;
	double dummy;

	ss << "grid_" << x_idx <<"_" << step << ".dat";

	temp_string = ss.str();
	myfile.open(temp_string.c_str());

    grid_vec::iterator it;
    int idx;
 //   std::getline(myfile,header);

    for (idx=0,it = m_vol_massfractions.begin();
       it!=m_vol_massfractions.end();++it,++idx){
        myfile >> *it;
    }


    for (idx=0,it = m_surf_coverages.begin();
       it!=m_surf_coverages.end();++it,++idx){
       myfile >> *it;
       if (!((idx+1) % (m_nx))) myfile >> dummy;
    }

    for (idx=0,it = m_temperature.begin();
       it!=m_temperature.end();++it,++idx){
       myfile >> *it;
    }

    myfile.close();
}

} // Namespace cantera

