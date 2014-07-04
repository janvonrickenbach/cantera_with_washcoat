/*
 * wcdata_container.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: vonrickenbach
 */

#include "cantera/kinetics/wcdata_container.h"
#include "cantera/kinetics/wcdata.h"

namespace Cantera{
   class SingleWc;


wcdata_container::wcdata_container(int istorf,Cantera::SingleWc* integ,wcdata* data):
	m_istorf(istorf){
	for (int ind=0;ind<istorf;++ind){
		m_data_vec.push_back(new wcdata(integ,data));
	}
}


wcdata* wcdata_container::get_wcdata(int iistr1_nb,int ii) {
	return m_data_vec.at(iistr1_nb+ii);
}

void wcdata_container::set_wcdata(int iistr1_nb, int ii,wcdata* data) {
	m_data_vec[iistr1_nb+ii] = data;
}

wcdata_container::~wcdata_container() {
	for (wcdata_vec::iterator it=m_data_vec.begin();
		 it!= m_data_vec.end();
		 ++it){
		delete(*it);
	}
	m_data_vec.clear();
}


}




