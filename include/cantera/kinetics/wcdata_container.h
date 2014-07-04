/*
 * wcdata_container.h
 *
 *  Created on: Mar 28, 2014
 *      Author: vonrickenbach
 */

#ifndef WCDATA_CONTAINER_H_
#define WCDATA_CONTAINER_H_

#include<vector>

namespace Cantera{
  class wcdata;
  class SingleWc;
}

namespace Cantera{

class wcdata_container {
public:
	wcdata_container(int istorf, Cantera::SingleWc* integ,wcdata* data);

	wcdata* get_wcdata(int iistr1_nb,int ii);
	void set_wcdata(int iistr1_nb,int ii,wcdata* data);
	virtual ~wcdata_container();
protected:
	typedef std::vector<wcdata*> wcdata_vec;
	wcdata_vec m_data_vec;
	int m_istorf;
};

}

#endif /* WCDATA_CONTAINER_H_ */
