/*
 * wcdata.h
 *
 *  Created on: Mar 28, 2014
 *      Author: vonrickenbach
 */

#ifndef WCDATA_H_
#define WCDATA_H_

#include <vector>
namespace Cantera{
   class ImplicitSurfChem_wc;

class wcdata {
public:
   wcdata(Cantera::ImplicitSurfChem_wc* wc_obj,wcdata* data);
   wcdata(Cantera::ImplicitSurfChem_wc* wc_obj);
   virtual ~wcdata();

   Cantera::ImplicitSurfChem_wc* m_wc_obj;

   typedef std::vector<double> grid_vec ;

   grid_vec& get_temperature();

   grid_vec& get_surf_massfractions();

   grid_vec& get_vol_massfractions();

   const grid_vec& get_temperature() const;

   const grid_vec& get_surf_massfractions() const;

   const grid_vec& get_vol_massfractions() const;

   void write_data() const;


protected:
grid_vec m_vol_massfractions;
   grid_vec m_surf_massfractions;
   grid_vec m_temperature;
   int m_vol_sp;
   int m_surf_sp;
   int m_nx;


};
}

#endif /* WCDATA_H_ */
