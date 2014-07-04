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
   class SingleWc;

class wcdata {
public:
   wcdata(Cantera::SingleWc* wc_obj,wcdata* data);
   wcdata(Cantera::SingleWc* wc_obj,int x_idx, int step);
   wcdata(Cantera::SingleWc* wc_obj);
   virtual ~wcdata();

   Cantera::SingleWc* m_wc_obj;

   typedef std::vector<double> grid_vec ;

   grid_vec& get_temperature();

   grid_vec& get_surf_coverages();

   grid_vec& get_vol_massfractions();

   const grid_vec& get_temperature() const;

   const grid_vec& get_surf_coverages() const;

   const grid_vec& get_vol_massfractions() const;

   void write_data(int x_idx, int step) const;
   void read_data(int x_idx, int step);


protected:
grid_vec m_vol_massfractions;
   grid_vec m_surf_coverages;
   grid_vec m_temperature;
   int m_vol_sp;
   int m_surf_sp;
   int m_nx;
   void initialize_arrays();

};
}

#endif /* WCDATA_H_ */
