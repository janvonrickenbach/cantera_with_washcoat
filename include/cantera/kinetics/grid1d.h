/*
 * grid1d.h
 *
 *  Created on: Mar 17, 2014
 *      Author: vonrickenbach
 */

#ifndef GRID1D_H_
#define GRID1D_H_

#include <vector>

class grid_1d {
public:
	grid_1d();
	virtual ~grid_1d();

protected:
	std::vector<double> dx;
	std::vector<double> corners;
	std::vector<double> finterw;
	void interpolate_to_face(const std::vector<double>& center,
			                       std::vector<double> face);
};

#endif /* GRID1D_H_ */
