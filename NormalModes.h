/*
 * NormalModes.h
 *
 *  Created on: Jul 5, 2012
 *      Author: gpolles
 */

#ifndef NORMALMODES_H_
#define NORMALMODES_H_
#include <vector>
#include "Vector3d.h"
#include "common.h"



class NormalModes
{
public:
  NormalModes();
  virtual ~NormalModes();

  void 		readFromFile(const std::string& basename, size_t n, size_t dim, bool printProgress);

  real_t 	getEigenvalue(size_t i) const;
  const std::vector<Vector3d>& getEigenvector(size_t i) const;

  size_t 	getNumModes() const { return _numModes; }

  void clear();

private:
  std::vector<real_t> _eigenvalues;
  std::vector<std::vector<Vector3d> > _eigenvectors;
  size_t _numModes;

};


#endif /* NORMALMODES_H_ */
