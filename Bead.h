/*
 * Bead.h
 *
 *  Created on: Jul 19, 2012
 *      Author: gpolles
 *
 *  Bead class.
 *  Each bead has 3 coordinates and pointers to its neighbors
 *  and to a domain.
 *  It have also some information from the pdb file
 *  as chainId, and betaFactor
 *
 */

#ifndef BEAD_H_
#define BEAD_H_

#include <vector>

#include "Vector3d.h"
#include "common.h"


class RigidDomain;

class Bead
{
public:
 Bead() : _index(0), _num_neighbors(0),_num_covBonded(0), _domain(NULL),
           _isEdge(false), _coordinates(Vector3d(0,0,0)),
    _atomType(""),_resNum(0),_resName(""),_hetatm(false),
           _betaFactor(0.0), _pdbPosition(0), _chainId(0) {}
  virtual ~Bead(){}

  double distance(Bead* b){
    return _coordinates.distance(b->_coordinates);
  }

  double distance(Bead& b){
    return _coordinates.distance(b._coordinates);
  }

  double distanceSQ(Bead* b){
    return _coordinates.distanceSQ(b->_coordinates);
  }

  double distanceSQ(Bead& b){
    return _coordinates.distanceSQ(b._coordinates);
  }

  RigidDomain* getDomain() const {
    return _domain;
  }

  void setDomain(RigidDomain* domain){
    _domain = domain;
  }

  int getIndex() const {
    return _index;
  }
  
  bool getHet() {
    return _hetatm;
  }

  void setHet(bool HET) {
    _hetatm=HET;
  }

  void setIndex(int index) {
    _index = index;
  }

  std::vector<Bead*>& getNeighbors() {
    return _neighbors;
  }

  void setNeighbors(std::vector<Bead*> neighbors) {
    _neighbors = neighbors;
  }

  int getNumNeighbors() const {
    return _num_neighbors;
  }

  void setNumNeighbors(int numNeighbors) {
    _num_neighbors = numNeighbors;
  }

  void addNeighbor(Bead* neighbor){
    _neighbors.push_back(neighbor);
    _num_neighbors = _neighbors.size();
  }



  /************************/
std::vector<Bead*>& getCovBonded() {
    return _covBonded;
  }

  void setCovBonded(std::vector<Bead*> covBonded) {
    _covBonded = covBonded;
  }

  int getNumCovBonded() const {
    return _num_covBonded;
  }

  void setNumCovBonded(int numCovBonded) {
    _num_covBonded = numCovBonded;
  }  

  void addCovBond(Bead* covbond){
    _covBonded.push_back(covbond);
    _num_covBonded = _covBonded.size();
  }

  /***********************/


  bool isEdge(){
    return _isEdge;
  }

  void setEdge(bool b){
    _isEdge =b;
  }

  double getBetaFactor() const {
    return _betaFactor;
  }

  void setBetaFactor(double betaFactor) {
    _betaFactor = betaFactor;
  }

  //#####################################
  std::string getAtomType() {
    return _atomType;
  }
  void setAtomType(std::string type) {
    _atomType=type;
  }

  int getBeadType() {
    return _beadType;
  }
  void setBeadType(int type) {
    _beadType=type;
  }
  
  std::string getResName() {
    return _resName;
  }
  void setResName(std::string name) {
    _resName=name;
  }

  int getResNum() {
    return _resNum;
  }
  void setResNum(int resnum) {
    _resNum=resnum;
  }
  //####################################

  


  const Vector3d& getCoordinates() const {
    return _coordinates;
  }

  void setCoordinates(const Vector3d& coordinates) {
    _coordinates = coordinates;
  }

  int getPdbPosition() const {
    return _pdbPosition;
  }

  void setPdbPosition(int pdbPosition) {
    _pdbPosition = pdbPosition;
  }

  int getChainId() const {
    return _chainId;
  }

  void setChainId(int chainId) {
    _chainId = chainId;
  }

private:
  int _index;
  int _num_neighbors;
  int _num_covBonded;
  RigidDomain* _domain;
  std::vector<Bead*> _neighbors;
  std::vector<Bead*> _covBonded;
  bool _isEdge;
  Vector3d _coordinates;
  double _betaFactor;
  int _pdbPosition;
  bool _hetatm;


  int _beadType; // 0 -> P 
                 // 1 -> sugar 
                 // 2 -> base
  std::string _atomType;
  std::string _resName;
  int _resNum;

  int _chainId;
};

#endif /* BEAD_H_ */
