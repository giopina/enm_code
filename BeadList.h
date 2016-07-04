/*
 * BeadList.h
 *
 * Created on Mar 25, 2014
 * Author: Giovanni Pinamonti 
 *
 * BeadList class.
 *
 */

#ifndef BEADLIST_H_
#define BEADLIST_H_

#include <vector>
#include "Bead.h"
#include <algorithm>

class BeadList
{

 private:
  std::vector<std::string> _atomTypes;
  std::vector<int> _resNumbers;
  bool _noH;
  bool _n_1_9;
  
 public:
  //creatori
 BeadList():_noH(true), _n_1_9(false) {};
  
 BeadList(std::vector<std::string> namelist):_noH(true), _n_1_9(false) {
    if(namelist.size()>0){   
      for (int i=0;i<namelist.size();++i){

	if((namelist.at(i)).compare("N1-9")==0){
	  _n_1_9=true;
	  namelist.erase(namelist.begin()+i);
	  cout<<")))))))))))) STRANO BEAD N1-9 ATTIVATO (((((((((((("<<endl;
	}
	if((namelist.at(i)).compare("ALL")==0){
	  namelist.erase(namelist.begin()+i);
	  cout<<")))))))))))) ALL ATOMS (((((((((((("<<endl;
	}
      }
    }
    _atomTypes=namelist;
    return;
  }

 BeadList(std::vector<int> numlist):_noH(true), _n_1_9(false) {
    _resNumbers=numlist;
  }
  /*
  //funzioni
  */
  void setResNum(std::vector<int> numlist){
    //    cout<<"cacca"<<endl;
    _resNumbers=numlist; 
  } 
  void setAtomTypes(std::vector<std::string> types){
    //  std::vector<std::string> 
    _atomTypes=types;
  } 
  
  
  inline bool check(Bead b){
    
    //    cout<<"cacca"<<endl;
    //    bool lacacca=true;
    if((_atomTypes.size()>0)||(_n_1_9)){
      bool lacacca=false;
      std::string s1=b.getAtomType();
      for (int i=0;i<_atomTypes.size();++i){
	std::string s2=_atomTypes.at(i);
	std::replace( s2.begin(), s2.end(), 'p', '\'');
	//	cout<<s2<<endl;
	if (s1.compare(s2)==0) lacacca=true;
      }
      if(_n_1_9){
	if( ( (b.getResName()).compare("G")==0) || 
	    ((b.getResName()).compare("A")==0) ){
	  if(s1.compare("N9")==0) lacacca=true;  
	}
	else{
	  if(s1.compare("N1")==0) lacacca=true;  
	}
      }
      if(!lacacca) return false;
    }
    
    if(_resNumbers.size()>0){
      bool lacacca=false;
      int r=b.getResNum();
      for (int i=0;i<_resNumbers.size();++i){
	if (r==_resNumbers.at(i)) {lacacca=true;}
      }
      if(!lacacca) return false;
    }
    
    //    cout<<_noH<<endl;
    if(_noH){
      std::string s=b.getAtomType();
      //      cout<<s<<endl;
      if (s.find("H") != std::string::npos) {
	//	cout<<"found!"<<endl;
	return false;
      } 
    }
    
    return true;
  }  
};
  
  
#endif /* BEAD_H_ */
