/* rewritten by Giovanni Pinamonti, Monday March 17th 2014 */
//============================================================================
// Name        : generateConfigs.cpp
// Author      : Guido Polles
// Version     :
// Copyright   : 
// Description : 
//============================================================================

#include <iostream>
#include <string>
#include <cmath>
#include <boost/random.hpp>

#include "Structure3d.h"
#include "PrincipalComp.h"
//#include "Matrix.h"
using namespace std;

void help_display(){
  cout<<"Help:"<<endl<<
    "-f --> input file  (required) format .pdb"<<endl<<
    "-temp --> normalized temperature (default =1)"<<endl<<
    "-nframes --> number of frames (default = 10000)"<<endl<<
    "-ntop --> number of modes to consider (default = all)"<<endl<<
    "-beads --> beads to consider (default = [ C1' C2 P ] )"<<endl<<
     "-h help"<<endl;;
  return;
}


int main(int argc, char* argv[]) {
  boost::mt19937 rng(time(0));
  int NTOP=-1;
  int NFRAMES=10000;
  double fattore=1.;
  char file_name[200];
  vector<string> bead_list;  
  bead_list.push_back("C1'");
  bead_list.push_back("C2");
  bead_list.push_back("P");

  cout<<" --- Program GenConf --- "<<endl<<endl;
  int i=1;
  int nopar=0;
  while(i<argc){
    if(!strncmp(argv[i],"-f",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      sprintf(file_name,argv[i]);
      nopar+=1;
      cout<<"will generate a trajectory for the file "<<file_name<<endl;
    }
    else if(!strncmp(argv[i],"-temp",5)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      if(sscanf(argv[i],"%lf",&fattore)!=1){cout<<"ERROR: in parameter -temp"<<endl; help_display(); return 0;};
      cout<<"normalized temperature = "<<fattore<<endl;
    }    
    else if(!strncmp(argv[i],"-nframes",8)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      if(sscanf(argv[i],"%d",&NFRAMES)!=1){cout<<"ERROR: in parameter -frames"<<endl; help_display(); return 0;};
      cout<<"N of frame = "<<NFRAMES<<endl;
    } 
    else if(!strncmp(argv[i],"-ntop",5)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      if(sscanf(argv[i],"%d",&NTOP)!=1){cout<<"ERROR: in parameter -ntop"<<endl; help_display(); return 0;};
      cout<<"N of modes to consider = "<<NTOP<<endl;
    } 
   
    else if(!strncmp(argv[i],"-beads",6)){
      bead_list.clear();
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      while(strncmp(argv[i],"[",1)!=0){cout<<"ERROR: in parameter -beads"<<endl; help_display(); return 0;};
      i=i+1;
      while(strncmp(argv[i],"]",1)!=0){
	char name[4];
	sprintf(name,argv[i]);
	bead_list.push_back(name);
	i=i+1;
      }
      if(bead_list.size()==0) {
	cout<<"WARNING: no beads declared!"<<endl;;
      }
    }    


    else if(!strncmp(argv[i],"-h",2)){
      help_display();
      return 0;
    }
    else{
      cout<<"WARNING: input n. "<<i<<" "<<argv[i]<<" ---> invalid parameter"<<endl;
    }
    i=i+1;
  }
  if(nopar<1){cout<<"ERROR: insufficient number of parameters"<<endl;
    help_display;
    return 0;
  }


  //  cout<<"cacca"<<endl;
  Structure3d ref_struc(0);

  cout<<"Beads to consider are: ";
  for(int i=0; i<bead_list.size();++i){
    cout<<bead_list.at(i)<<" ";}
  cout<<endl<<endl;
  ref_struc.setBeadList(bead_list);


  cout<<"cacca"<<endl;
  ref_struc.readFromPDBFile(file_name);
  cout<<"cacca"<<endl;
  Structure3d new_struc=ref_struc;

  int n_beads=ref_struc.getSize();
  int nmodes=3*n_beads;
  //  if(NTOP>0) nmodes=NTOP;

  cout<<"N. beads = "<<n_beads<<endl;

  PrincipalComp nm(n_beads);
  nm.readFromFile(file_name,nmodes);

  vector<Vector3d> coords(ref_struc.getSize());

  // +++ OPEN the OUTPUT file +++
  char oname[1024]; 
  ofstream dumpfile;
  sprintf(oname,"%s_trajectory.pdb",file_name);
  dumpfile.open(oname);
  // +++ +++ +++ +++ +++ +++ +++

  int NFPRINT=NFRAMES/10;
  if(NFPRINT<1) NFPRINT=1;

  for (size_t f = 0; f < NFRAMES; ++f){
    if(f%NFPRINT==0) cout<<f<<" frames generated"<<endl;
    //    cout<<f<<" frames generated"<<endl;
    for (size_t i = 0; i < ref_struc.getSize(); ++i){
      coords[i]=ref_struc.getBead(i).getCoordinates();
    }
    //cout<<"cacca"<<endl;
    if((NTOP<1)||(NTOP>nmodes-6)) NTOP=nmodes-6;
    for (int m = 0; m < NTOP-6; ++m) {
      double s=nm.CovGetEigenval(nmodes-m-1);
      //      double s=nm.CovGetEigenval(m);
      if(s<0.0001){ cout<<"genConf *** WARNING: negative eigenvalue -> "<<m<<endl; s=abs(s);}
      s=sqrt(s);
      boost::normal_distribution<> nd(0.0,s);
      
      boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >  rnd(rng, nd);
      double width = rnd()/5.0*fattore;
      for (size_t i = 0; i < ref_struc.getSize(); ++i){
	coords[i] += nm.CovGetEigenvec(nmodes-m-1)[i]*width;
      }
    }
    
    //printf("%ld\ncoordinates generated from normal modes\n",ref_struc.getSize());
    for (size_t i = 0; i < ref_struc.getSize(); ++i){
      char atomo[4];
      
      new_struc.getBead(i).setCoordinates(coords[i]);
    }
    new_struc.dumpPDBFile(dumpfile);

  }//enddo frames

  return 0;
}
