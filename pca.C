/* Program written by Giovanni Pinamonti
  PhD student at 
  Scuola Internazionale Superiori di Studi Avanzati, Trieste, Italy
   begin june 11th 2013 */

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>   
#include <iostream>

//#include "my_malloc.h" //libraries to allocate vectors and matrices
//#include "io.h" //definitions of input and output subroutines
#include "Matrix.h" //definition of the class Matrix
#include "PrincipalComp.h"


using namespace std;

void help_display(){
  cout<<"Help:"<<endl<<
    "-f --> input file (required)"<<endl<<
    "-ntop --> number of eigenvectors to write on file"<<endl<<
    "-proj --> number of components on which project on"<<endl<<
    //    "-dist whether or not to compute the distance fluctuations (will read the traj again)"<<endl<<
    "-beads --> beads to consider (default = ALL )"<<endl<<
    "-nolast --> skip the first/last residues (default = 0 )"<<endl<<
    "-dpf --> dump partial mean square fluctuations (default = false)"<<endl<<
    "-int --> compute effective interaction matrix? (default = false)"<<endl<<
    "-dumpcov --> dump the covariance matrix? (default = false)"<<endl<<
    "-dumpdist --> dump the fluctuations of distances? (default = false)"<<endl<<
    "-ntmax --> maximum number of time steps to read (default = -1 i.e. all the trajectory)"<<endl<<
    "-ntmin --> minimum number of time steps to start from (default = -1 i.e. all the trajectory)"<<endl<<
    "-covfor --> file to compute forces covariance (DEFAULT = none)"<<endl<<
    "-name --> (DEFAULT = intput file name)"<<endl<<
    "-h help"<<endl;;
  return;
}


int main(int argc, char **argv){//PROGRAM MAIN
  char file_name[200];
  char name[200];
  char covfor_name[200];
  int n_comp=0;
  int ntop=-1;
  int nlast=0;
  int NTMAX=-1;
  int NTMIN=-1;
  bool COMPINT=false;
  bool DUMPCOV=false;
  bool DUMPDIST=false;
  bool DPF=false;
  bool COVFOR=false;
  bool has_name=false;
  vector<string> bead_list;  


  cout<<" --- Program PCA --- "<<endl<<endl;
  int i=1;
  bool nopar=true;
  while(i<argc){
    if(!strncmp(argv[i],"-f",2)){
      //    if(argv[i]=="-f"){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(file_name,argv[i]);
      nopar=false;
      cout<<"will analize the file "<<file_name<<endl;
    }
    else if(!strncmp(argv[i],"-proj",5)){
      //    else if(argv[i]=="-proj"){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sscanf(argv[i],"%d",&n_comp);
      cout<<"will project on "<<n_comp<<" directions"<<endl;
    }
    else if(!strncmp(argv[i],"-ntop",5)){
      //    else if(argv[i]=="-proj"){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sscanf(argv[i],"%d",&ntop);
      cout<<"will dump the first "<<ntop<<" eigenvectors"<<endl;
    }
    else if(!strncmp(argv[i],"-int",4)){
      COMPINT=true;
      cout<<"I'm going to compute the effective interaction matrix"<<endl;
    }
    else if(!strncmp(argv[i],"-dumpcov",8)){
      DUMPCOV=true;
      cout<<"Dump the covariance matrix"<<endl;
    }
    else if(!strncmp(argv[i],"-dumpdist",9)){
      DUMPDIST=true;
      cout<<"Dump the distance fluctuations"<<endl;
    }
    else if(!strncmp(argv[i],"-dpf",4)){
      DPF=true;
      cout<<"Dump partial fluctuations"<<endl;
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


    else if(!strncmp(argv[i],"-nolast",7)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      sscanf(argv[i],"%d",&nlast);
      cout<<"won't consider the first and last "<<nlast<<" residues"<<endl;
    }    

    else if(!strncmp(argv[i],"-ntmax",6)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      sscanf(argv[i],"%d",&NTMAX);
      if (NTMAX>0) cout<<"stop after "<<NTMAX<<" steps"<<endl;
    }    
    
    else if(!strncmp(argv[i],"-ntmin",6)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      sscanf(argv[i],"%d",&NTMIN);
      if (NTMIN>0) cout<<"stop after "<<NTMIN<<" steps"<<endl;
    }    
    else if(!strncmp(argv[i],"-covfor",7)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(covfor_name,argv[i]);
      COVFOR=true;
      cout<<"will analize the fforces from file "<<covfor_name<<endl;
    }
    else if(!strncmp(argv[i],"-name",5)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(name,argv[i]);
      cout<<"\"My name is "<<name<<".\""<<endl;
      has_name=true;
    }


    else if(!strncmp(argv[i],"-h",2)){
      //    else if(argv[i]=="-h"){
      help_display();
      return 0;
    }
    else{
      cout<<"WARNING: input parameter n. "<<i<<" "<<argv[i]<<" ---> invalid parameter"<<endl;
    }
    i=i+1;
  }
  if(nopar){cout<<"ERROR: insufficient number of parameters"<<endl;
    help_display;
    return 0;
  }

  if(!has_name) sprintf(name,file_name);

  cout<<"Beads to consider are: ";
  for(int i=0; i<bead_list.size();++i){
    cout<<bead_list.at(i)<<" ";}
  cout<<endl<<endl;

  BeadList lista(bead_list);

  if(nlast>0){
    Structure3d temp_structure(file_name);
    vector<int> res_list;  
    int temp_size=temp_structure.getSizeRes();
    if( temp_size<2*nlast+1 ){ cout<<"ERROR: the structure in"<<file_name<<" has less than "<<2*nlast+1<<" residues! I can't apply the option -nolast"<<endl; return 0;}
    else{
      for(int i=0;i<temp_size;i++){
	if((i>=nlast)&&(i<temp_size-nlast))
	  res_list.push_back( (*(temp_structure.getBeads().begin())).getResNum()+i );
      }
      
      
      lista.setResNum(res_list);

      cout<<"Residues to consider are: ";
      for(int i=0; i<res_list.size();++i){
	cout<<res_list.at(i)<<" ";}
      cout<<endl<<endl;
      

    }
  }

  ///TODO still have to set the bead list from here
  //PrincipalComp cacca(file_name);
  PrincipalComp cacca;
  cacca.setBeadList(bead_list);
  if (NTMAX>0) cacca.set_ntmax(NTMAX);
  if (NTMIN>0) cacca.set_ntmin(NTMIN);
  cacca.readTrajectory(file_name);
  //  PrincipalComp cacca(file_name,lista);
  //PrincipalComp cacca(file_name,bead_list);
  //  cacca.readTrajectory(file_name,bead_list);
  //  cacca(file_name);
  cout<<"traj was read from -> "<<file_name<<endl;

  cacca.Diagonalize();
  cacca.setNTOP(ntop);
  cacca.Dump(name);

  if(DUMPCOV){
    cacca.DumpCov(name);
  }
  if(DUMPDIST){
    cacca.dumpDistFluc(name);
  }

  if(DPF) cacca.dumpMSF_partial(name);

  if(COMPINT){
    cout<<endl<<endl<<"### NOW I DO STUFFS WITH THE INTERACTION MATRIX ###"<<endl<<file_name<<endl;
    cacca.ComputeInt();
    cacca.DumpInt(name);
  }

  if(n_comp>0){
    cout<<endl<<"### PROIETTIAMO "<<endl;
    cacca.Project(n_comp,file_name,name);
  }


  if(COVFOR){
    cacca.readForces(covfor_name);
  }//endif covfor

  return 0;

}//END PROGRAM MAIN
 
