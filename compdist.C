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

#include "Matrix.h" //definition of the class CMatrix
#include "Structure3d.h"

using namespace std;


void help_display(){
  cout<<"Help:"<<endl<<
    "-f --> input trajectory file  (required) format .pdb"<<endl<<
    //    "-n --> number of frames (default = 10000)"<<endl<<
    //    "-beads --> beads to consider (default = [ C1' C2 P ] )"<<endl<<
    "-of --> output file for fluctuations (DEFAULT = none)"<<endl<<
    "-ribo --> is a riboswitch? (do I need to take away the last res because it's an ADA?) (default=false)"<<endl<<
     "-h help"<<endl;;
  return;
}

int main(int argc, char **argv){//PROGRAM MAIN
  char file_name[200];
  char oname[200];
  bool out_has_name=false;
  bool ribo=false;

  cout<<" --- Program CompDist --- "<<endl<<endl;
  cout<<" Computed the distance    "<<endl;
  cout<<" fluctuations between C2 (or C6!!!!!!)  "<<endl;
  cout<<" from a .pdb trajectory   "<<endl<<endl;
  cout<<" ------------------------ "<<endl<<endl;
  int i=1;
  int nopar=0;
  while(i<argc){
    if(!strncmp(argv[i],"-f",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      sprintf(file_name,argv[i]);
      nopar+=1;
      cout<<"will analize the file "<<file_name<<endl;
    }
    else if(!strncmp(argv[i],"-of",3)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(oname,argv[i]);
      cout<<"will write the flluctuations on the file "<<oname<<endl;
      out_has_name=true;
    }
    else if(!strncmp(argv[i],"-ribo",5)){
      cout<<" $ $ $ Riboswitch $ $ $ "<<endl;
      ribo=true;
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

  if(!out_has_name)  sprintf(oname,"%s_C2_fluc.dat",file_name);
  
  double *dist_mean, *d_sq_mean;    
  //  double *dist_mean, *d_sq_mean;    
  Structure3d structure(0);
  
  //open input file
  FILE* fin;
  if ((fin = fopen (file_name, "r")) == NULL) {
    fprintf (stderr,"CompDist:: Could not open file %s.\n.Exiting.\n", file_name);
    exit (1);    }    
  
  //check if the input trajectory is correct
  int error=0;
  int Nsteps=0;
  error=structure.readFromPDBFile(fin);
  Structure3d ref_struc=structure;
  int n_beads=structure.getSize();


  int* is_C2=i1t(n_beads);
  vector<int> indexC2;
  for (size_t i = 0; i < ref_struc.getSize(); ++i){
    if(ref_struc.getBead(i).getHet()) continue;
    if ((ref_struc.getBead(i).getAtomType()=="C2") ||
	(ref_struc.getBead(i).getAtomType()=="C6"))
      {is_C2[i]=1; indexC2.push_back(i);}// cout<<"eccolo: "<<i<<endl;}
  }
  
  cout<<"N. beads = "<<n_beads<<endl;

  if(indexC2.size()<2) {cout<<"ERROR: INSUFFICIENT NUMBER OF C2 BEADS!!!"<<endl; return 0;}

  vector<Vector3d> coords(indexC2.size());
    
  dist_mean=d1t(n_beads);
  d_sq_mean=d1t(n_beads);
  
  cout<<"PrincipalComp: reading trajectory..."<<endl;
  while(error==0){
    ++Nsteps;
    if(Nsteps%1000==0) cout<<"Frame: "<<Nsteps<<endl;

     for (size_t i = 0; i < indexC2.size(); ++i){
       int i_bead=indexC2.at(i);
       coords[i] = structure.getBead(i_bead).getCoordinates();
     }
     //  $  $  $  $  $  $  $  $
     
    // #  #  #  #  #  #  # 
     for (size_t i = 0; i < indexC2.size(); ++i){ //### qui calcolo le distanze
       double temp_dsq=coords[i].distanceSQ(coords[i+1]);
       d_sq_mean[i]+=temp_dsq;
       dist_mean[i]+=sqrt(temp_dsq);
     }
    
    error=structure.readFromPDBFile(fin);
  }//end while

  //close input file
  fclose(fin);
  cout<<"@@@ chiudo file input"<<endl;
  cout<<"@@@ "<<Nsteps<<" frame letti."<<endl;

  vector<double> fluct;
  int n_C2=0;
  if(ribo) n_C2=indexC2.size()-2;//NB: -2 per escludere ADA del riboswitch, altrimenti -1 !!
  else n_C2=indexC2.size()-1;

  for (size_t i = 0; i < n_C2; ++i){ 
    d_sq_mean[i]/=Nsteps;
    dist_mean[i]/=Nsteps;
    fluct.push_back(d_sq_mean[i]-(dist_mean[i]*dist_mean[i]));
  }
  /// +++ OPEN the OUTPUT file +++
  //  if(out_has_name){
    ofstream fout;
    fout.open(oname);
    // ### qui calcolo le distanze e le printo ###
    for (size_t i = 0; i < n_C2; ++i){ //NB: -2 per escludere ADA del riboswitch, altrimenti -1 !!
      fout<<"X "<<i<<" "<<fluct.at(i)<<" "<<indexC2.at(i)<<" "<<ref_struc.getBead(indexC2.at(i)).getResName()<<endl;
    }
    fout.close();
    //  }// +++ +++ +++ +++ +++ +++ +++  
  
  cout<<"...done!"<<endl<<endl;  
  
  return 0;

}//END PROGRAM MAIN
