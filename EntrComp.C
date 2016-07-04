/* Program EntrComp:                */
/* reads a pdb trajectory and       */
/* for each frame constructs        */
/* the ENM and compute the entropy  */
/* PhD student at SISSA, Trieste    */
/* March 27th, 2014                 */


#include <iostream>

#include "Structure3d.h"
#include "ElasticNet.h"
#include <cmath>


#define K_BOLT 1.3806488e-23 // J K^-1
#define N_AVOG 6.02214129e23 // mol^-1
#define TEMP 300             // K

using namespace std;
void help_display(){
  cout<<"Help:"<<endl<<
    "-f --> input file  (required) format .pdb"<<endl<<
    "-o --> output file (DEFAULT = entropy.dat)"<<endl<<
    "-p --> parameters file (DEFAULT = PARAMS_RIBOGM.DAT)"<<endl<<
    //    "-beads --> beads to consider (default = C1' C2 P )"<<endl<<
    "-h help"<<endl;;
  return;
}


int main(int argc, char*argv[]){
  // ++ DEFAULT PARAMETERS
  char file_name[200], outf_name[200], par_name[200];
  sprintf(outf_name,"entropy.dat");
  sprintf(par_name,"PARAMS_RIBOGM.DAT");
  //  vector<string> bead_list;
  // bead_list.push_back("P");
  // bead_list.push_back("C2");
  // bead_list.push_back("C4'");

  // ++ READ PARAMETERS
  cout<<" --- Program ENTRCOMP --- "<<endl<<endl;
  int i=1;
  int nopar=0;
  while(i<argc){
    if(!strncmp(argv[i],"-f",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(file_name,argv[i]);
      nopar+=1;
      cout<<"will analize the trajectory "<<file_name<<endl;
    }
    else if(!strncmp(argv[i],"-o",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(outf_name,argv[i]);
      cout<<"will write the output on the file "<<outf_name<<endl;
    }

    else if(!strncmp(argv[i],"-p",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(par_name,argv[i]);
      cout<<"will read parameters from the file "<<par_name<<endl;
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
  if(nopar<1){cout<<"ERROR: insufficient number of parameters"<<endl;
    help_display;
    return 0;
  }

  // +++++++++++++++++
  // +++ START PROGRAM
  // +++++++++++++++++

  // +++ OPEN the INPUTfile +++
  FILE* fp;
  if ((fp = fopen (file_name, "r")) == NULL) {
    fprintf (stderr,"Could not open file %s.\n.Exiting.\n", file_name);
    exit (1);    }

  // +++ OPEN the OUTPUTfile +++
  ofstream fout;
  fout.open(outf_name);

  int error=0;
  int Nsteps=0;

  ElasticNet ENM;
  ENM.setFast();
  ENM.readParameters(par_name);
  error=ENM.readPDBFile(fp);
  int N=ENM.getN();

  cout<<"#taglia : "<<N<<endl;
  cout<<"#error : "<<error<<endl;
  
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  while(error==0){                                     // * * *   
    ++Nsteps;					       // * * * 
    if(Nsteps%1000==0) cout<<"Step "<<Nsteps<<endl;    // * * * 
    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    //    ENM.setStructure(structure);
    ENM.constructContactMap();
    ENM.constructIntMat();
    ENM.Solve();
    double s=-0.5*ENM.getIntMatrix().GetLogDet()*N_AVOG*K_BOLT;
    fout<<Nsteps<<" "<<s<<" "<<TEMP*s<<endl;
      
    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    error=ENM.readPDBFile(fp);  	       // * * *
  }//end while  			 	       // * * *
  //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  fclose(fp);
  fout.close();
  
  cout<<"Done!"<<endl<<endl;
  
  return 0;
}

#undef TEMP
#undef K_BOLT
#undef N_AVOG
