/* Program fast_rwsip:                  */
/* compute a gaussian network model */
/* for a given RNA molecule.        */
/* then compute RWSIP with a given  */
/* PCA set.                         */
/* Created by Giovanni Pinamonti    */
/* PhD student at SISSA, Trieste    */
/* October, 2014                    */


#include <iostream>

#include "Matrix.h"
#include "Structure3d.h"
#include "ElasticNet.h"
#include "PrincipalComp.h"

using namespace std;


void help_display(){
  cout<<"Help:"<<endl<<
    "-f --> input file  (required) format .pdb"<<endl<<
    "-pca --> PCA file  (required) format .pdb"<<endl<<
    "-par --> parameters file (DEFAULT = PARAMS_RIBOGM.DAT)"<<endl<<
    "-o --> output file (DEFAULT = $input_file+_rwsip.dat)"<<endl<<
    //    "-beads --> beads to consider (default = [ C1' C2 P ] )"<<endl<<
     "-h help"<<endl;;
  return;
}


int main(int argc, char*argv[]){
  char file_name[200];
  char par_name[200];
  char pca_name[200];
  char oname[1024]; 
  sprintf(par_name,"PARAMS_RIBOGM.DAT");
  sprintf(oname,"%s_rwsip.dat",file_name);

  int i=1;
  int nopar=0;
  while(i<argc){
    if(!strncmp(argv[i],"-f",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      sprintf(file_name,argv[i]);
      nopar+=1;
    }
    else if(!strncmp(argv[i],"-pca",4)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(pca_name,argv[i]);
    }
    else if(!strncmp(argv[i],"-par",4)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(par_name,argv[i]);
    }
    else if(!strncmp(argv[i],"-o",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(oname,argv[i]);
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

  cout<<"will read parameters from the file "<<par_name<<endl;
  cout<<"will create the ENM on the file "<<file_name<<endl;
  cout<<"will read PCA from the file "<<pca_name<<endl;
  cout<<"will write the output on the file "<<oname<<endl;

  ElasticNet ENM;
  cout<<"Reading the parameters..."<<endl;
  ENM.readParameters(par_name);
  cout<<"Reading the pdb file..."<<endl;
  ENM.readPDBFile(file_name);
  cout<<"Constructing the gaussian model..."<<endl;
  cout<<"...the contact map..."<<endl;
  ENM.constructContactMap();
  cout<<"...and the interaction matrix!"<<endl<<endl;
  ENM.constructIntMat();
  cout<<"Let's go!"<<endl;
  ENM.Solve();
  cout<<"...done!"<<endl<<endl;
  
  Structure3d ref_struc=ENM.getStructure();
  int n_beads=ref_struc.getSize();
  int nmodes=3*n_beads;

  PrincipalComp PCA(n_beads);
  PCA.readFromFile(pca_name,nmodes);
  //  int NTOP=nmodes;

  cout<<"### RWSIP ###"<<endl;
  //### prima opzione: leggo e calcolo in contemporanea ###
  double RWSIP=0.0;
  double NORM=0.0;
  for (int m = 6; m < nmodes; ++m) {
    double eval_ENM=ENM.getCovMatrix().GetEigenval(m);
    //    cout<<"ENM: "<<eval_ENM<<endl;
    vector<Vector3d> evec_ENM=ENM.getCovMatrix().GetEigenvec(m);
    for (int n = 6; n < nmodes; ++n) {
      double eval_PCA=PCA.getCovMat().GetEigenval(n);
      vector<Vector3d> evec_PCA=PCA.getCovMat().GetEigenvec(n);
      double scal_prod=0;
      for(int i=0; i<n_beads; ++i){
	scal_prod+=evec_ENM[i].dot(evec_PCA[i]);
      }
      double cacca=eval_PCA*eval_ENM;
      RWSIP+=cacca*scal_prod*scal_prod;
    }
    double eval_PCA_i=PCA.getCovMat().GetEigenval(m);
    //  cout<<"PCA: "<<eval_PCA_i<<endl;
    NORM+=eval_PCA_i*eval_ENM;
  }
  RWSIP/=NORM;
  RWSIP=sqrt(RWSIP);
  cout<<endl<<endl
      <<"#####################################"<<endl
      <<"#     The RWSIP is "<<RWSIP<<endl
      <<"#####################################"<<endl<<endl;
  cout<<"===============Done!================="<<endl;


  // +++ OPEN the OUTPUT file +++
  ofstream fout;
  fout.open(oname);
  // +++ +++ +++ +++ +++ +++ +++
  fout<<RWSIP<<endl;
  
  fout.close();

  return 0;
}
