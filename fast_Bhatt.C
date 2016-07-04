/* Program fast_Bhatt:                  */
/* compute a gaussian network model */
/* for a given RNA molecule.        */
/* then compute Bhatt with a given  */
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


  //++++++++++++Bhatt+++++++++++++
  //  PCA.getCovMat().dumpFullMatrix(string(oname)+"_cov");

  double norm_C_ENM=ENM.getCovMatrix().GetTrace();
  double norm_C_PCA=PCA.getCovMat().GetTrace();
    ENM.getCovMatrix().Times(1./norm_C_ENM);
  PCA.getCovMat().Times(1./norm_C_PCA);
  cout<<"creo D"<<endl;
  CovMatrix D=PCA.getCovMat();
  D.Plus(ENM.getCovMatrix());
  D.Times(0.5);
  D.Decompose(); // NB: tutto e' ordinato con gli eval crescenti! (ma perche'??!)
 
  // calcolo il numero di autovalori per avere il *nu* per cento di RMSD
  double nu_frac=0.95; //fraction of RMSD to be reproduced by the eigenvalues
  double somma=0.0;
  for(int i=0; i<nmodes; i++){
    somma+=D.GetEigenval(nmodes-i-1);
  }
  double partsomma=0.0;
  int NTOP=0;
  for(int i=0;i<nmodes;i++){
    NTOP++;
    partsomma+=D.GetEigenval(nmodes-i-1);
    if(partsomma>nu_frac*somma) break;
  }
  cout<<"NTOP="<<NTOP<<endl;
  NTOP/=3;
  NTOP*=3;
  cout<<"NTOP corretto="<<NTOP<<endl;
  
  cout<<"Reducing"<<endl;
  ENM.getCovMatrix().Reduce(NTOP,D);
  PCA.getCovMat().Reduce(NTOP,D);
  ENM.getCovMatrix().Decompose();
  PCA.getCovMat().Decompose();
  //  cout<<"DIM1="<<PCA.getCovMat().GetSize()<<endl;
  CovMatrix Dnew=PCA.getCovMat();
  //  cout<<"DIM2="<<Dnew.GetSize()<<endl;
  Dnew.Plus(ENM.getCovMatrix());
  Dnew.Times(0.5);
  Dnew.Decompose();

  double log_det_C_ENM=0.0;
  double log_det_C_PCA=0.0;
  double log_det_D_new=0.0;
  for(int i=0;i<NTOP;++i){
    //    cout<<<<endl;
    cout<<ENM.getCovMatrix().GetEigenval(NTOP-i-1)<<" "<<
      PCA.getCovMat().GetEigenval(NTOP-i-1)<<" "<<
      Dnew.GetEigenval(NTOP-i-1)<<endl;

    log_det_C_ENM+=log(ENM.getCovMatrix().GetEigenval(NTOP-i-1));
    log_det_C_PCA+=log(PCA.getCovMat().GetEigenval(NTOP-i-1));
    log_det_D_new+=log(Dnew.GetEigenval(NTOP-i-1));
    //    if(ENM.getCovMatrix().GetEigenval(i)<0.0) cout<<"ERROR"<<endl;
    //    if(PCA.getCovMat().GetEigenval(i)<0.0) cout<<"ERROR"<<endl;
    //    if(Dnew.GetEigenval(i)<0.0) cout<<"ERROR"<<endl;
  }
  cout<<log_det_D_new<<" "<<log_det_C_ENM/2.<<" "<<log_det_C_PCA/2.<<" "<<2*NTOP<<endl;
  double D_B=(log_det_D_new - log_det_C_ENM/2. - log_det_C_PCA/2.)/(2*NTOP);

  double Bhatt=exp(-D_B);
  cout<<endl<<endl
      <<"#####################################"<<endl
      <<"#     The Bhatt is "<<Bhatt<<endl
      <<"#####################################"<<endl<<endl;
  cout<<"===============Done!================="<<endl;


  // +++ OPEN the OUTPUT file +++
  ofstream fout;
  fout.open(oname);
  // +++ +++ +++ +++ +++ +++ +++
  fout<<Bhatt<<endl;
  
  fout.close();

  return 0;
}
