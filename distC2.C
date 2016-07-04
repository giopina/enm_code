/* written by Giovanni Pinamonti, October 2014 */
//============================================================================
// originally based on
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
#include "omp.h"
#include <vector>

#include "Structure3d.h"
#include "PrincipalComp.h"
#include "ElasticNet.h"
//#include "Matrix.h"
using namespace std;

void help_display(){
  cout<<"Help:"<<endl<<
    "-f --> input file  (required) format .pdb"<<endl<<
    "-temp --> normalized temperature (default =0.1)"<<endl<<
    "-nframes --> number of frames (default = -1)"<<endl<<
    "-ntop --> number of modes to consider (default = all)"<<endl<<
    "-par --> parameters file (DEFAULT = PARAMS_RIBOGM.DAT)"<<endl<<
    "-pca --> PCA file  (DEFAULT = none) format .pdb"<<endl<<
    "-of --> output file for fluctuations (DEFAULT = none)"<<endl<<
    "-ocorr --> output file for correlations (DEFAULT = none)"<<endl<<
    "-shape --> input file with shape data [data in third column] (DEFAULT = none)"<<endl<<
    "-bhatt --> compute the Bhattacharya distance (DEFAULT = false)"<<endl<<
    "-ribo --> is a riboswitch? (do I need to take away the last res because it's an ADA?) (default=false)"<<endl<<
    "-dump --> dump the ENM results on file (default=false)"<<endl<< 
    "-enm --> are there pre-constructed ENM files? (DEFAULT = false)"<<endl<<
    "-name --> (DEFAULT = intput file name)"<<endl<<
    "-C6 --> use C6 instead of C2 (default=false)"<<endl<<
    "-h help"<<endl;;
  return;
}


int main(int argc, char* argv[]) {
  boost::mt19937 rng(time(0));
  int NTOP=-1;
  int NFRAMES=-1;
  //  double fattore=1.;
  double fattore=0.1;
  char name[200];
  char file_name[200];
  char par_name[200];
  char pca_name[200];
  char shape_name[200];
  char corr_name[200];
  char oname[1024]; 
  sprintf(par_name,"PARAMS_RIBOGM.DAT");
  vector<string> bead_list;  
  bool out_has_name=false;
  bool shape_has_name=false;
  bool corr_has_name=false;
  bool pca_has_name=false;
  bool comp_bhatt=false;
  bool ribo=false;
  bool dump=false;
  bool make_enm=true;
  bool has_name=false;
  bool use_C6=false;
  //  bead_list.push_back("C1'");
  //  bead_list.push_back("C2");
  //  bead_list.push_back("P");

  cout<<" --- Program GenConf --- "<<endl;
  cout<<" Computed the distance   "<<endl;
  cout<<" fluctuations between C2 (or C6) "<<endl;
  cout<<" creating a gaussian     "<<endl;
  cout<<" trajectory based on an  "<<endl;
  cout<<" ENM (intput)            "<<endl;

  cout<<endl;


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
    else if(!strncmp(argv[i],"-pca",4)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(pca_name,argv[i]);
      pca_has_name=true;
    }
    else if(!strncmp(argv[i],"-par",4)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(par_name,argv[i]);
      cout<<"will read parameters from the file "<<par_name<<endl;
    }
    else if(!strncmp(argv[i],"-shape",6)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(shape_name,argv[i]);
      cout<<"will read shape data from the file "<<shape_name<<endl;
      shape_has_name=true;
    }
    else if(!strncmp(argv[i],"-bhatt",6)){
      comp_bhatt=true;
      cout<<"will compute the bhatt"<<endl;
    }
    else if(!strncmp(argv[i],"-ocorr",6)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(corr_name,argv[i]);
      cout<<"will write correlations on the file "<<corr_name<<endl;
      corr_has_name=true;
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
    else if(!strncmp(argv[i],"-dump",5)){
      cout<<" $ $ $ DUMP $ $ $ "<<endl;
      dump=true;
    }
    else if(!strncmp(argv[i],"-enm",4)){
      cout<<" $ $ $ read from previous ENM $ $ $ "<<endl;
      make_enm=false;
    }
    else if(!strncmp(argv[i],"-name",5)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(name,argv[i]);
      cout<<"\"My name is "<<oname<<".\""<<endl;
      has_name=true;
    }
    else if(!strncmp(argv[i],"-C6",3)){
      cout<<" $ $ $ C6?!?! $ $ $ "<<endl;
      use_C6=true;
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

  if(!has_name) sprintf(name,file_name);

  //#####################################
  //creo il modello elastico (o lo leggo)
  ElasticNet ENM;
  if(make_enm){
    cout<<" +++ Construct ENM +++"<<endl;
    cout<<"     read parameters..."<<endl;
    ENM.readParameters(par_name);
    ENM.setFast();
    cout<<"     read file..."<<endl;
    int  error=ENM.readPDBFile(file_name);
    if(error!=0) {cout<<"ERROR READING THE FILE"<<endl;}// return 0;}
    cout<<"     contact map..."<<endl;
    ENM.constructContactMap();  
    cout<<"     interaction map..."<<endl;
    ENM.constructIntMat();
    cout<<"     solve..."<<endl;
    ENM.Solve();
    if(dump) ENM.Dump(name);
  }//endif make_enm
  else{
    cout<<" +++ Read old ENM +++"<<endl;
    ElasticNet ENM;
    cout<<"     read parameters..."<<endl;
    ENM.readParameters(par_name);
    //    ENM.setFast();
    cout<<"     read file..."<<endl;
    int  error=ENM.readPDBFile(file_name);
    if(error!=0) {cout<<"ERROR READING THE FILE"<<endl;}// return 0;}
    cout<<endl<<"$ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $"<<endl;
    cout<<"$ $ $ PARTE DEL PROGRAMMA IN PREPARAZIONE $ $ $"<<endl;
    cout<<"$ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $"<<endl<<endl;
    return 0;
  }//endelse make_enm
  //#########################################
  cout<<" +++++++++++++++++++++"<<endl;
  
  Structure3d ref_struc=ENM.getStructure();
  int n_beads=ref_struc.getSize();
  int nmodes=3*n_beads;

  //leggo quanti C2 ci sono, memorizzo le loro posizioni e cosi' via
  int* is_C2=i1t(n_beads);
  vector<int> indexC2;
  for (size_t i = 0; i < ref_struc.getSize(); ++i){
    if(ref_struc.getBead(i).getHet()) continue;
    //    if ((ref_struc.getBead(i).getAtomType()=="C2") ||
    //	(ref_struc.getBead(i).getAtomType()=="C6")
    if(use_C6){
      if (ref_struc.getBead(i).getAtomType()=="C6")
	{is_C2[i]=1; indexC2.push_back(i);}// cout<<"eccolo: "<<i<<endl;}
    }
    else{
      if (ref_struc.getBead(i).getAtomType()=="C2")
	{is_C2[i]=1; indexC2.push_back(i);}// cout<<"eccolo: "<<i<<endl;}
    }
  }
  cout<<"N. beads = "<<n_beads<<endl;
  if(indexC2.size()<2) {cout<<"ERROR: INSUFFICIENT NUMBER OF C2 BEADS!!!"<<endl; return 0;}

  int n_C2;
  if(ribo) n_C2=indexC2.size()-2;//NB: -2 per escludere ADA del riboswitch, altrimenti -1 !! //obsloleta come cosa
  else n_C2=indexC2.size()-1;

  vector<double> fluct; //vettore su cui salvero' le fluttuazioni C2-C2

  if(NFRAMES<3){ //se il numero di frame e' 2 o meno non genero traiettorie e uso la formula "analitica"
    if(NFRAMES==2) cout<<"Be careful! 2 frames doesn't make sense. I'm using the exact approximation"<<endl;
    for (size_t i = 0; i < n_C2; ++i){ //ciclo sui C2
      double new_fluc=ENM.getDistFluc(indexC2.at(i),indexC2.at(i+1)); //calcolo le fluttuazioni C2-C2 con la formula
      fluct.push_back(new_fluc); //salvo nel vector
    }//enddo i
  }//endif NFRAMES

  else{ //se il numero di frame e' maggiore di 2 allora genero la traiettoria
    Structure3d new_struc=ref_struc;
    double *dist_mean, *d_sq_mean; // definisco gli array
    dist_mean=d1t(indexC2.size()); // in cui salvero'
    d_sq_mean=d1t(indexC2.size()); // media e media quadrata
    
    boost::normal_distribution<> nd(0.0,1); // definisco la distribuzione normale da cui estrarro i numeri gaussiani
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >  rnd(rng, nd);
    
    if((NTOP<1)||(NTOP>nmodes-6)) NTOP=nmodes-6;
    vector<double> sqrt_eval;
    for (int m = 6; m < NTOP+6; ++m) {
      double s=ENM.CovGetEigenval(nmodes-m-1);
      if(s<0.0){cout<<"genConf *** WARNING: negative eigenvalue -> "<<s<<" on mode: "<<m<<endl; s=abs(s); }//check
      if(s<0.000001){cout<<"genConf *** WARNING: null eigenvalue on mode: "<<m<<endl; continue; } //non fanno male
      //      s=sqrt(s)/5.0*fattore;
      s=sqrt(s)*fattore;
      sqrt_eval.push_back(s);
    }
    vector<Vector3d> coords(indexC2.size());
    
    int NFPRINT=NFRAMES/10;
    if(NFPRINT<1) NFPRINT=1;
    for (size_t f = 0; f < NFRAMES; ++f){
      if(f%NFPRINT==0) cout<<f<<" frames generated"<<endl;
      for (size_t i = 0; i < indexC2.size(); ++i){
	int i_bead=indexC2.at(i);
	coords[i]=ref_struc.getBead(i_bead).getCoordinates();
      }
      
      vector<double> width;
      for (int m = 6; m < NTOP+6; ++m) {
	//double s=ENM.CovGetEigenval(nmodes-m-1);       
	double s=sqrt_eval[m-6];
	width.push_back(rnd()*s);
      }
      //  $  $  $  $  $  $  $  $
#pragma omp parallel for
      for (size_t i = 0; i < indexC2.size(); ++i){
	int i_bead=indexC2.at(i);
	for (int m = 6; m < NTOP+6; ++m) {
	  coords[i] += ENM.CovGetEigenvec(nmodes-m-1)[i_bead]*width[m-6]; //sbagliavo qui...
	}
      }
      //  $  $  $  $  $  $  $  $
      
      // #  #  #  #  #  #  # 
      for (size_t i = 0; i < indexC2.size(); ++i){ //### qui calcolo le distanze
	double temp_dsq=coords[i].distanceSQ(coords[i+1]);
	d_sq_mean[i]+=temp_dsq;
	dist_mean[i]+=sqrt(temp_dsq);
      }

    }//enddo frames    
   
    for (size_t i = 0; i < n_C2; ++i){ 
      //  for (size_t i = 0; i < indexC2.size()-2; ++i){ //NB: -2 per escludere ADA del riboswitch, altrimenti -1 !!
      d_sq_mean[i]/=NFRAMES;
      dist_mean[i]/=NFRAMES;
      fluct.push_back(d_sq_mean[i]-(dist_mean[i]*dist_mean[i]));
    }
  }//endelse

  /// +++ OPEN the OUTPUT file +++
  if(out_has_name){
    ofstream fout;
    fout.open(oname);
    //    ENM.setFast(false);
    //    ENM.computeCovar();

    //    ENM.Dump("prova");
    //    ENM.dump_top_modes(10,5,"prova");
    // ### qui calcolo le distanze e le printo ###
    //    ENM.dumpDistFluc("prova");
    for (size_t i = 0; i < n_C2; ++i){ 
      //    for (size_t i = 0; i < indexC2.size()-2; ++i){ //NB: -2 per escludere ADA del riboswitch, altrimenti -1 !!
      fout<<i<<" "<<fluct.at(i)<<" "<<indexC2.at(i)<<" "<<ref_struc.getBead(indexC2.at(i)).getResName()<<" "<<ref_struc.getBead(indexC2.at(i)).getResNum()<<endl;
    }
    fout.close();
  }// +++ +++ +++ +++ +++ +++ +++
  double kend=0.0;
  double corr=0.0;
  if(shape_has_name){
    vector<double> shape;
    ifstream shapefile(shape_name);
    char cdum[2];
    int idum;
    double ddum;
    double tmpshape;
    ///    while(shapefile >> cdum >> idum >> tmpshape >> cdum >> idum >>ddum >> cdum >> idum >>ddum >> cdum >> idum >>ddum >> cdum >> idum >>ddum){
    string line;
    while(getline(shapefile,line)){
      std::istringstream iss(line);
      iss >> cdum >> idum >> tmpshape;
      if(!strncmp(cdum,"#",1)) continue;
      shape.push_back(tmpshape);
    }
    cout<<"### corr ###"<<endl;
    int ndata=shape.size();
    if(!(ndata==fluct.size())){cout<<"ERROR: lunghezza dati shape diversa da fluttuazioni. Non posso correlarle. "<<ndata<<" "<<fluct.size()<<endl;}
    else{
      double C=0;
      double D=0;
      double C_n_2=ndata*(ndata-1)/2.;
      
      double sum_X=0.0;
      double sum_Y=0.0;
      double sum_X2=0.0;
      double sum_Y2=0.0;
      double sum_XY=0.0;
      
      //#pragma omp parallel for
      for(int i=0;i<ndata;++i){
	sum_X+=fluct[i];
	sum_X2+=fluct[i]*fluct[i];
	sum_Y+=shape[i];
	sum_Y2+=shape[i]*shape[i];
	sum_XY+=fluct[i]*shape[i];
	//#pragma omp parallel for
	for(int j=i+1;j<ndata;++j){
	  if((fluct[i]-fluct[j])*(shape[i]-shape[j])>0.0)
	    C+=1;
	  else
	    if((fluct[i]-fluct[j])*(shape[i]-shape[j])<0.0)
	      D+=1;
	}//enddo j
      }//enddo i
      kend=(C-D)/C_n_2;
      corr=(sum_XY-(sum_X*sum_Y)/ndata)/(sqrt((sum_X2-(sum_X*sum_X)/ndata)*(sum_Y2-(sum_Y*sum_Y)/ndata)));
      cout<<endl<<endl
	  <<"#####################################"<<endl
	  <<"#     The correlations are "<<corr<<" "<<kend<<endl
	  <<"#####################################"<<endl<<endl;
      cout<<"===============Done!================="<<endl;
    }//endelse
  }//endif shape_has_name

  // ##################################################################################################
  // ##################################  R.W.S.I.P.  ##################################################
  // ##################################################################################################
  double RWSIP=0.0;
  double Bhatt=0.0;
  if(pca_has_name){
    PrincipalComp PCA(n_beads);
    PCA.readFromFile(pca_name,nmodes);
    cout<<"### RWSIP ###"<<endl;
    RWSIP=0.0;
    double NORM=0.0;
    //#pragma omp parallel for
    for (int m = 6; m < nmodes; ++m) {
      double eval_ENM=ENM.CovGetEigenval(m-6);
      //      cout<<"ENM: "<<eval_ENM;
      vector<Vector3d> evec_ENM=ENM.CovGetEigenvec(m-6);
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
      //      cout<<" PCA: "<<eval_PCA_i<<endl;
      NORM+=eval_PCA_i*eval_ENM;
    }
    RWSIP/=NORM;
    RWSIP=sqrt(RWSIP);
    cout<<endl<<endl
	<<"#####################################"<<endl
	<<"#     The RWSIP is "<<RWSIP<<endl
	<<"#####################################"<<endl<<endl;
    cout<<"===============Done!================="<<endl;
    //}
  // ##################################################################################################
  // ##################################################################################################
  // ##################################################################################################




  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++     Bhatt     ++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //  if((pca_has_name)&&(comp_bhatt)){
    if(comp_bhatt){
      cout<<"+++ Bhatt man +++"<<endl;
      ENM.computeCovar();
      double norm_C_ENM=ENM.getCovMatrix().GetTrace();
      double norm_C_PCA=PCA.getCovMat().GetTrace();
      ENM.getCovMatrix().Times(1./norm_C_ENM);
      PCA.getCovMat().Times(1./norm_C_PCA);
      cout<<"creo D"<<endl;
      CovMatrix D=PCA.getCovMat();
      cout<<"creata D"<<endl;
      D.Plus(ENM.getCovMatrix());
      D.Times(0.5);
      D.Decompose(); // NB: tutto e' ordinato con gli eval crescenti! (ma perche'??!)
      
      // calcolo il numero di autovalori per avere il *nu* per cento di RMSD
      double nu_frac=0.95; //fraction of RMSD to be reproduced by the eigenvalues
      double somma=0.0;
      for(int i=0; i<nmodes; i++){somma+=D.GetEigenval(nmodes-i-1); }
      double partsomma=0.0;
      int NTOP=0;
      for(int i=0;i<nmodes;i++){
	NTOP++;
	partsomma+=D.GetEigenval(nmodes-i-1);
	if(partsomma>nu_frac*somma) break;
      }
      cout<<"NTOP="<<NTOP<<endl;
      NTOP/=3;  NTOP*=3;
      cout<<"NTOP corretto="<<NTOP<<endl;
      
      cout<<"Reducing"<<endl;
      ENM.getCovMatrix().Reduce(NTOP,D);
      PCA.getCovMat().Reduce(NTOP,D);
      ENM.getCovMatrix().Decompose();
      PCA.getCovMat().Decompose();
      CovMatrix Dnew=PCA.getCovMat();
      Dnew.Plus(ENM.getCovMatrix());
      Dnew.Times(0.5);
      Dnew.Decompose();
      
      double log_det_C_ENM=0.0;
      double log_det_C_PCA=0.0;
      double log_det_D_new=0.0;
      for(int i=0;i<NTOP;++i){
	//    cout<<<<endl;
	//	cout<<ENM.getCovMatrix().GetEigenval(NTOP-i-1)<<" "<<
	//	  PCA.getCovMat().GetEigenval(NTOP-i-1)<<" "<<
	//	  Dnew.GetEigenval(NTOP-i-1)<<endl;
	
	log_det_C_ENM+=log(ENM.getCovMatrix().GetEigenval(NTOP-i-1));
	log_det_C_PCA+=log(PCA.getCovMat().GetEigenval(NTOP-i-1));
	log_det_D_new+=log(Dnew.GetEigenval(NTOP-i-1));
      }
      //      cout<<log_det_D_new<<" "<<log_det_C_ENM/2.<<" "<<log_det_C_PCA/2.<<" "<<2*NTOP<<endl;
      double D_B=(log_det_D_new - log_det_C_ENM/2. - log_det_C_PCA/2.)/(2*NTOP);
      
      Bhatt=exp(-D_B);
      cout<<endl<<endl
	  <<"#####################################"<<endl
	  <<"#     The Bhatt is "<<Bhatt<<endl
	  <<"#####################################"<<endl<<endl;
      cout<<"===============Done!================="<<endl;
    }
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  if(corr_has_name){
    ofstream corrout;
    corrout.open(corr_name);
    if(shape_has_name){
      corrout<<corr<<" "<<kend<<" ";
    }
    if(pca_has_name){
      corrout<<RWSIP<<" "<<Bhatt<<" ";
    }
    corrout<<endl;
    corrout.close();
  }

  return 0;
}

