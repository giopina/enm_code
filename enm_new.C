
//============================================================================
/* written by Giovanni Pinamonti, January 2015 */
//============================================================================

#include <iostream>
#include <string>
#include <cmath>
#include <boost/random.hpp>
#include "omp.h"
#include <vector>
#include <ctime>


#include "Structure3d.h"
#include "PrincipalComp.h"
#include "ElasticNet.h"
//#include "Matrix.h"
using namespace std;

#include "correlation.cc"

void help_display(){
  cout<<"Help:"<<endl<<
    "-f --> input file  (required) format .pdb"<<endl<<
    "-ntop --> number of modes to consider (default = all)"<<endl<<
    "-par --> parameters file (DEFAULT = PARAMS_RIBOGM.DAT)"<<endl<<
    "-beads --> beads to consider (default = C1' C2 P )"<<endl<<
    "-cutoff --> cutoff value in angstrom (default=10)"<<endl<<
    "-multicutoff --> multiple cutoff value in angstrom (default= single cutoff)"<<endl<<
    "-rprot --> cutoff value in angstrom (default=15), for protein residues"<<endl<<
    "-outbeads --> beads to consider for output(default = same as beads)"<<endl<<
    "-pca --> PCA file  (DEFAULT = none) format .pdb"<<endl<<
    "-ocorr --> output file for correlation and RWSIP (DEFAULT = none)"<<endl<<
    "-bhatt --> compute the Bhattacharya distance (DEFAULT = false)"<<endl<<
    "-ribo --> is a riboswitch? (do I need to take away the last res because it's an ADA?) (default=false)"<<endl<<
    "-dump --> dump the ENM results on file (default=false)"<<endl<< 
    "-enm --> are there pre-constructed ENM files? (DEFAULT = false)"<<endl<<
    "-name --> (DEFAULT = intput file name)"<<endl<<
    "-covfor --> compute force covariance"<<endl<<
    "-covmat --> compute covariance matrix"<<endl<<
    "-h help"<<endl;;
  return;
}


int main(int argc, char* argv[]) {
  time_t start_time=time(0);
  clock_t start_clock=clock();

  boost::mt19937 rng(time(0));
  int NTOP=-1;
  double cutoff=10;
  double rprot=15;
  char name[200];
  char file_name[200];
  char par_name[200];
  char pca_name[200];
  char corr_name[200];
  sprintf(par_name,"PARAMS_RIBOGM.DAT");
  vector<string> bead_list;  
  vector<string> outbead_list;
  vector<double> cutoff_list;
  bool corr_has_name=false;
  bool pca_has_name=false;
  bool comp_bhatt=false;
  bool ribo=false;
  bool dump=false;
  bool make_enm=true;
  bool has_name=false;
  bool compforce=false;
  bool compcov=false;
  bool par_has_name=false;
  bool ENM_ERROR=false;
  bool dump_nn=false;
  bool multi_cutoff=false;
  bool use_exp=false;
  bead_list.push_back("C1'");
  bead_list.push_back("C2");
  bead_list.push_back("P");

  cout<<" --- Program ENM_New --- "<<endl;
  cout<<"                         "<<endl;
  cout<<"                         "<<endl;
  cout<<"                         "<<endl;
  cout<<"                         "<<endl;
  cout<<"                         "<<endl;

  cout<<endl;


  int i=1;
  int nopar=0;
  while(i<argc){
    if(!strncmp(argv[i],"-f",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      sprintf(file_name,argv[i]);
      nopar+=1;
      cout<<"will compute ENM from the file "<<file_name<<endl;
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
      par_has_name=true;
    }
    else if(!strncmp(argv[i],"-beads",6)){
      bead_list.clear();
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      while(strncmp(argv[i],"[",1)!=0){cout<<"ERROR: in parameter -beads"<<endl; help_display(); return 0;}
      if(strncmp(argv[i],"[]",2)!=0){
	i=i+1;
	if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
	while(strncmp(argv[i],"]",1)!=0){
	  char bname[4];
	  if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
	  sprintf(bname,argv[i]);
	  bead_list.push_back(bname);
	  i=i+1;
	}
      }
      if(bead_list.size()==0) {
	cout<<"WARNING: no beads declared!"<<endl;;
      }
    }
    else if(!strncmp(argv[i],"-outbeads",9)){
      outbead_list.clear();
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      while(strncmp(argv[i],"[",1)!=0){cout<<"ERROR: in parameter -outbeads"<<endl; help_display(); return 0;}
      if(strncmp(argv[i],"[]",2)!=0){
	i=i+1;
	if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
	while(strncmp(argv[i],"]",1)!=0){
	  char bname[4];
	  if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
	  sprintf(bname,argv[i]);
	  outbead_list.push_back(bname);
	  i=i+1;
	}
      }
      if(outbead_list.size()==0) {
	cout<<"WARNING: no outbeads declared!"<<endl;;
      }
    }
    else if(!strncmp(argv[i],"-cutoff",7)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      if(sscanf(argv[i],"%lf",&cutoff)!=1){cout<<"ERROR: in parameter -cutoff"<<endl; help_display(); return 0;};
      cout<<"cutoff = "<<cutoff<<endl;
    }    
    else if(!strncmp(argv[i],"-multicutoff",12)){
      double min_rc, max_rc;
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      if(sscanf(argv[i],"%lf",&min_rc)!=1){cout<<"ERROR: in parameter -cutoff"<<endl; help_display(); return 0;};
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      if(sscanf(argv[i],"%lf",&max_rc)!=1){cout<<"ERROR: in parameter -cutoff"<<endl; help_display(); return 0;};
      double dum_rc=min_rc;
      while(dum_rc<=max_rc){
	cutoff_list.push_back(dum_rc);
	dum_rc+=1.;
      }
      if(cutoff_list.size()<1) cout<<"!!!ERROR IN MULTIPLE CUTOFF LIST!!!"<<endl;
      cout<<"multiple cutoff = ";
      for(int i=0;i<cutoff_list.size();i++) cout<<cutoff_list[i]<<" ";
      cout<<endl;
    }    
    else if(!strncmp(argv[i],"-rprot",6)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
      if(sscanf(argv[i],"%lf",&rprot)!=1){cout<<"ERROR: in parameter -rprot"<<endl; help_display(); return 0;};
      cout<<"rprot = "<<rprot<<endl;
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
      cout<<"\"My name is "<<name<<".\""<<endl;
      has_name=true;
    }
    else if(!strncmp(argv[i],"-covfor",7)){
      cout<<" Computing the force covariance matrix!!! "<<endl;
      compforce=true;
    }
    else if(!strncmp(argv[i],"-covmat",7)){
      cout<<" Computing the covariance matrix!!! "<<endl;
      compcov=true;
    }
    else if(!strncmp(argv[i],"-nn",3)){
      cout<<" Dump number of nearest neighbours "<<endl;
      dump_nn=true;
    }
    else if(!strncmp(argv[i],"-h",2)){
      help_display();
      return 0;
    }
    else if(!strncmp(argv[i],"-exp",4)){
      cout<<" _____| Using an exponential model |____ "<<endl;
      use_exp=true;
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

  cout<<"Beads to consider are: ";
  for(int i=0; i<bead_list.size();++i){
    cout<<bead_list.at(i)<<" ";}
  cout<<endl<<endl;

  if(outbead_list.size()>0){
    cout<<"Out Beads to consider are: ";
    for(int i=0; i<outbead_list.size();++i){
      cout<<outbead_list.at(i)<<" ";}
    cout<<endl<<endl;
  }

  if(!has_name) sprintf(name,file_name);

  ofstream corrout;
  if(corr_has_name) corrout.open(corr_name);
  /* # # # fake enm just for the reference structure! # # # */
  Structure3d ref_struc;
  vector<int> res_list;
  if(par_has_name){ //qui leggo l'eventuale file parametri (lista di OUTRES)
    char stringa[1000],line[1000];
    ifstream file_in;
    file_in.open(par_name);
    if(!file_in.is_open()){ cout<<"Problem opening parameters file: "<<par_name<<endl;}
    while(file_in.getline(line,1000)!=NULL){
      istringstream iss(line);
      if(!(iss>>stringa)) continue;
      if(!strncmp(stringa,"OUTRES:",6)){
	cout<<stringa<<" ";
	int num;
	while (iss >> num){cout<<num<<" "; res_list.push_back(num);}
	cout<<endl<<"N of res for output = "<<res_list.size()<<endl;
      }//end if OUTRES
    }//END LOOP ON FILE LINES
    file_in.close();
  }
  BeadList lista;
  if(outbead_list.size()>0)  lista.setAtomTypes(outbead_list);
  else lista.setAtomTypes(bead_list); // setto i numeri dei residui di output
  lista.setResNum(res_list);
  ref_struc.setBeadList(lista);
  cout<<ref_struc.getSize()<<endl;
  ref_struc.readFromPDBFile(file_name);
  cout<<ref_struc.getSize()<<endl;
  /* # # # fake enm just for the reference structure! # # # */

  /* //ho tolto sta parte solo momentaneamente. Quando sistemi il problema dei terminali usa questa che e' sicuro piu' pulita
  Structure3d ref_struc;//=ENM.getStructure();
  if(outbead_list.size()>0)  ref_struc.setBeadList(outbead_list);
  else ref_struc.setBeadList(bead_list);
  ref_struc.readFromPDBFile(file_name);
*/
  int n_beads=ref_struc.getSize();
  int nmodes=3*n_beads;
  PrincipalComp PCA(n_beads);
  if(pca_has_name){
    cout<<"reading pca"<<endl;
    cout<<"TIME= "<<difftime(time(0), start_time)<<endl;

    time_t pca_read_start_t=time(0);
    clock_t pca_read_start_c=clock();
    PCA.readFromFile(pca_name,nmodes);
    cout<<"             PCA READ TIME  = "<<difftime(time(0),pca_read_start_t)<<endl;
    cout<<"             PCA READ CLOCK  = "<<(double) (clock()-pca_read_start_c)/CLOCKS_PER_SEC<<endl;
  }

  if(cutoff_list.size()<1){
    cutoff_list.push_back(cutoff);
  }

  for(int i_rc=0; i_rc<cutoff_list.size();++i_rc){

    cutoff=cutoff_list.at(i_rc);

  //#####################################
  //creo il modello elastico (o lo leggo)
  ElasticNet ENM;
  ENM_ERROR=false;
  if(make_enm){
    cout<<"TIME= "<<difftime(time(0), start_time)<<endl;
    cout<<" +++ Construct ENM +++"<<endl;
    cout<<"     read parameters..."<<endl;
    if(par_has_name) ENM.readParameters(par_name);
    ENM.setBeadList(bead_list);
    if(outbead_list.size()>0) {cout<<"setto outbead"<<endl; ENM.setOutBeadList(outbead_list);}
    cout<<cutoff<<" rc"<<endl;
    ENM.setCutOff(cutoff);
    ENM.setRprot(rprot);
    if(!compcov) ENM.setFast();
    ENM.setExp(use_exp);
    cout<<"     read file..."<<endl;
    int  error=ENM.readPDBFile(file_name);
    if(error>0) {cout<<"ERROR READING THE FILE"<<endl; return 0;}
    cout<<"     contact map..."<<endl;
    ENM.constructContactMap();  
    cout<<"     interaction map..."<<endl;
    ENM.constructIntMat();
    cout<<"     solve..."<<endl;
    time_t enm_solve_start_t=time(0);
    clock_t enm_solve_start_c=clock();
    if(ENM.Solve()==0){
      cout<<"             ENM TIME  = "<<difftime(time(0),enm_solve_start_t)<<endl;
      cout<<"             ENM CLOCK  = "<<(double) (clock()-enm_solve_start_c)/CLOCKS_PER_SEC<<endl;
      if(dump_nn) ENM.dumpNearestNeighbours(name);
      if(dump){
	cout<<"     dump..."<<endl;
	ENM.setNtopVectors(NTOP);
	ENM.Dump(name);}
      if(compforce){
	ENM.computeForCovMatrix();
	ENM.dumpForCovMatrix("prova_for");}
      cout<<" ENM ready"<<endl;
    cout<<"TIME= "<<difftime(time(0), start_time)<<endl;
    }
    else ENM_ERROR=true;
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
  

  double msfcorr=-1.;
  double rmsip=-1.;
  double rwsip=-1.;
  double Bhatt=-1.;

  if(ENM_ERROR==false){
    // ##################################################################################################
    // ##################################################################################################
    // ##################################################################################################
    if(pca_has_name){
      // ##################################  Correlation ###############################################
      cout<<"compute corr"<<endl;
      cout<<"TIME= "<<difftime(time(0), start_time)<<endl;
      vector<double> MSF_pca, MSF_enm;
      for(int i=0; i<n_beads;++i){
	MSF_pca.push_back(PCA.getCovMat().GetMSF(i));
	MSF_enm.push_back(ENM.getIntMatrix().GetMSF(i));
      }
      msfcorr=correlation(MSF_pca,MSF_enm);
      cout<<endl<<endl
	  <<"#####################################"<<endl
	  <<"#     The correlation with MD's MSFs is "<<msfcorr<<endl
	  <<"#####################################"<<endl<<endl;
      cout<<"===============Done!================="<<endl;
      // ################################################################################################
      
      // ##################################  R.W.S.I.P.  ################################################
      
      //    cout<<"### RWSIP ###"<<endl;
      //      cout<<"TIME= "<<difftime(time(0), start_time)<<endl;
      time_t rwsip_start_t=time(0);
      clock_t rwsip_start_c=clock();
      //      for(int idum=0;idum<10;idum++)
      rwsip=RWSIP(PCA.getCovMat(),ENM.getCovMatrix());
      //      double rwsip2=RWSIP(PCA.getCovMat(),ENM.getIntMatrix());
      cout<<"         RWSIP TIME  = "<<difftime(time(0),rwsip_start_t)<<endl;
      cout<<"         RWSIP CLOCK  = "<<(double) (clock()-rwsip_start_c)/CLOCKS_PER_SEC<<endl;

      //      cout<<"TIME= "<<difftime(time(0), start_time)<<endl;
      rmsip=RMSIP(PCA.getCovMat(),ENM.getIntMatrix(),10);
      cout<<endl<<endl
	  <<"#####################################"<<endl
	  <<"#     The RWSIP is "<<rwsip<<endl
	  <<"#     The RMSIP is "<<rmsip<<endl
	  <<"###########################4##########"<<endl<<endl;
      cout<<"===============Done!================="<<endl;
      // ################################################################################################
      
      //++++++++++++++++++++++++++++++++++++++     Bhatt     ++++++++++++++++++++++++++++++++++++++++++++
      if(comp_bhatt){
	cout<<"+++ Bhatt man +++"<<endl;
	cout<<"TIME= "<<difftime(time(0), start_time)<<endl;
	ENM.computeCovar();
	
	Bhatt=ComputeBhatt(ENM.getCovMatrix(), PCA.getCovMat());
	
	cout<<endl<<endl
	    <<"#####################################"<<endl
	    <<"#     The Bhatt is "<<Bhatt<<endl
	    <<"#####################################"<<endl<<endl;
	cout<<"===============Done!================="<<endl;
      }
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }//endif PCAHASNAME
  }//endif ENMERROR

  if(corr_has_name){
    corrout<<cutoff<<" ";
    if(pca_has_name){
      if(rwsip<0) {msfcorr=-999;Bhatt=-999;}
      corrout<<msfcorr<<" "<<rwsip<<" "<<rmsip<<" ";
      if(comp_bhatt) corrout<<Bhatt<<" ";
    }
    corrout<<endl;
  }
  }    
  cout<<"TIME= "<<difftime(time(0), start_time)<<endl;
  if(corr_has_name){
    corrout.close();
  }
  
  return 0;
}

