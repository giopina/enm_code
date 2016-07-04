/* Created by Giovanni Pinamonti */
/* PhD student @ SISSA, Trieste  */
/* November 20th, 2013           */

#include <time.h>
#include "omp.h"
#include "ElasticNet.h"
#include "io.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "lapack_matrix_routines_wrapper_v3.cc"
#include "BeadList.h"

#define BROKEN_CHAIN_CUTOFF 7.0

#define DEFAULT_CUTOFF 10.0
#define DEFAULT_NTOP 10

//###### COSTRUTTORI DELLA CLASSE ######
ElasticNet::ElasticNet():
  //  _covMatrix(0)
  _N(0),
  _ntop_vectors(DEFAULT_NTOP),
  _cutoff(DEFAULT_CUTOFF),

  _R_BB(-1),
  _R_B(-1),
  _R_PP(-1),
  _R_prot(-1),
  _R_N7(-1),

  _TRACCIAMENTO(false),
  _Kcov(1.0),
  _KBB(1.0),
  _KB(1.0),
  _KPP(1.0),
  _COM(false),
  _OUTRES(false),
  //  _PROTEIN(false),
  _FAST(false),
  _DPF(false),
  _RANDOM(false),
  _RARAND(false),
  _EXP(false),
  _DUMPMODES(false),
  _DUMPEIGENVALUES(true),
  _DUMPCOVMAT(false),
  _DUMPDISTFLUC(false),
  _DUMPREDCOVMAT(false),
  _DUMPMSD(true)
{
  _COM=false;
}

ElasticNet::ElasticNet(const char *filename):
  _ntop_vectors(DEFAULT_NTOP),
  _cutoff(DEFAULT_CUTOFF),

  _R_BB(-1),
  _R_B(-1),
  _R_PP(-1),
  _R_prot(-1),
  _R_N7(-1),

  _TRACCIAMENTO(false),
  _Kcov(1.0),
  _KBB(1.0),
  _KB(1.0),
  _KPP(1.0),
  _COM(false),
  _OUTRES(false),
  //  _PROTEIN(false),
  _FAST(false),
  _DPF(false),
  _RANDOM(false),
  _RARAND(false),
  _EXP(false),
  _DUMPMODES(false),
  _DUMPEIGENVALUES(true),
  _DUMPCOVMAT(false),
  _DUMPDISTFLUC(false),
  _DUMPREDCOVMAT(false),
  _DUMPMSD(true)
{
  //  cout<<"ElasticNet::reading file "<<filename<<endl;
  readPDBFile(filename);
}



ElasticNet::ElasticNet(Structure3d structure):
  _ntop_vectors(DEFAULT_NTOP),
  _cutoff(DEFAULT_CUTOFF),

  _R_BB(-1),
  _R_B(-1),
  _R_PP(-1),
  _R_prot(-1),
  _R_N7(-1),

  _TRACCIAMENTO(false),
  _Kcov(1.0),
  _KBB(1.0),
  _KB(1.0),
  _KPP(1.0),
  _COM(false),
  _OUTRES(false),
  //  _PROTEIN(false),
  _FAST(false),
  _DPF(false),
  _RANDOM(false),
  _RARAND(false),
  _EXP(false),
  _DUMPMODES(false),
  _DUMPEIGENVALUES(true),
  _DUMPCOVMAT(false),
  _DUMPDISTFLUC(false),
  _DUMPREDCOVMAT(false),
  _DUMPMSD(true)
{
  cout<<"ElasticNet: setting the structure"<<endl;
  setStructure(structure);
  // cout<<"ElasticNet: creating the contact map"<<endl;
  // constructContactMap();
  // cout<<"ElasticNet: creating the interaction matrix"<<endl;
  // constructIntMat();
}
#undef DEFAULT_CUTOFF
#undef DEFAULT_NTOP

ElasticNet::~ElasticNet() {
  cout<<"deleting"<<endl;
  //Be careful: Cmap, in principle, could be not allocated.
  free_i2t(_cmap);
}


// ### ### METODI DI ELASTICNET ### ###
//set the structure
void ElasticNet::setStructure(Structure3d structure){
  _structure=structure;
  _N=_structure.getSize();
  //  _size=3*_N;
  cout<<"ElasticNet:: costruita struttura di dimensione " <<_N<<endl;
}

int ElasticNet::readPDBFile(const char* fname){
  FILE* fp;
  if ((fp = fopen (fname, "r")) == NULL) {
    fprintf (stderr,"Could not open file %s.\n.Exiting.\n", fname);
    exit (1);    }
  int idum= readPDBFile(fp);
  fclose(fp);
  return idum;
}

//read a PDB file
int ElasticNet::readPDBFile(FILE *fp){
  int error;
  if(_COM) {error=_structure.centersFromPDBFile(fp);cout<<"***********8COM"<<endl;}
  else error=_structure.readFromPDBFile(fp);
  _N=_structure.getSize();
  //  _size=3*N;
  cout<<"ElasticNet:: costruita struttura di dimensione "<<_N<<endl;
  return error;
};


/* decompose the interaction matrix */
int ElasticNet::Solve(){
  //_intMatrix->Decompose()
  if(_TRACCIAMENTO) {
    if(tracciamento()<0){cout<<"ERROR IN TRACCIAMENTO"<<endl; return -1;} 
  }
  
  _intMatrix.Decompose();

  time_t begin_time2=time(0);
  if (_FAST) {cout<<"ElasticNet::I don't compute the covariance matrix to save time!!!"<<endl; return 0;}
  else computeCovar();
  cout << "Solve: tempo nel ciclo 2 = "<<difftime(time(0),begin_time2)<<endl;    
  //  IntMatrix mat_prova=_intMatrix;  
  //  cout<<"ElasticNet:: Computing the covariance matrix..."<<endl;
  //  cout<<"ElasticNet:: invert..."<<endl;
  //  _covMatrix=_intMatrix.GetInverse();
  //  cout<<"ElasticNet:: decompose..."<<endl;
  //  _covMatrix.Decompose();

  return 0; //success
}

void ElasticNet::computeCovar(){
  cout<<"ElasticNet:: Computing the covariance matrix..."<<endl;
  cout<<"ElasticNet:: invert..."<<endl;
  _covMatrix=_intMatrix.GetInverse();
  //cout<<"ElasticNet:: CovMatrix -> decompose..."<<endl;
  //  _covMatrix.Decompose();
  return;
};


/* dump the important stuff of the covariance matrix
 AND the mean square displacements */
void ElasticNet::Dump(const char *name){
  //qui puoi dumpare quello che vuoi. Non e' male in realta'


  string fname(name);

  cout<<"dumpmodes "<<_DUMPMODES<<endl;  
  if(_FAST) {
    cout<<"          FASTASHELL!!"<<endl;
    _intMatrix.dumpCovEigenvalues(fname);
    _intMatrix.dumpTopCovVectors(_ntop_vectors,fname);
    cout<<"          FASTASHELL!!"<<endl;
  }
  else {    
    cout<<"ElasticNet::Dump -> _covMatrix.dumpEigenvalues()"<<endl;
    _covMatrix.dumpEigenvalues(name);
    cout<<"ElasticNet::Dump -> _covMatrix.dumpEigenvectors()"<<endl;
    cout<<"vec="<<_ntop_vectors<<endl;
    _covMatrix.dumpTopVectors(_ntop_vectors,name);
    
    if(_DUMPCOVMAT)    _covMatrix.dumpFullMatrix(fname+"_cov");
    if(_DUMPREDCOVMAT) {_covMatrix.dumpReducedMatrix(fname+"_cov");
      _covMatrix.dumpNormalizedReducedMatrix(fname+"_cov");}
    _intMatrix.dumpReducedMatrix(fname+"_int");
    
    if(_DUMPDISTFLUC) dumpDistFluc(fname);

    cout<<"ElasticNet::Dump -> Structure"<<endl;
    _structure.dumpPDBFile(fname+"_structure.pdb");
  }

  _intMatrix.dumpFullMatrix(fname);
  cout<<"ElasticNet::Dump -> dumpMSF()()"<<endl;
  //  _intMatrix.dumpMSF(name);
  dumpMSF(name);

  if(_DPF){
    char filename[200];
    ofstream fp;
    sprintf(filename,"%s_mean_square_fluc_part.dat",name);
    fp.open(filename);
    fp<<"# n. bead & atom type & n. res & xx & xy & xz & yy & yz & zz"<<endl;
    for(int i=0; i < _intMatrix.GetN(); i++){
      fp <<i <<" "
	 <<_structure.getBead(i).getAtomType()<<" "
	 <<_structure.getBead(i).getResNum()<<" ";
	for(int mu=0;mu<3;mu++){
	for(int nu=mu;nu<3;nu++){
	  //	  string coord1,coord2;
	  // if(mu==0)  coord1="x";
	  // if(mu==1)  coord1="y";
	  // if(mu==2)  coord1="z";
	  // if(nu==0)  coord2="x";
	  // if(nu==1)  coord2="y";
	  // if(nu==2)  coord2="z";
	  double temp;
	  if(_FAST) temp=_intMatrix.GetMSF(i,mu,nu);
	  else temp=_covMatrix.GetMSF(i,mu,nu);
	  fp<<temp<<" ";
	    }}
      fp<<endl;
    }
    fp.close();
    cout<<"Beads' mean square partial fluctuations written to file "<<filename<<endl;
    
  }//endif_DPF

  if(_DUMPMODES){
    cout<<"dumping modes"<<endl;
    if(_FAST){ cout<<"ElasticNet -> PROBLEM: won't dump principal modes in FAST configuration!"<<endl; return; }
    //    if(_ntop_vectors<11) dump_top_modes(_ntop_vectors,5,name);
    //    else dump_top_modes(10,5,name);
    dump_top_modes(_ntop_vectors,5,name);
  }

  cout<<"dumping the number of NN"<<endl;
  dumpNearestNeighbours(name);

}


void ElasticNet::constructContactMap(){
  //TODO una cosa del genere?
  // 1 -> contatto normale
  // 3 -> legame covalente
  // 10 -> contatto proteina-proteina
  // 11 -> contatto proteina-RNA

  //  double R_BB=5.0;

  if(_R_B<0)    _R_B   =_cutoff;
  if(_R_PP<0)   _R_PP  =_cutoff;
  if(_R_BB<0)   _R_BB  =_cutoff;
  if(_R_prot<0) _R_prot=_cutoff;
  if(_R_N7<0)   _R_N7  =_cutoff;

  cout<<" *** RADIUS OF CUT-OFF *** "<<endl;
  cout<<" *** Rc  ="<<_cutoff<<endl;
  cout<<" *** R_B ="<<_R_B   <<endl;
  cout<<" *** R_PP="<<_R_PP  <<endl;
  cout<<" *** R_BB="<<_R_BB  <<endl;
  cout<<" *** R_N7="<<_R_N7  <<endl;
  cout<<" ************************* "<<endl;

  _cmap=i2t(_N,_N);

  for(int i=0; i < _N-1; i++){
    for(int j= i+1; j < _N; j++){      
      /* Get the distance of beads i and j */
      double d = _structure.getDistance(i,j);
      /* Check if they are closer than Int_Range */
      Bead bi=_structure.getBead(i);
      Bead bj=_structure.getBead(j);
      double Rc=_cutoff;
      int kind_o_int=1;
      //22 -> base-base
      //20 -> base other
      //99 -> P-P
      //77 -> prot-prot
      //77 -> prot-RNA
      

      if (bi.getAtomType()=="N7"){
	//	if (bj.getAtomType()=="N7"){
	//	  Rc=_R_BB;  kind_o_int=22;}
	//	else
	{Rc=_R_N7;  kind_o_int=20;}
      }
      if (bj.getAtomType()=="N7"){
	//	if (bi.getAtomType()=="N7"){
	//	  Rc=_R_BB;  kind_o_int=22;}
	//	else
	{Rc=_R_N7;  kind_o_int=20;}
      }


      // if ((bi.getAtomType()=="C2")||(bi.getAtomType()=="N7")){
      // 	if ((bj.getAtomType()=="C2")||(bj.getAtomType()=="N7")){
      if (bi.getAtomType()=="C2"){
	if (bj.getAtomType()=="C2"){
	  Rc=_R_BB;  kind_o_int=22;}
	else
	  {Rc=_R_B;  kind_o_int=20;}
      }
      // if ((bj.getAtomType()=="C2")||(bj.getAtomType()=="N7")){
      // 	if ((bi.getAtomType()=="C2")||(bi.getAtomType()=="N7")){
      if (bj.getAtomType()=="C2"){
	if (bi.getAtomType()=="C2"){
	  Rc=_R_BB;  kind_o_int=22;}
	else
	  {Rc=_R_B;  kind_o_int=20;}
      }


      
      if (bi.getAtomType()=="P")
	if (bj.getAtomType()=="P"){
	  Rc=_R_PP; kind_o_int=99;}
      
      double Rhybrid=(_cutoff+_R_prot)*0.5;
      //protein RNA interaction
      if (bi.getAtomType()=="CA"){
	if (bj.getAtomType()=="CA"){
	  Rc=_R_prot;kind_o_int=77;}
	else{ Rc=Rhybrid;kind_o_int=70;}
      }
      else{
	if (bj.getAtomType()=="CA"){
	  Rc=Rhybrid;kind_o_int=70;}
      }
      
      
      if (d < Rc) _cmap[i][j]=_cmap[j][i]=kind_o_int;
    }//end do j
  }//end do i


  for(int i=0; i < _N-1; i++){
    /* ### stiffest springs between bonded nodes ### */
    Bead bi=_structure.getBead(i);
    if (bi.getAtomType()=="P"){
      int j=i+1;
      Bead bj=_structure.getBead(j);
      if (bj.getAtomType()=="C1'"){
	double d=_structure.getDistance(i,j);
	if (d < BROKEN_CHAIN_CUTOFF)
	  _cmap[i][j]=_cmap[j][i]=3;
      }//endif
    }//endif P

    else if (bi.getAtomType()=="C1'"){
      int j=i+1;
      Bead bj=_structure.getBead(j);
      if (bj.getAtomType()=="C2"){
       	double d=_structure.getDistance(i,j);
       	if (d < BROKEN_CHAIN_CUTOFF)
       	  _cmap[i][j]=_cmap[j][i]=3;
      }//endif
       if (bj.getAtomType()=="N7"){
        	double d=_structure.getDistance(i,j);
        	if (d < BROKEN_CHAIN_CUTOFF)
        	  _cmap[i][j]=_cmap[j][i]=3;
       }//endif

      j=i+2;
      if(j>_N-1) continue;
      bj=_structure.getBead(j);
      if (bj.getAtomType()=="P"){
	double d=_structure.getDistance(i,j);
	if (d < BROKEN_CHAIN_CUTOFF)
	  _cmap[i][j]=_cmap[j][i]=3;
      }//endif

       j=i+3;
       if(j>_N-1) continue;
       bj=_structure.getBead(j);
       if (bj.getAtomType()=="P"){
       	double d=_structure.getDistance(i,j);
       	if (d < BROKEN_CHAIN_CUTOFF)
       	  _cmap[i][j]=_cmap[j][i]=3;
       }//endif

    }//endif C1'

    else if (bi.getAtomType()=="N7"){
      int j=i+1;
      Bead bj=_structure.getBead(j);
      if (bj.getAtomType()=="C2"){
	double d=_structure.getDistance(i,j);
	if (d < BROKEN_CHAIN_CUTOFF)
	  _cmap[i][j]=_cmap[j][i]=3;
      }//endif
    }//endif N7

  }//enddo I

}//end function


void ElasticNet::constructRandMat(double **intmat){
  cout<<"ElasticNet: creating a random interaction matrix"<<endl;
  for(int i=0; i < _N; i++){
    for(int j=0; j < _N; j++){
      if (i==j) continue;
      double k_int=((double) rand() / ((double) RAND_MAX));
      double d = _structure.getDistance(i,j);
      /* Let's build the distance vector */
      Vector3d dvec= _structure.getDistanceVector(i,j);      
      for (int mu=0; mu < 3; mu++){
	for (int nu=0; nu < 3; nu++){
	  double temp = 0.5*k_int*dvec[mu]*dvec[nu]/(d*d);
	  intmat[3*i+mu][3*i+nu] += temp;
	  intmat[3*i+mu][3*j+nu] += -temp;
	  intmat[3*j+mu][3*i+nu] += -temp;
	  intmat[3*j+mu][3*j+nu] += temp;
	}//end for nu
      }// end for mu      
    }//i
  }//j
  return;
}

void ElasticNet::constructDemiRandMat(){
  cout<<"ElasticNet: shuffling the connectivity matrix to obtain a demi-random interaction matrix"<<endl;
  
  //  int n_conn=0;

  vector <Connection> connessioni;
  for(int i=0;i<_N;++i){
    for (int j=i+1;j<_N;++j){
      if (i==j) continue;
      if(_cmap[i][j]!=0){
	//n_conn++;
      if(_cmap[i][j]>1) cout<<i<<" "<<j<<endl;
	Connection temp_conn(i,j);
	connessioni.push_back(temp_conn);
      }
    }
  }

  int N_CONN=connessioni.size();
  cout<<"ci sono "<<N_CONN<<" connessioni"<<endl;

  int N_STEP_MAX=100*N_CONN;

  int n2swap=2*N_CONN;

  int n_swapped=0;
  for(int n=0;n<N_STEP_MAX;++n){
    //    cout<<n<<endl;  
    int i= rand()%N_CONN;
    int j= rand()%N_CONN;

    Connection tmp_connA=connessioni.at(i);
    Connection tmp_connB=connessioni.at(j);

    int old_i1_A=tmp_connA.GetI1();
    int old_i2_A=tmp_connA.GetI2();
    int old_i1_B=tmp_connB.GetI1();
    int old_i2_B=tmp_connB.GetI2();

    cout<<"Scelgo queste due connessioni: "<<endl
	<<"["<<i<<"] "<<tmp_connA.GetI1()<<"-"<<tmp_connA.GetI2()<<endl
	<<"["<<j<<"] "<<tmp_connB.GetI1()<<"-"<<tmp_connB.GetI2()<<endl;    
    
    if(!(RandSwap(tmp_connA, tmp_connB)) ) {cout<<"Non le scambio!"<<endl; continue;}

    int i1_A=tmp_connA.GetI1();
    int i2_A=tmp_connA.GetI2();
    int i1_B=tmp_connB.GetI1();
    int i2_B=tmp_connB.GetI2();
    cout<<"ora connetterebbero i seguenti nodi: "<<endl
	<<"["<<i<<"] "<<i1_A<<"-"<<i2_A<<endl
	<<"["<<j<<"] "<<i1_B<<"-"<<i2_B<<endl;

    if(_cmap[i1_A][i2_A]>0){cout<<"[i] were already connected"<<endl; continue;}
    if(_cmap[i1_B][i2_B]>0){cout<<"[j] were already connected"<<endl; continue;}
    
    _cmap[i1_A][i2_A]++;
    _cmap[i1_B][i2_B]++;
    _cmap[i2_A][i1_A]++;
    _cmap[i2_B][i1_B]++;

    _cmap[old_i1_A][old_i2_A]=0;
    _cmap[old_i1_B][old_i2_B]=0;
    _cmap[old_i2_A][old_i1_A]=0;
    _cmap[old_i2_B][old_i1_B]=0;

    connessioni.at(i)=tmp_connA;
    connessioni.at(j)=tmp_connB;

    n_swapped++;
    cout<<"SCAMBIATE!"<<endl;

    if(n_swapped>n2swap) break;
  }//end for
  cout<<"Ho scambiato "<<n_swapped<<" connessioni"<<endl;

  int n_conn=0;
  for(int i=0;i<_N;++i){
    for (int j=i+1;j<_N;++j){
      if (i==j) continue;
      if(_cmap[i][j]>1) cout<<i<<" "<<j<<endl;
      if(_cmap[i][j]==1){
	n_conn++;
      }
    }
  }
  cout<<"ci sono "<<n_conn<<" connessioni"<<endl;
  return;
}

void ElasticNet::constructExpMat(double **intmat){
  cout<<"ElasticNet: creating an exponential interaction matrix"<<endl;
  for(int i=0; i < _N; i++){
    for(int j=0; j < _N; j++){
      if (i==j) continue;
      double k_int=1;
      double d = _structure.getDistance(i,j);
      /* Let's build the distance vector */
      Vector3d dvec= _structure.getDistanceVector(i,j);
      k_int *= exp(-((d/_cutoff)*(d/_cutoff)) ) ;    
      for (int mu=0; mu < 3; mu++){
	for (int nu=0; nu < 3; nu++){
	  double temp = 0.5*k_int*dvec[mu]*dvec[nu]/(d*d);
	  intmat[3*i+mu][3*i+nu] += temp;
	  intmat[3*i+mu][3*j+nu] += -temp;
	  intmat[3*j+mu][3*i+nu] += -temp;
	  intmat[3*j+mu][3*j+nu] += temp;
	}//end for nu
      }// end for mu
      
    }//i
  }//j
  return;
}

void ElasticNet::constructIntMat(){
  double **intmat=d2t(3*_N,3*_N);
  /* Include contact interactions i-j*/

  if(_RANDOM)    constructRandMat(intmat);
  else if(_EXP) constructExpMat(intmat);
  else {
    
    if(_RARAND) constructDemiRandMat();
    cout<<"ElasticNet: start computing the matrix"<<endl;
    for(int i=0; i < _N; i++){
      for(int j=0; j < _N; j++){
	double k_int=1;
	if (i==j) continue;
	if (_cmap[i][j]==0) continue;
	double d = _structure.getDistance(i,j);
	/* Let's build the distance vector */
	Vector3d dvec= _structure.getDistanceVector(i,j);
	/* Are you covalent? */
	if(_cmap[i][j]==3) {k_int*=_Kcov;}
	else if(_cmap[i][j]==22) k_int*=_KBB;
	else if(_cmap[i][j]==20) k_int*=_KB;
	else if(_cmap[i][j]==99) k_int*=_KPP;
	
	for (int mu=0; mu < 3; mu++){
	  for (int nu=0; nu < 3; nu++){
	    double temp = 0.5*k_int*dvec[mu]*dvec[nu]/(d*d);
	    intmat[3*i+mu][3*i+nu] += temp;
	    intmat[3*i+mu][3*j+nu] += -temp;
	    intmat[3*j+mu][3*i+nu] += -temp;
	    intmat[3*j+mu][3*j+nu] += temp;
	  }//end for nu
	}// end for mu
      }//end for j
    }//end for i
  }//endif
  cout<<"ElasticNet: saving it as a 'Matrix'"<<endl;
  _intMatrix.SetMatrix(_N,intmat);
  free_d2t(intmat);
  printf("End construction of matrix with effective quadratic interactions among beads's\n");
  return;
}


void ElasticNet::readParameters(const char* filename){
  cout<<"ElasticNet::redParameters"<<endl;
  //  FILE *fp;
  char stringa[1000];
  char argument[1000];
  char line[1000];
  char name[200];
  ifstream file_in;

  sprintf(name,"%s",filename);
  //  fp = open_file_r(name);                      
  file_in.open(name);
  if(!file_in.is_open()){
    cout<<"ElasticNet::readParameters-> Problem opening file: "<<name<<endl;
    return;
  }
  
  //  while(fgets(line,1000,fp)!=NULL){
  while(file_in.getline(line,1000)!=NULL){

    istringstream iss(line);

    //    sscanf(line,"%s %s", string, argument);
    if(!(iss>>stringa)) continue;

    if(!strncmp(stringa,"Int_Range",9)){iss>>argument; sscanf(argument,"%lf",&_cutoff);}
    if(!strncmp(stringa,"R_BB",4)){iss>>argument; sscanf(argument,"%lf",&_R_BB);}
    if(!strncmp(stringa,"R_Ba",4)){iss>>argument; sscanf(argument,"%lf",&_R_B);}
    if(!strncmp(stringa,"R_PP",4)){iss>>argument; sscanf(argument,"%lf",&_R_PP);}
    if(!strncmp(stringa,"R_prot",6)){iss>>argument; sscanf(argument,"%lf",&_R_prot);}
    if(!strncmp(stringa,"R_N7",4)){iss>>argument; sscanf(argument,"%lf",&_R_N7);}


    if(!strncmp(stringa,"N_SLOWEST_EIGENVECT",19)){iss>>argument; sscanf(argument,"%d",&_ntop_vectors);}


    // if(!strncmp(stringa,"Smooth_Cutoff",13)) sscanf(argument,"%d",&Smooth_Cutoff);
    // if(!strncmp(stringa,"nodesXnucl",10)) sscanf(argument,"%d",&nodesXnucl);
    if(!strncmp(stringa,"K_COV",5)) {iss>>argument; sscanf(argument,"%lf",&(_Kcov)); cout<<"K_COV="<<_Kcov<<endl;}
    if(!strncmp(stringa,"K_BB",4)) {iss>>argument; sscanf(argument,"%lf",&(_KBB)); cout<<"K_BB="<<_KBB<<endl;}
    if(!strncmp(stringa,"K_PP",4)) {iss>>argument; sscanf(argument,"%lf",&(_KPP)); cout<<"K_PP="<<_KPP<<endl;}
    if(!strncmp(stringa,"K_Ba",4)) {iss>>argument; sscanf(argument,"%lf",&(_KB)); cout<<"K_Ba="<<_KB<<endl;}
    // if(!strncmp(stringa,"SPECIAL_K",9)) sscanf(argument,"%d",&SPECIAL_K);
    // if(!strncmp(stringa,"BASE_CUTOFF",11)) sscanf(argument,"%lf",&BASE_CUTOFF);
    // if(!strncmp(stringa,"RangeSmooth",11)) sscanf(argument,"%lf",&RangeSmooth);
    // if(!strncmp(stringa,"COMPUTE_DIST_MATRIX",19)) sscanf(argument,"%d",&COMPUTE_DIST_MATRIX);

    if(!strncmp(stringa,"BEADS:",6)){
      cout<<"ElasticNet::   "<<stringa<<" ";
      char name[4];
      _beadList.clear();
      while (iss >> name)
	{
	  cout<<name<<" ";
	  _beadList.push_back(name);
	}
      cout<<endl<<"ElasticNet::  N of beads = "<<_beadList.size()<<endl;
      _structure.setBeadList(_beadList);
    }
    
    //STA COSA NON SERVE, TANTO NON USO MAI STA VARIABILE: SE HO CA FRA I BEAD BENE< ALTRIMENTI CICCE
    // if(!strncmp(stringa,"PROTEIN",7)){
    //   cout<<"ElasticNet:: I want your proteins baby..."<<endl;
    //   _PROTEIN=true;
    // }

    if(!strncmp(stringa,"OUTBEADS:",9)){
      cout<<"ElasticNet::   "<<stringa<<" ";
      char name[4];
      _outBeadList.clear();
      while (iss >> name)
	{
	  cout<<name<<" ";
	  _outBeadList.push_back(name);
	}
      cout<<endl<<"ElasticNet::  N of beads for output = "<<_outBeadList.size()<<endl;
      if((_outBeadList.size()>0)&&(_outBeadList.size()<_beadList.size())) _TRACCIAMENTO=true;
      if((_outBeadList.size()>0)&&(_beadList.size()<1)){ _TRACCIAMENTO=true; cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;}
      //      _structure.setBeadList(_outBeadList);
    }
    
    if(!strncmp(stringa,"OUTRES:",6)){
      cout<<"ElasticNet::   "<<stringa<<" ";
      int num;
      _outResList.clear();
      while (iss >> num)
	{
	  cout<<num<<" ";
	  _outResList.push_back(num);
	}
      cout<<endl<<"ElasticNet::  N of res for output = "<<_outResList.size()<<endl;
      _TRACCIAMENTO=true;
      _OUTRES=true;
    }
    
    if(!strncmp(stringa,"COM",3)){
      _COM=true;
      cout<<" >>> Beads in the centers of mass for sugar and base <<< "<<endl;
    }

    if(!strncmp(stringa,"FAST",4)){
      _FAST=true;
      cout<<" >>> I'm gonna be faster than light! <<< "<<endl;
    }

    if(!strncmp(stringa,"DPF",3)){
      _DPF=true;
      cout<<" >>> Dump Partial Fluctuations <<< "<<endl;
    }

    if(!strncmp(stringa,"DUMPCOV",7)){
      _DUMPCOVMAT=true;
      _DUMPREDCOVMAT=true;
      cout<<" >>> Dump covariance matrix <<< "<<endl;
    }

    if(!strncmp(stringa,"DISTFLUC",7)){
      _DUMPDISTFLUC=true;
      cout<<" >>> Dump distance fluctuations <<< "<<endl;
    }

    if(!strncmp(stringa,"RANDOM",6)){
      _RANDOM=true;
      double seed;
      iss>>argument; sscanf(argument,"%lf",&seed);
      srand(seed);
      cout<<" @!#$$$&^%$!  rAnDOm NeTWork  &*#$&#%&^@## (seed = "<<seed<<")"<<endl;
    }

    if(!strncmp(stringa,"RARAND",6)){
      _RARAND=true;
      double seed;
      iss>>argument; sscanf(argument,"%lf",&seed);
      srand(seed);
      cout<<" @!#$$$&^%$!  rAnDOm NeTWork  &*#$&#%&^@## (seed = "<<seed<<")"<<endl;
    }

    if(!strncmp(stringa,"EXP",3)){
      _EXP=true;
      cout<<" >>> exponential <<< "<<endl;
    }


    if(!strncmp(stringa,"DUMP_MODES",10)){
      _DUMPMODES=true;
      cout<<" >>> Dump top modes <<< "<<endl;
    }

    // if(!strncmp(stringa,"DUMP_EIGENVALUES",16)) sscanf(argument,"%d",&DUMP_EIGENVALUES);
    // if(!strncmp(stringa,"DUMP_FULL_COVMAT",16)) sscanf(argument,"%d",&DUMP_FULL_COVMAT);
    // if(!strncmp(stringa,"DUMP_REDUCED_COVMAT",19)) sscanf(argument,"%d",&DUMP_REDUCED_COVMAT);
    // if(!strncmp(stringa,"DUMP_NORMALISED_REDUCED_COVMAT",30)) sscanf(argument,"%d",&DUMP_NORMALISED_REDUCED_COVMAT);
    // if(!strncmp(stringa,"DUMP_MEAN_SQUARE_DISPL",22)) sscanf(argument,"%d",&DUMP_MEAN_SQUARE_DISPL);


    
  }//END LOOP ON FILE LINES

  //  fclose(fp);
  file_in.close();
  return;
}


int ElasticNet::tracciamento(){
  // MODIFICO PERCHE' CAMBI ANCHE LA STRUTTURA DI RIFERIMENTO!!!
  int *cacca=i1t(_N);
  int *prog=i1t(_N);
  int Na=0;

  //*********


  cout<<" find the out-beads "<<endl;

  if(_outBeadList.size()>_beadList.size()) {cout<<"WARNING: OUTBEADS are more than beads!"<<endl;}
  if(_outBeadList.size()==0) {cout<<"WARNING: OUTBEADS is empty!"<<endl;}

  for(int i=0;i<_N;++i){
    string b=_structure.getBead(i).getAtomType();
    //    cout<<b<<endl;
    for (int k=0;k<_outBeadList.size();++k)
      {
	string obead_k=_outBeadList.at(k);
	//	cout<<obead_k<<endl;
	replace(obead_k.begin(), obead_k.end(), 'p', '\'');
	//	cout<<obead_k<<endl;
	if (obead_k.compare(b)==0) { cacca[i]=1;}
      }
  }
  //*********


  if(_OUTRES){
    cout<<"LOOKING FOR THE NODES IN OUTRES"<<endl;
    for(int i=0;i<_N;++i){
      //    cout<<_structure.getBead(i).getResNum()<<endl;
      if(cacca[i]==0) continue; //se non e' un outbead non mi preoccupo nemmeno
      if(cacca[i]==1) cacca[i]=0; //se si lo metto uguale a zero e...
      int n=_structure.getBead(i).getResNum();
      //      cout<<" "<<n<<" "<<endl;
      for (int k=0;k<_outResList.size();++k)
	{                                       //...lo rimetto uguale a uno 
	  if (_outResList.at(k)==n) {cacca[i]=1;cout<<"###"<<n<<endl; /*cout<<"ok"<<endl;*/} // solo se e' anche nella lista 
	}                                       // dei residui per l'output
    }
  }
  

  /* save the number of out-beads progressively on an array */
  prog[0]=cacca[0]-1;
  for(int i=1;i<_N;++i){
    prog[i]=prog[i-1]+cacca[i];
  }
  Na=prog[_N-1]+1;
  int Nb=_N-Na;

  if(Nb==0) {cout<<endl<<"!!!"<<endl<<"!!! ElasticNet:: tracciamento -> no need for that"<<endl<<"!!!"<<endl; return 0; }

  if(Na==0) {cout<<endl<<"!!!"<<endl<<"!!! ElasticNet:: tracciamento -> ERROR with the OUTBEADS list"<<endl<<"!!!"<<endl; return -1; }

  cout<<"Na="<<Na<<" Nb="<<Nb<<endl;

  double **Ma=d2t(3*Na,3*Na);
  double **Mb=d2t(3*Nb,3*Nb);
  double **V=d2t(3*Na,3*Nb);

   for(int i=0;i<_N;++i){
     for(int j=0;j<_N;++j){
       ////////////////////////////////
       /* Wanted-wanted => matrix Ma */
       if((cacca[i]==1) && (cacca[j]==1)){ 
	 for(int mu=0;mu<3;++mu){
	   for(int nu=0;nu<3;++nu){
	     Ma[3*prog[i]+mu][3*prog[j]+nu]=_intMatrix.GetElement(3*i+mu,3*j+nu);
	   }//nu
	 }//mu
       }//if
       /* not wanted - not wanted => matrix Mb */
       if((cacca[i]==0) && (cacca[j]==0)){ 
	 for(int mu=0;mu<3;++mu){
	   for(int nu=0;nu<3;++nu){
	     Mb[3*(i-prog[i]-1)+mu][3*(j-prog[j]-1)+nu]=_intMatrix.GetElement(3*i+mu,3*j+nu);
	   }//nu
	 }//mu
       }//if       
       /* wanted - not wanted => matrix V */
       if((cacca[i]==1)&&(cacca[j]==0)){ 
	 for(int mu=0;mu<3;++mu){
	   for(int nu=0;nu<3;++nu){
	     V[3*prog[i]+mu][3*(j-prog[j]-1)+nu]=_intMatrix.GetElement(3*i+mu,3*j+nu);
	   }//nu
	 }//mu
       }//if       
       ////////////////////////////////
     }//j
   }//i
  
   double**invMb=d2t(3*Nb,3*Nb);

   cout<<"ElasticNet::Tracciamento - inverto la matrice"<<endl;
   if(invert_symmetric_matrix_lapack(Mb,3*Nb,invMb)<0){cout<<"ERRORE NELL'INVERSIONE"<<endl;return -1;}
   
   const time_t begin_time=time(0);
   //   cout<<"0^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
   IntMatrix schifo(Na,Ma);
   cout<<"inizio lo schifo"<<endl;
   double**VinvMb;
   VinvMb=d2t(3*Na,3*Nb);
   cout<<"cacca"<<endl;
#pragma omp parallel for
   for(int i=0;i<3*Na;++i){
     for(int l=0; l<3*Nb; ++l){
       //#pragma omp parallel for
       for(int k=0; k<3*Nb; ++k){    
	 VinvMb[i][l]+=V[i][k]*invMb[k][l];
       }
     }
   }
   cout<<"cacca"<<endl;
   double**merda;
   merda=d2t(3*Na,3*Na);
   cout << "tempo nel ciclo = "<<difftime(time(0),begin_time)<<endl;
   clock_t new_time=clock();
   //   cout<<"1^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
   #pragma omp parallel for
   for(int i=0;i<3*Na;++i){
     for(int j=0; j<3*Na; ++j){
       for(int l=0; l<3*Nb; ++l){    
	 merda[i][j]+=VinvMb[i][l]*V[j][l];
       }
     }
   }
   //   cout << float( clock () - new_time ) /  CLOCKS_PER_SEC<<endl;
   new_time=clock();
   //   cout<<"2^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
   for(int i=0;i<3*Na;++i){
     for(int j=0;j<3*Na;++j){
       Ma[i][j]-=merda[i][j];
     }
   }
   //   cout << float( clock () - new_time ) /  CLOCKS_PER_SEC<<endl;

   //   cout<<"3^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
   cout<<"finito lo schifo"<<endl;

   _intMatrix=(IntMatrix(Na,Ma));

  //*** Creo la nuova struttura di riferimento: ***
  cout<<"*** Creo la nuova struttura di riferimento: ***"<<endl;
  Structure3d new_struc;
  BeadList lista(_outBeadList);
  lista.setResNum(_outResList);
  new_struc.setBeadList(lista);
  for(int i=0;i<_N;++i){
    //    if(cacca[i]==0) cout<<"no! ";
    new_struc.addBead(_structure.getBead(i));
    //    cout<<i<<" cacca "<<new_struc.getSize()<<endl;   
  }
  cout<<new_struc.getSize()<<endl;
  setStructure(new_struc);
  cout<<"*** *** *** *** *** *** *** *** *** *** *** ***"<<endl;
  //*** *** *** *** *** *** *** *** *** *** *** ***
  _structure.dumpPDBFile("strutture.pdb");


   free_d2t(Ma);
   free_d2t(Mb);
   free_d2t(V);

   free_d2t(invMb);
   free_d2t(merda);
   free_d2t(VinvMb);

   free_i1t(cacca); 
   free_i1t(prog);

   return 0; //success
}
   
void ElasticNet::dumpDistVecFluc(const char*fname){
  /*************************/
  char filename[200];
  sprintf(filename,"%s_dist_vec_fluc.dat",fname);
  ofstream fp;
  fp.open(filename);
  for(int i=0; i < _N; i++){
    for(int j=0; j < _N; j++){
      double Ctilde_ii=0.;
      double Ctilde_jj=0.;
      double Ctilde_ij=0.;
      for(int mu=0;mu<3;mu++){
	Ctilde_ii+=_covMatrix.GetElement(3*i+mu,3*i+mu);
	Ctilde_jj+=_covMatrix.GetElement(3*j+mu,3*j+mu);
	Ctilde_ij+=_covMatrix.GetElement(3*i+mu,3*j+mu);
      }      
      double sigma=Ctilde_ii+Ctilde_jj-2*Ctilde_ij;
      fp<<i<<" "<<j<<" "<<sigma<<" "
	 <<_structure.getBead(i).getAtomType()<<" "
	 <<_structure.getBead(j).getAtomType()<<" "
	 <<_structure.getBead(i).getResNum()<<" "
	 <<_structure.getBead(j).getResNum()<<" "
<<endl;
    }
  }
fp.close();
printf("Fluctuations of distance vectors between beads written on file %s\n",filename);

}//ENDFUNCTION

void ElasticNet::dumpDistFluc(const char*fname){
  /*************************/
  char filename[200];
  sprintf(filename,"%s_dist_fluc.dat",fname);
  ofstream fp;
  fp.open(filename);
  for(int i=0; i < _N; i++){
    for(int j=0; j < _N; j++){
      double sigma=getDistFluc(i,j);
      fp<<i<<" "<<j<<" "<<sigma<<" "
	 <<_structure.getBead(i).getAtomType()<<" "
	 <<_structure.getBead(j).getAtomType()<<" "
	 <<_structure.getBead(i).getResNum()<<" "
	 <<_structure.getBead(j).getResNum()<<" "
<<endl;
    }
  }
fp.close();
printf("Fluctuations of distance vectors between beads written on file %s\n",filename);

}//ENDFUNCTION


double ElasticNet::getDistFluc(int i, int j){
  if(i==j) return 0;
  double sigma=0;

  Vector3d d_ij=_structure.getDistanceVector(i,j);
  d_ij/=_structure.getDistance(i,j);

  if(_covMatrix.GetSize()>0){//Se ho calcolato covmatrix la uso ...
    for(int mu=0;mu<3;mu++){
      for(int nu=0;nu<3;nu++){
	sigma+=d_ij[mu]*d_ij[nu]*(
				  _covMatrix.GetElement(3*i+mu,3*i+nu)+
				  _covMatrix.GetElement(3*j+mu,3*j+nu)-
				  _covMatrix.GetElement(3*i+mu,3*j+nu)-
				  _covMatrix.GetElement(3*j+mu,3*i+nu));
      }
    }
    return sigma;
  }
  else{                     // ... altrimenti ricalcolo solo gli elem che mi servono!
    for(int mu=0;mu<3;mu++){
      for(int nu=0;nu<3;nu++){
	double C_ii_munu=0;
	double C_jj_munu=0;
	double C_ij_munu=0;
	double C_ji_munu=0;
	for(int alpha=0;alpha<3*_N-6;alpha++){
	  double temp=_intMatrix.GetEigenval(alpha);
	  if (temp<0.000001) cout<<"WARNING: eigenvalue very small"<<endl;
	  double lambda=1./temp;
	  Vector3d v_i=_intMatrix.GetEigenvec(alpha).at(i);
	  Vector3d v_j=_intMatrix.GetEigenvec(alpha).at(j);
	  C_ii_munu+=lambda*v_i[mu]*v_i[nu];
	  C_jj_munu+=lambda*v_j[mu]*v_j[nu];
	  C_ij_munu+=lambda*v_i[mu]*v_j[nu];
	  C_ji_munu+=lambda*v_j[mu]*v_i[nu];
	}
	sigma+=d_ij[mu]*d_ij[nu]*(
				  C_ii_munu+
				  C_jj_munu-
				  C_ij_munu-
				  C_ji_munu);
	
      }
    }
    return sigma;
  }
}



void ElasticNet::dump_top_modes(int ntop, double amp,const char *name){
  //this subroutine prints on files the trajectories correspondent to the first ntop modes
  char filename[200];
  ofstream modeout;
  double PI=3.1415;
  int tmax=50;
  Structure3d temp_struc=_structure;
  
  if (ntop<0)ntop=_covMatrix.GetSize();

  for (int n=0; n<ntop; n++){
    sprintf(filename,"%s_mode_%d.pdb",name,n);
    
    cout<<"Writing mode on file "<<filename<<endl;
    modeout.open(filename);
    
    for(int t=0;t<tmax;t++){
      
      double prefac=amp*cos((2.*PI*t)/tmax);
      int size=_structure.getSize();
      for(int i=0; i<size;i++){
	Vector3d pos = _structure.getBead(i).getCoordinates();
	pos+=prefac*_covMatrix.GetEigenvec(3*size-1-n).at(i);
	temp_struc.getBead(i).setCoordinates(pos);
      }
      temp_struc.dumpPDBFile(modeout);
      // 	for(int mu=0; mu<3;mu++){
      // 	  modeout<<nodes[i].pos[mu]+prefac*eigenvectors[3*mol_len-n-1][i+mu*mol_len]<<" ";}
    }
    modeout.close();
  }
  return;
}
  
void ElasticNet::dumpNearestNeighbours(const char*name){
  char filename[200];
  ofstream fout;
  sprintf(filename,"%s_near_neigh.dat",name);
  fout.open(filename);
  for(int i=0;i<_N;++i){
    int NN =0;
    for (int j=0;j<_N;++j){
      if(i==j) continue;
      if(_cmap[i][j]!=0) NN++;
    }
    double iNN=1./((double) NN);
    fout<<i<<" "<<iNN<<" "<<NN<<" "<<_structure.getBead(i).getResName()<<endl;
  }
  fout.close();
  return;
}


void ElasticNet::computeForCovMatrix(){
  // cout<<"ElasticNet::computeForCovMatrix *** inizio la function. Ora quadro"<<endl;
  // Matrix tempMatSQ=dot_product(_intMatrix,_intMatrix);
  // cout<<"ElasticNet::computeForCovMatrix *** Ho quadrato, ora cubo"<<endl;
  // Matrix tempMatCU=dot_product(tempMatSQ,_intMatrix);
  // cout<<"ElasticNet::computeForCovMatrix *** Ho cubato"<<endl;
  double **forcovmat=d2t(_N,_N);
  for(int i=0;i<_N;i++){
    for(int j=i;j<_N;j++){
      for(int mu=0;mu<3;mu++){
	forcovmat[i][j]+=_intMatrix.GetElement(3*i+mu,3*j+mu);
      }//enddo mu
      forcovmat[i][j]*=4;
      if (i!=j) forcovmat[j][i]=forcovmat[i][j];
    }//endo i
  }//endo j
  cout<<"ElasticNet::computeForCovMatrix *** Ora setto"<<endl;
  _forCovMatrix.SetMatrix(_N,forcovmat);
  _forCovMatrix.Decompose();
  //  _forCovMatrix.dumpEigenvalues(
}

void ElasticNet::dumpForCovMatrix(const char* fname){
  _forCovMatrix.dumpEigenvalues(fname);
  _forCovMatrix.dumpTopVectors(_ntop_vectors,fname);
  _forCovMatrix.dumpMSF(fname);
  _forCovMatrix.dumpMatrix(fname);
  cout<<"FINISH COVFOR"<<endl;
}

void ElasticNet::dumpMSF(const char*name){
  char filename[200];
  ofstream fout;
  sprintf(filename,"%s_mean_square_displ.dat",name);
  fout.open(filename);
  fout<<"# index | MSF | atom type | residue name | res number"<<endl;
  for(int i=0; i < _N; i++){
    double temp=_intMatrix.GetMSF(i);
    //    fprintf(fp,"%4d %e\n",i,temp);
    Bead b=_structure.getBead(i);
    fout<<i<<" "<<temp<<" "<<b.getAtomType()<<" "<<b.getResName()<<" "<<b.getResNum()<<endl;
  }
  fout.close();
  printf("Beads' mean square displacement written to file %s\n",filename);
}

