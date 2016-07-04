/*
  Created by Giovanni Pinamonti
  PhD student @
  Scuola Internazionale Superiore di Studi Avanzati, Trieste, Italy
  November 15th 2013
 */

//extern "C" {
//#include <mkl.h>
//}


#include <fstream>
#include <iostream>
#include <sstream>
#include <ctime>

#include "PrincipalComp.h"
#include "my_malloc.h"
#include "io.h"
#include "Vector3d.h"


using namespace std;

//########################################
//Constructors for the class PrincipalComp
//########################################
PrincipalComp::PrincipalComp():_ntmax(0),_ntmin(0){_N=0;_size=0;}

PrincipalComp::PrincipalComp(int n):_ntmax(0),_ntmin(0){
  _N=n; _size=3*n;}

PrincipalComp::PrincipalComp(string fname):_ntmax(0),_ntmin(0){
  readTrajectory(fname);}

PrincipalComp::PrincipalComp(string fname, std::vector<std::string> beadlist):_ntmax(0),_ntmin(0){
  setBeadList(beadlist);
  readTrajectory(fname.c_str());}

PrincipalComp::PrincipalComp(string fname, BeadList beadlist):_ntmax(0),_ntmin(0){
  setBeadList(beadlist);
  readTrajectory(fname.c_str());}


PrincipalComp::~PrincipalComp() { }

//########################################


void PrincipalComp::readTrajectory(const char* fname){
//   vector<string> list;
//   return readTrajectory(fname, list);
// }

// void PrincipalComp::readTrajectory(const char* fname, std::vector<std::string> beadlist){
  
  Structure3d structure(0);
  
  structure.setBeadList(_beadList);

  cout<<"PrincipalComp: Opening file"<<endl;
  FILE* fp;
  if ((fp = fopen (fname, "r")) == NULL) {
    fprintf (stderr,"Could not open file %s.\n.Exiting.\n", fname);
    exit (1);    }

  cout<<"PrincipalComp: ...reading..."<<endl;
  int error=0;
  int Nsteps=0;
  int Nsteps_read=0;

  error=structure.readFromPDBFile(fp);
  _N=structure.getSize();
  if((_size!=3*_N)&&(_size>0)) cout<<"WARNING: old size ("<<_size<<") different from the new one ("<<3*_N<<")"<<endl;
  _size=3*_N;
  cout<<"taglia : "<<_size<<endl;
  cout<<"error : "<<error<<endl;
  
  _refStructure=structure;
  _mean=structure;

  
 double **matrix=d2t(_size,_size);
 double *mean=d1t(_size);

 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  while(error==0){                                     // * * *   
    ++Nsteps;					       // * * * 
    //cout<<"Step "<<Nsteps<<endl;    // * * * 
    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    if(Nsteps<_ntmin) {
      error=structure.readFromPDBFile(fp);  	       // * * *
      continue;}
    if(Nsteps%100==0) cout<<"Step "<<Nsteps<<endl;    // * * * 
    //    cout<<_ntmax<<">"<<Nsteps<<endl;
    if((Nsteps>_ntmax)&&(_ntmax>0)) break;
    Nsteps_read++;
    for(int i=0; i<_N; ++i){
      Vector3d temp_vec_i= structure.getBead(i).getCoordinates();      
      for(int mu=0; mu<3; ++mu){
	mean[3*i+mu]+= temp_vec_i[mu];
      }
      for(int j=0; j<_N; ++j){
	Vector3d temp_vec_j=structure.getBead(j).getCoordinates();
	for(int mu=0; mu<3; ++mu){
	  for(int nu=0; nu<3; ++nu){
	    matrix[3*i+mu][3*j+nu]+= temp_vec_i[mu]*temp_vec_j[nu];
	  }//enddo nu
	}//enddo mu
      }//enddo j
    }//enddo i


    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    error=structure.readFromPDBFile(fp);  	       // * * *
  }//end while  			 	       // * * *
  //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  fclose(fp);

  //normalizing mean and matrix
  for(int i=0; i<_size; ++i){
    mean[i]/=Nsteps_read;
    for(int j=0; j<_size; ++j){
      matrix[i][j]/=Nsteps_read;
    }//enddo
   

  }//enddo
  //"centering" the covariance matrix
  // if (_REF){
  //   for(int i=0; i<_size; i++){
  //     for(int j=0; j<_size; j++){
  // 	_matrix[i][j]+=-_mean[i]*_refer[j]-_mean[j]*_refer[i]+_refer[i]*_refer[j];
  // 	//refer[i]*refer[j];
  //     }//enddo
  //   }//enddo
  // }
  
  for(int i=0; i<_N; ++i){
    double x=mean[3*i+0];
    double y=mean[3*i+1];
    double z=mean[3*i+2];
    Vector3d tempvec(x,y,z);
    _mean.getBead(i).setCoordinates(tempvec);
  }
  
  for(int i=0; i<_size; ++i){
    for(int j=0; j<_size; ++j){
      matrix[i][j]-=mean[i]*mean[j];
    }//enddo
  }//enddo
  
  _covMat=(CovMatrix(_N,matrix));

  free_d1t(mean);
  free_d2t(matrix);

  cout<<"...done!"<<endl<<endl;
}

void PrincipalComp::Diagonalize(){
  cout<<"Diagonalizing the covariance matrix..."<<endl;
  //  cout<<_intMat.GetN();
  _covMat.Decompose();
  cout<<"...done!"<<endl<<endl;
}

void PrincipalComp::Dump(const char*fname){
  cout<<"OUTPUT:"<<endl;
  //_covMat.Dump(fname);
  _covMat.dumpTopVectors(_NTOP,fname);
  _covMat.dumpEigenvalues(fname);
  dumpMSF(fname);
}

void PrincipalComp::DumpCov(const char*fname){
  string covname(fname);
    _covMat.dumpFullMatrix(covname+"_cov");
    _covMat.dumpReducedMatrix(covname+"_cov");
    _covMat.dumpNormalizedReducedMatrix(covname+"_cov");
}


//void PrincipalComp::DumpDistC2(const char*fname){
//

void PrincipalComp::dumpDistFluc(const char*fname){
  /*************************/
  char filename[200];
  sprintf(filename,"%s_dist_fluc.dat",fname);
  ofstream fp;
  fp.open(filename);
  for(int i=0; i < _N; i++){
    for(int j=0; j < _N; j++){
      //	  if(i==j) fp<<3*i+mu<<" "<<3*j+nu<<" "<<0.0<<endl;
      //	  else fp<<3*i+mu<<" "<<3*j+nu<<" "<<_intMat.GetElement(3*i+mu,3*j+nu)<<endl;
      
      double Ctilde_ii=0.;
      double Ctilde_jj=0.;
      double Ctilde_ij=0.;
      for(int mu=0;mu<3;mu++){
	Ctilde_ii+=_covMat.GetElement(3*i+mu,3*i+mu);
	Ctilde_jj+=_covMat.GetElement(3*j+mu,3*j+mu);
	Ctilde_ij+=_covMat.GetElement(3*i+mu,3*j+mu);
      }      
      double sigma=Ctilde_ii+Ctilde_jj-2*Ctilde_ij;
      fp<<i<<" "<<j<<" "<<sigma<<" "
       <<_refStructure.getBead(i).getAtomType()<<" "
       <<_refStructure.getBead(j).getAtomType()<<" "
       <<_refStructure.getBead(i).getResNum()<<" "
       <<_refStructure.getBead(j).getResNum()<<" "
<<endl;
    }
  }
fp.close();
printf("Fluctuations of distances between bases written on file %s\n",filename);

}//ENDFUNCTION


void PrincipalComp::dumpMSF(const char *name){
  char filename[200];
  ofstream fp;
  sprintf(filename,"%s_mean_square_displ.dat",name);
  fp.open(filename);
  for(int i=0; i < _covMat.GetN(); i++){
    double temp=_covMat.GetMSF(i);
    fp <<i <<" "
       <<temp<<" "
       <<_refStructure.getBead(i).getAtomType()<<" "
       <<_refStructure.getBead(i).getResNum()<<" "
       <<endl;
  }
  fp.close();
  cout<<"Beads' mean square displacement written to file "<<filename<<endl;
}

void PrincipalComp::dumpMSF_partial(const char *name){
  char filename[200];
  ofstream fp;
  
  sprintf(filename,"%s_mean_square_fluc_part.dat",name);
  fp.open(filename);
  for(int i=0; i < _covMat.GetN(); i++){
    fp <<i <<" "
       <<_refStructure.getBead(i).getAtomType()<<" "
       <<_refStructure.getBead(i).getResNum()<<" ";
    for(int mu=0;mu<3;mu++){
      for(int nu=mu;nu<3;nu++){    
	double temp=_covMat.GetMSF(i,mu,nu);
	fp<<temp<<" ";
      }}
    fp<<endl;
    }
    fp.close();
    cout<<"Beads' mean square partial fluctuations written to file "<<filename<<endl;
  
}


void PrincipalComp::DumpInt(const char*fname){

   /*************************/
  /* print the full matrix */
  char filename[200];
  sprintf(filename,"%s_int_full_matrix.dat",fname);
  ofstream fp;
  fp.open(filename);
  for(int i=0; i < _N; i++){
    for(int j=0; j < _N; j++){
      for (int mu=0; mu < 3; mu++){
	for (int nu=0; nu < 3; nu++){
	  //	  if(i==j) fp<<3*i+mu<<" "<<3*j+nu<<" "<<0.0<<endl;
	  //	  else fp<<3*i+mu<<" "<<3*j+nu<<" "<<_intMat.GetElement(3*i+mu,3*j+nu)<<endl;
	  fp<<3*i+mu<<" "<<3*j+nu<<" "<<_intMat.GetElement(3*i+mu,3*j+nu)<<endl;
	}
      }
    }//enddo j
  }//enddo i
  fp.close();
  printf("Full interaction matrix written to file %s\n",filename);
  
   /****************************/
  /* print the reduced matrix */
  sprintf(filename,"%s_int_reduced_matrix.dat",fname);
  fp.open(filename);
  for(int i=0; i < _N; i++){
    for(int j=0; j < _N; j++){
      if (i==j) {      fp<<i<<" "<<j<<" "<<0.0<<endl; continue;}
      double temp=0;
      for (int mu=0; mu < 3; mu++){
        temp += _intMat.GetElement(3*i+mu,3*j+mu);
      }
      fp<<i<<" "<<j<<" "<<temp<<endl;
    }//enddo j
  }//enddo i
  fp.close();
  printf("Reduced interaction matrix written to file %s\n",filename);
  
   /**************/
  /* print shit */
  sprintf(filename,"%s_int_shit.dat",fname);
  fp.open(filename);
  for(int i=0;i<_N;++i){
    for(int j=i+1; j<_N;++j){
      double d=_refStructure.getDistance(i,j);
      double temp=0;
      for(int mu=0; mu < 3; mu++){
        temp += _intMat.GetElement(3*i+mu,3*j+mu);
      }
      string itype=_refStructure.getBead(i).getAtomType();
      string jtype=_refStructure.getBead(j).getAtomType();
      fp<<d<<" "<<temp<<" "<<itype<<" "<<jtype<<endl;
    }    
  }

  cout<<"Other random stuff written on file "<<filename<<endl;

}//ENDFUNCTION

//#############################################

// double CMatrix:: Converge(double *coordinates, int n){
//   double scal_prod=0;
//   for(int i=0;i<_N;i++){
//     for(int mu=0;mu<3;mu++){
//       scal_prod+=_eigenvectors[_size-n-1][3*i+mu]*coordinates[3*i+mu];
//     }//ENDDO mu
//   }//ENDDO i
//   return scal_prod;
// }

// void CMatrix:: SetReference(double *coordinates){
//   _refer=d1t(_size);
//   for(int i=0;i<_size;i++){
//     _refer[i]=coordinates[i];
//   }
//   _REF=true;
// }


void PrincipalComp::ComputeInt(){
  _intMat=_covMat.GetInverse();
  //  _covMat.GetInverse(_intMat);
}

//###################################################################

void PrincipalComp::Project(int n_comp,const char*iname, const char*oname){
  
  cout<<"PrincipalComp: reading trajectory..."<<endl;
  Structure3d structure(0);
  structure.setBeadList(_beadList);

  //open input file
  FILE* fin;
  if ((fin = fopen (iname, "r")) == NULL) {
    fprintf (stderr,"Could not open file %s.\n.Exiting.\n", iname);
    exit (1);    }  
    
  //open outputfile
  char filename[200];
  sprintf(filename,"%s_projection.dat",oname);
  ofstream fout;
  fout.open(filename);
  
  //check if the input trajectory is correct
  int error=0;
  int Nsteps=0;
  error=structure.readFromPDBFile(fin);
  int n_beads=structure.getSize();
  if (_N!=n_beads) { cout<<"PrincipalComp.Project:: ERROR!  ---  the number of beads in the traj I'm reading is not consistent with the principal components I've stored"<<endl; return;}
  
  while(error==0){
    ++Nsteps;
    if(Nsteps%1000==0) cout<<"Frame: "<<Nsteps<<endl;
    fout<<Nsteps<<" ";
    /* qui stampo la sequenza */
    int oldres=-1;
    for(int i=0;i<_N;++i){
      if(structure.getBead(i).getResNum()!=oldres){
	fout<<structure.getBead(i).getResName();
	oldres=structure.getBead(i).getResNum();
      }
    }
    fout<<" ";
    /* qui proietto e stampo */
    for (int n=0; n<n_comp; ++n){
      double scal=0;
      for(int i=0; i<_N; ++i){
	Vector3d temp_pos_i= structure.getBead(i).getCoordinates();
	temp_pos_i-= _mean.getBead(i).getCoordinates();
	scal+=temp_pos_i.dot(_covMat.GetEigenvec(_size-1-n).at(i));
      }
      fout<<scal<<" ";
    }//enddo i
    /* fatto, vado a capo */
    fout<<structure.getTitle();
    //    fout<<endl;
    error=structure.readFromPDBFile(fin);
  }//end while
  
  fclose(fin);
  fout.close();

  cout<<"...done!"<<endl<<endl;  
  
}


//  ATTENZIONE QUI STA LEGGENDO EVAL E EVEC COME SE FOSSERO DA UNA PCA!!!! (giusto?)
void PrincipalComp::readFromFile(const char *basename, int numModes) {
  // questa funzione dovrebbe farla forse Matrix.cc per avere una gestione piu' sensata dell'allocazione della memoria dei vector
  clock_t start_c=clock();
  time_t start_t=time(0);
  char fname[1024];
  vector<double> eigenvalues;
  vector<vector <Vector3d> > eigenvectors;
  cout<<"ATTENZIONE QUI STA LEGGENDO EVAL E EVEC COME SE FOSSERO DA UNA PCA!!!! (giusto?)"<<endl;

  if(_size==0) {
    _N=numModes/3; 
    _size=3*_N; 
    cout<<"Old size was 0, setting it to "<<_size<<endl;
    if(_size!=numModes) {cout<<"  ******** ERROR numModes should be a multiple of 3!! ***********"<<endl; exit(1);}
  }
  else if(_size!=numModes) {
    _N=numModes/3; 
    _size=3*_N; 
    cout<<" - - - WARNING: Old size was different from numModes, setting it to "<<_size<<" - - - "<<endl;
    if(_size!=numModes) {cout<<"  ******** ERROR numModes should be a multiple of 3!! ***********"<<endl; exit(1);}
  }

  cout<<"%%%"<<basename<<endl;
  sprintf(fname,"%s%s", basename, "_eigenvalues.dat");
  ifstream eigvalFile(fname);
  if (!eigvalFile.is_open()){
    cerr << CLRLINE << "PrincipalComp.readFromFile:: FATAL ERROR. Could not open " << fname << endl;
    exit(1);
  }  
  try{
    eigenvalues.reserve(numModes);
  }catch(std::bad_alloc const&){
    cerr << CLRLINE << "NormalMode eigenvalues memory allocation failed!" << endl << flush;
    exit(1);
  }  
  while(eigvalFile.good()){
    if (eigenvalues.size() == numModes) break;
    int index;
    double val;
    eigvalFile >> index >> val;
    eigenvalues.push_back(val);
  }  
  if(eigenvalues.size() < numModes){
    cerr << CLRLINE << "Error. Could read only " <<  eigenvalues.size() << " eigenvalues. I wanted "<<numModes<< endl;
    exit(1);
  }
  try{
    eigenvectors.resize(numModes);
    for (int i = 0; i < numModes; ++i) {
      
      sprintf(fname, "%s%s%ld%s",basename,"_eigenvector_",i,".dat");
      ifstream eigvecFile(fname);
      if (!eigvecFile.is_open()){
	cerr << CLRLINE << "Fatal error. Could not open " << fname << endl<< flush;
	exit(1);
      }   
      eigenvectors[i].reserve(_N);
      for (int j = 0; j < _N; ++j) {
	int index;
	Vector3d p3d;
	eigvecFile >> index >> p3d;
	eigenvectors[i].push_back(p3d);	
      }      
      if(eigvecFile.fail()){
	cerr << CLRLINE << "Fatal error reading " << fname << endl << flush;
	exit(1);
      }
      //	if(printProgress) printProgressBar(i,numModes,std::cout);      
    }    
  }catch(std::bad_alloc const&){
    cerr << CLRLINE << "NormalMode eigenvector memory allocation failed!" << endl << flush;
    exit(1);
  }  
  cout<<"PrincipalComp:ReadFromFile->Compose"<<endl;
  //  _covMat=(CovMatrix(_N));
  _covMat.SetN(_N);
  _covMat.Compose(eigenvalues,eigenvectors);

}


void PrincipalComp::dump_top_modes(int ntop, double amp,const char *name){
  //this subroutine prints on files the trajectories correspondent to the first ntop modes
  char filename[200];
  ofstream modeout;
  double PI=3.1415;
  int tmax=50;
  Structure3d temp_struc=_refStructure;
  
  for (int n=0; n<ntop; n++){
    sprintf(filename,"%s_mode_%d.pdb",name,n);
    
    cout<<"Writing mode on file "<<filename<<endl;
    modeout.open(filename);
    
    for(int t=0;t<tmax;t++){
      
      double prefac=amp*cos((2.*PI*t)/tmax);
      int size=_refStructure.getSize();
      for(int i=0; i<size;i++){
	Vector3d pos = _refStructure.getBead(i).getCoordinates();
	pos+=prefac*_covMat.GetEigenvec(size-1-n).at(i);
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


void PrincipalComp::readForces(const char* fname){

  cout<<"PrincipalComp: Opening file"<<endl;
  FILE* fp;
  if ((fp = fopen (fname, "r")) == NULL) {
    fprintf (stderr,"Could not open file %s.\n.Exiting.\n", fname);
    exit (1);    }

  cout<<"PrincipalComp: ...reading..."<<endl;
  int error=0;
  int nframes=0;
  int nframes_read=0;
  int nlines=0;
  int step;
  int nforces, temp_nforces;
  char line[1000];

  double **matrix=d2t(3*_N,3*_N);
  //  vector<Vector3d> mean_forces;
  double *mean_forces=d1t(3*_N);
  // for(int i=0;i<_N;i++){
  //   mean_forces.push_back(Vector3d(0,0,0));
  // }
  vector<Vector3d> forces;
  while(fgets(line,1000,fp)!=NULL){
    ++nlines;		 //numero di linee scannate
    //    cout<<"lines="<<nlines<<endl; 
    
    istringstream iss(line); // variabili che 
    char stringa[1000];      // usero' per leggere
    char argument[1000];

    if (!(iss>>stringa)) continue; //butta in stringa la prima merda di iss, se iss e' vuoto passa alla riga dopo
    if (!strncmp(stringa,"natoms=",7)){ //se iss inizia con natoms controlla di avere il giusto numero
      iss>>argument; 
      sscanf(argument,"%d",&nforces);
      if((_N!=nforces)&&(_N>0)) cout<<"WARNING: old dimension ("<<_N<<") different from the one of forces file ("<<nforces<<")"<<endl;
    }//endif strcmp
    if (!strncmp(stringa,"f",1)){// se iss inizia con f allora sto leggendo cose di forze
      iss>>argument;

      if(!strncmp(argument,"(",1)){                           // se dopo f c'e' una parentesi 
	if(sscanf(argument,"(%dx3)",&temp_nforces)!=1)        // allora e' la riga che 
	  cerr<<"WARNING: can't read line"<<endl<<line<<endl; // inizia un nuovo frame con tutte le forze
	if(nforces==0) nforces=temp_nforces;   //se non ho ancora letto forze aggiorno nforces vecchio  
	if(temp_nforces!=nforces) cout<<"WARNIGN: problem with numbers"<<endl; //check nforces
	if(forces.size()>0){ //se non ho gia' letto forze allora aggiorno la matrice
	  //$$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$
	  //$$$ $$$ ---    Qui aggiorno la matrice (meta')  --- $$$ $$$	  
	  if(forces.size()!=nforces) {cerr<<"ERROR: problem with numbers"<<endl; break;} //check nforces
	  // cout<<"F25= "<<forces.at(25)<<endl;
	  // cout<<"F5= "<<forces.at(5)<<endl;
	  // cout<<"F70= "<<forces.at(70)<<endl;
	  for (int i=0;i<nforces;++i){
	    for (int mu=0;mu<3;++mu){
	      mean_forces[3*i+mu]+=forces.at(i).getCoord(mu);
	      for (int j=i;j<nforces;++j){
		for (int nu=0;nu<3;++nu){
	      // double temp=forces.at(i).dot(forces.at(j));
	      // matrix[i][j]+=temp;
		  matrix[3*i+mu][3*j+nu]+=forces.at(i).getCoord(mu)*forces.at(j).getCoord(nu);
		    }
	      }
	    }//endo j	  
	  }//endo i
	  forces.clear();
	  //$$$ $$$                                             $$$ $$$
	  //$$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$ $$$
	}//endif forces.size()>0
	nframes++; //aumento di uno il numero di frame
	//	if((nframes>=_ntmax)&&(_ntmax>0)) break;	
	if(nframes%100==0) cout<<"frame="<<nframes<<endl; //printo il numero di frame letto
	continue;
      }//endif strncmp(
      //      if(nframes<_ntmin) continue;
      //%%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
      //%%% %%% --- LEGGO LE FORZE! --- %%% %%%
      int ibead;
      double fx,fy,fz;
      if(sscanf(line,"      f[%5d]={%12le, %12le, %12le}",&ibead,&fx,&fy,&fz)!=4) 
	{fprintf(stderr,"Fatal error. Could not read forces.\nLine %s\n",line); exit(1);}
      Vector3d tempfor(fx,fy,fz);
      forces.push_back(tempfor);
      //      cout<<"prova: frame="<<nframes<<" ibead="<<ibead<<" x="<<fx<<" y="<<fy<<" z="<<fz<<endl;      
      //%%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
    }//endif strncmp f

    //    ++nframes_read;                                   // * * * 
    //    cout<<"frames read="<<nframes_read<<endl;         // * * *     


  }//endwhile
  fclose(fp);

  nframes_read=nframes;
  // *** normalizing mean and matrix
  for(int i=0; i<3*_N; ++i){
    mean_forces[i]/=nframes_read;
    cout<<"MEAN FORCE: "<<i<<" "<<mean_forces[i]<<endl;
    for(int j=i; j<3*_N; ++j){
      matrix[i][j]/=nframes_read;
      if(j!=i) matrix[j][i]=matrix[i][j];
    }//enddo
  }//enddo

  for(int i=0; i<_size; ++i){
    for(int j=0; j<_size; ++j){
      matrix[i][j]-=mean_forces[i]*mean_forces[j];
    }//enddo
  }//enddo

  double **for_matrix=d2t(_N,_N);

  for(int i=0; i<_N; ++i){
    for(int j=0; j<_N; ++j){
      for (int mu=0;mu<3;++mu){
	for_matrix[i][j]+=matrix[3*i+mu][3*j+mu];
      }
    }//enddo
  }//enddo
  

  cout<<"cacca"<<endl;
  //######################################################################
  //######################################################################
  //######################################################################

  //  double **cacca=d2t(_N,_N);  
  //_forCovMat=(CovMatrix1d(_N,cacca));
  _forCovMat=(CovMatrix1d(_N,matrix));
  //  free_d1t(mean);
  free_d2t(matrix);
  cout<<"...done reading!"<<endl<<endl;

  _forCovMat.Decompose();
  _forCovMat.dumpEigenvalues(fname);
  //  _forCovMat.dumpEigenvectors(fname);
  _forCovMat.dumpTopVectors(_NTOP,fname);
  _forCovMat.dumpMSF(fname);
  _forCovMat.dumpFullMatrix(fname);
  cout<<"FINISH COVFOR"<<endl;
  }//end function readForces

