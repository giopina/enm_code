/*
 * Structure3d.cpp
 *
 *  Created on: Jul 3, 2012
 *      Author: gpolles
 */

#define DEFAULT_CUTOFF 10.0 // default cutoff for neighbors search
#define CHAIN_BREAK_DIST 7.0
#define CHAIN_BREAK_DIST_SQ (CHAIN_BREAK_DIST*CHAIN_BREAK_DIST)

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
//#include "quaternion.h"
//#include "kabsch.cc"
//#include "quaternion.c++"

/* subroutines per l'allineamento con i quaternioni :) :) :) */
#include "qcprot.h"
#include "qcprot.cpp"


#include "Structure3d.h"
//#include "BeadList.h"

using namespace std;

Structure3d::Structure3d() :_size(0),_sizeRes(0), _cutoff(DEFAULT_CUTOFF), _title(""){
}


Structure3d::Structure3d(std::string fname) : _size(0), _sizeRes(0), _cutoff(DEFAULT_CUTOFF), _title(""){  
  readFromPDBFile(fname);
}

Structure3d::Structure3d(double cutoff = DEFAULT_CUTOFF) : _size(0), _sizeRes(0),_cutoff(cutoff), _title(""){
}
#undef DEFAULT_CUTOFF
Structure3d::Structure3d(double cutoff, std::string fname) : _size(0), _sizeRes(0),_cutoff(cutoff),_title(""){
  readFromPDBFile(fname);
}

Structure3d::Structure3d(vector<Bead> beads){
  _beads=beads;
  _size=beads.size();
}

// Structure3d::Structure3d(int N,double *pos) :_cutoff(cutoff){
//   _size=N
//} 


Structure3d::~Structure3d() { }

int Structure3d::readFromPDBFile(const char* fname){
  FILE* fp;
  if ((fp = fopen (fname, "r")) == NULL) {
    fprintf (stderr,"Could not open file %s.\n.Exiting.\n", fname);
    exit (1);    }
  int idum= readFromPDBFile(fp);
  fclose(fp);
  return idum;
}

//return 0 if read succesfully, -1 when file is finished (1 if no atoms in the frame read)
int Structure3d::readFromPDBFile(FILE *fp) {
  // For better performance, the vector capacity is preallocated
  // and expanded when necessary.
  const size_t expandSize = 1000; // size of expansion of capacity

  int discard_alternative_position=0;
  int currentChainId = 0;
  char line[200], *a;
  char field[10];
  double x,y,z,bFactor;
  int old_resNum=-999;


  /* PDB FORMAT

  COLUMNS        DATA TYPE     FIELD        DEFINITION
  -------------------------------------------------------------------------------------
   1 -  6        Record name   "ATOM  "

   7 - 11        Integer       serial       Atom serial number.

  13 - 16        Atom          name         Atom name.

  17             Character     altLoc       Alternate location indicator.

  18 - 20        Residue name  resName      Residue name.

  22             Character     chainID      Chain identifier.

  23 - 26        Integer       resSeq       Residue sequence number.

  27             AChar         iCode        Code for insertion of residues.

  31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.

  39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.

  47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.

  55 - 60        Real(6.2)    occupancy     Occupancy.

  61 - 66        Real(6.2)    tempFactor    Temperature factor.

  77 - 78        LString(2)   element       Element symbol, right-justified.

  79 - 80        LString(2)   charge        Charge on the atom.

   */



  // ### here I clear the beads to cancel the previous frame (if some)
  if(_beads.size()>0) {_beads.clear(); _size=0; _sizeRes=0;}
  // TODO: Do I need to reset also the capacity?     //
  //       It will only be usefull if the structure  //
  //       size changes from one frame to the other  //

  while (fgets(line, 200, fp)){
    if ((strncmp (line, "ENDMDL", 6)==0)||(strncmp (line, "END", 3)==0)){
      if(_size==0) {cout<<"WARNING: no beads read from file"<<endl; return 1;}
      updateNeighbors();
      updateCovBonded();
      return 0;
    }
    

    /* read structure title (if any) */
    if(!strncmp (line, "TITLE", 5)){
      a=line+5;
      _title=std::string(a);
      //      cout<<_title<<endl;
    }

    bool ISHET=false;
    /* check if line buffer begins with "ATOM" */
    if ( (!strncmp (line, "ATOM", 4))||(!strncmp (line, "HETATM", 6)) ){
      if (!strncmp (line, "HETATM", 6)) ISHET=true;
      /* check if character at position 16 is different from ' '
           or 'A'. If so, we are facing an alternative position
           and we will not consider it */
      discard_alternative_position=1;
      if (line[16]==' ') discard_alternative_position=0;
      if (line[16]=='A') discard_alternative_position=0;

      if (discard_alternative_position==1) continue;

      a = line +12; /* advance line buffer and read
      /* Atom type */
      char cacca[4];
      strncpy(field,a,4);field[4]='\0';
      if (sscanf(field,"%s",cacca)!=1){fprintf(stderr,"Fatal error. Could not read atom type in molecule.\nLine %s\n",line); exit(1);}
      std::string type(cacca);
      /* Residue name */
      char name[3];
      a=line+17;
      strncpy(field,a,3);field[3]='\0';
      if (sscanf(field,"%s",name)!=1){fprintf(stderr,"WARNING. Could not read residue name in molecule.\nLine %s\n",line);}
      std::string resName(name);
      /* Residue number */
      a=line+22;
      int resNum;
      strncpy(field,a,4); field[4]='\0';
      if (sscanf(field,"%d",&resNum)!=1){fprintf(stderr,"WARNING. Could not read residue number in molecule.\nLine %s\n",line);}// exit(1);}
      /* Coordinates */
      a=line+30;
      strncpy(field,a,8);field[8]='\0';
      if (sscanf(field,"%lf",&x)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}
      a=line+38;
      strncpy(field,a,8);field[8]='\0';
      if (sscanf(field,"%lf",&y)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}
      a=line+46;
      strncpy(field,a,8);field[8]='\0';
      if (sscanf(field,"%lf",&z)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}
      /* B factor */
      a= line + 60;
      strncpy(field,a,6);field[6]='\0';
      if ( sscanf(field,"%lf",&bFactor)!=1)  {bFactor=0.0; }
      /**********************************/
      
      
      if(resNum!=old_resNum){ old_resNum=resNum; _sizeRes++; }
      Bead temp_bead;
      temp_bead.setCoordinates(Vector3d(x,y,z));
      temp_bead.setBetaFactor(bFactor);
      temp_bead.setPdbPosition(_size);
      temp_bead.setIndex(_size);
      temp_bead.setAtomType(type);
      temp_bead.setResName(resName);
      temp_bead.setResNum(resNum);
      temp_bead.setHet(ISHET);

      if(_beadList.check(temp_bead)){
	bool cacca=_beadList.check(temp_bead);
	//      	cout<<cacca<<endl;
	/* check capacity */
	if(_beads.capacity() <= _beads.size()) _beads.reserve(_beads.capacity() + expandSize);
	/* add record */
	_beads.push_back(temp_bead);	
	++_size;

      }//endif (isABead)
    }
  }
  if(_beads.size()>0){
    updateNeighbors();
    updateCovBonded();}
  return -1; //ritorna -1 se e' finito il file;
}//end function readFromPDBFile()


void Structure3d::addBead(Bead temp_bead){
  const size_t expandSize = 1000; // size of expansion of capacity
  if(_beadList.check(temp_bead)){
    //    bool cacca=_beadList.check(temp_bead);
    //      	cout<<cacca<<endl;
    /* check capacity */
    if(_beads.capacity() <= _beads.size()) _beads.reserve(_beads.capacity() + expandSize);
    /* add record */
    _beads.push_back(temp_bead);	
    ++_size;
  }//endif (isABead)
}


double Structure3d::getDistance(size_t i, size_t j) const{
  assert(i<_size && j < _size);
  return _beads[i].getCoordinates().distance(_beads[j].getCoordinates());
}



Vector3d Structure3d::getDistanceVector(size_t i, size_t j) {
  return _beads[j].getCoordinates() - _beads[i].getCoordinates();
}

double Structure3d::getDistanceSQ(size_t i, size_t j) const {
  assert(i<_size && j < _size);
  return _beads[i].getCoordinates().distanceSQ(_beads[j].getCoordinates());
}



void Structure3d::updateNeighbors() {
  
  _size = _beads.size();
  double cutoffSQ = _cutoff*_cutoff;
  
  for (size_t i = 0; i < _size; ++i) {
    for (size_t j = i+1; j < _size; ++j){
      if(_beads[j].getCoordinates().distanceSQ(_beads[i].getCoordinates())< cutoffSQ){
	_beads[i].addNeighbor(&_beads[j]);
	_beads[j].addNeighbor(&_beads[i]);
      }
    }
  }
}



void Structure3d::updateCovBonded() {

  _size=_beads.size();
  _beads[0].getAtomType();
  /* ### stiffest springs between bonded nodes ### */
  for(int i=0; i < _size-1; i++){
    if (_beads[i].getAtomType()=="P"){
      int j=i+1;
      if (_beads[j].getAtomType()=="C1'"){
	double dsq=_beads[j].getCoordinates().distanceSQ(_beads[i].getCoordinates());
	if (dsq < CHAIN_BREAK_DIST_SQ){
	  _beads[i].addCovBond(&_beads[j]);
	  _beads[j].addCovBond(&_beads[i]);
	}
      }//endif
    }//endif P
    else if (_beads[i].getAtomType()=="C1'"){
      int j=i+1;
      if (_beads[j].getAtomType()=="C2"){
	double dsq=_beads[j].getCoordinates().distanceSQ(_beads[i].getCoordinates());
	if (dsq < CHAIN_BREAK_DIST_SQ){
	  _beads[i].addCovBond(&_beads[j]);
	  _beads[j].addCovBond(&_beads[i]);
	}
      }//endif
      j=i+2;
      if(j>_size-1) continue;
      if (_beads[j].getAtomType()=="P"){
	double dsq=_beads[j].getCoordinates().distanceSQ(_beads[i].getCoordinates());
	if (dsq < CHAIN_BREAK_DIST_SQ){
	  _beads[i].addCovBond(&_beads[j]);
	  _beads[j].addCovBond(&_beads[i]);
	}
      }//endif
    }//endif C4'
  }//enddo I
}

bool Structure3d::isBond(int i, int j){
  for(int k=0;k<_beads[i].getNumCovBonded();++k){
    int n=(*(_beads[i].getCovBonded().at(k))).getIndex();
    if(n==j) return true;
  }
  return false;
}

// bool Structure3d::_isABead(char *a){

//   bool cacca=false;

//   if(_beadList.size()==0) return true;

//   stringstream ss;
//   string s;
//   ss<<a;
//   ss>>s;
  
//   for(int i=0;i<_beadList.size();++i)
//     {
//        string b;
//        b=_beadList.at(i);
//        if (s.compare(b)==0) cacca=true;
//        //       cout<<b<<endl;
//     }
//   return cacca;
// }


int Structure3d::centersFromPDBFile(const char* fname){
  cout<<"Structure3d:: reading centers of masses"<<endl;
  FILE* fp;
  if ((fp = fopen (fname, "r")) == NULL) {
    fprintf (stderr,"Could not open file %s.\n.Exiting.\n", fname);
    exit (1);    }
  int idum= centersFromPDBFile(fp);
  fclose(fp); 
  return idum;
}

int Structure3d::centersFromPDBFile(FILE *fp) {
  const size_t expandSize = 1000; // size of expansion of capacity

  int discard_alternative_position=0;
  int currentChainId = 0;
  char line[200], *a;
  char field[10];
  double x=0,y=0,z=0,bFactor;

  vector<Vector3d> sugar, base;


  // ### here I clear the beads to cancel the previous frame (if some)
  if(_beads.size()>0) _beads.clear();

  while (fgets(line, 200, fp)){
    if (strncmp (line, "ENDMDL", 6)==0){
      updateNeighbors();
      updateCovBonded();
      return 0;
    }
    /* check if line buffer begins with "ATOM" */
    if ( (!strncmp (line, "ATOM", 4))||(!strncmp (line, "HETATM", 6)) ){
      a = line +12; /* advance line buffer and check the atom's type */
      /* check if char at 16 is diff from ' ' or 'A' */
      discard_alternative_position=1;
      if (line[16]==' ') discard_alternative_position=0;
      if (line[16]=='A') discard_alternative_position=0;

      if (discard_alternative_position==1) continue;

      a = line +12; /* advance line buffer and check if
                                 we have a suitable atom */
      
      bool civa=false;
      int resNum;
      stringstream ss;
      string s, type;
      ss<<a;
      ss>>s;

      /* Phosphates */
      if (s.compare("P")==0){ civa=true; 
	type=s;
	/* Residue number */
	a=line+22;
	strncpy(field,a,4); field[4]='\0';
	if (sscanf(field,"%d",&resNum)!=1){fprintf(stderr,"Error. Could not read residue number in molecule.\nLine %s\n",line); exit(1);}
	/* Coordinates */
        a=line+30;
        strncpy(field,a,8);field[8]='\0';
        if (sscanf(field,"%lf",&x)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}
        a=line+38;
        strncpy(field,a,8);field[8]='\0';
        if (sscanf(field,"%lf",&y)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}
        a=line+46;
        strncpy(field,a,8);field[8]='\0';
        if (sscanf(field,"%lf",&z)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}
        /* B factor */
       	bFactor=0.0; 
      }
      
      /* sugars */ 
      bool issugar=false;
      if(s.compare("C2'")==0) issugar=true; 
      if(s.compare("C3'")==0) issugar=true; 
      if(s.compare("C4'")==0) issugar=true; 
      if(s.compare("O4'")==0) issugar=true; 
      if(s.compare("C1'")==0) issugar=true;
      
      if(issugar){
	/* store coordinates of sugar atoms */
	double tempx, tempy, tempz;

	a=line+30;	strncpy(field,a,8);field[8]='\0';
        if (sscanf(field,"%lf",&tempx)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}	

        a=line+38;        strncpy(field,a,8);field[8]='\0';
        if (sscanf(field,"%lf",&tempy)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}
	
        a=line+46;        strncpy(field,a,8);field[8]='\0';
        if (sscanf(field,"%lf",&tempz)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}
	
	sugar.push_back(Vector3d(tempx,tempy,tempz));
	
	/* if the ring is completely stored, compute the center */
	if(sugar.size()==5)
	  {  civa=true; bFactor=0.0; type="S";
	    x=0;y=0;z=0;
	    for (int i=0; i<5;++i){ 
	      x+=sugar.at(i).X;
	      y+=sugar.at(i).Y;
	      z+=sugar.at(i).Z;
	    }
	    x*=0.2;
	    y*=0.2;
	    z*=0.2;
	    sugar.clear();
	    /* Residue number */ a=line+22; strncpy(field,a,4); field[4]='\0';
	    if (sscanf(field,"%d",&resNum)!=1){fprintf(stderr,"Error. Could not read residue number in molecule.\nLine %s\n",line); exit(1);}      }
	
      }//end issugar
    

      /*  bases */
      
      bool isbase=false;
      if(s.compare("C5")==0) isbase=true;// x=0; y=0; z=0;}
      if(s.compare("C2")==0) isbase=true;

      if(isbase){
	/* store coordinates of base atoms */
	double tempx, tempy, tempz;
	
	a=line+30;	strncpy(field,a,8);field[8]='\0';
        if (sscanf(field,"%lf",&tempx)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}	
	
        a=line+38;        strncpy(field,a,8);field[8]='\0';
        if (sscanf(field,"%lf",&tempy)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}
	
        a=line+46;        strncpy(field,a,8);field[8]='\0';
        if (sscanf(field,"%lf",&tempz)!=1){fprintf(stderr,"Fatal error. Could not read bead coordinate in molecule.\nLine %s\n",line); exit(1);}
	
	base.push_back(Vector3d(tempx,tempy,tempz));
	
	/* if the ring is completely stored, compute the center */
	if(base.size()==2){
	  civa=true; bFactor=0.0; type="B";
	  x=0;y=0;z=0;
	  for (int i=0; i<2;++i){ 
	    x+=base.at(i).X;
	    y+=base.at(i).Y;
	    z+=base.at(i).Z;
	  }
	  x*=0.5;
	  y*=0.5;
	  z*=0.5;
	  base.clear();
	  /* Residue number */ a=line+22; strncpy(field,a,4); field[4]='\0';
	  if (sscanf(field,"%d",&resNum)!=1){fprintf(stderr,"Error. Could not read residue number in molecule.\nLine %s\n",line); exit(1);}      
	}
      }//end isbase
	
      
      /**********************************/
      if (civa) {   	
        /* check capacity */
        if(_beads.capacity() <= _beads.size()) _beads.reserve(_beads.capacity() + expandSize);
        /* add record */
        _beads.push_back(Bead());
        _beads.back().setCoordinates(Vector3d(x,y,z)); x=0; y=0; z=0;
        _beads.back().setBetaFactor(bFactor);
        _beads.back().setPdbPosition(_size);
        _beads.back().setIndex(_size);
	_beads.back().setAtomType(type);
	_beads.back().setResNum(resNum);
        ++_size;
	if (type=="B") ++_sizeRes;
      }//endif (civa)
    }//endif ATOM
  }//END WHILE
  if(_beads.size()>0){
  updateNeighbors();
  updateCovBonded();}
  return -1;
}




  int Structure3d::dumpPDBFile(const char* fname){
   char name[200];
   ofstream dumpfile;
   sprintf(name,"%s",fname);
   dumpfile.open(name);
   return dumpPDBFile(dumpfile); 
 }

int Structure3d::dumpPDBFile(std::ofstream &dumpfile){
  char line[200], *a;
  char field[10];


  // PDB FORMAT

//   COLUMNS        DATA TYPE     FIELD        DEFINITION

  char record[6];  //  1 2 3 4 5 6        Record name   "ATOM  "
  char serial[5];  //  7 8 9 10 11      Integer       serial       Atom serial number.
  char name[4];    // 13 14 15 16        Atom          name         Atom name.
  char altLoc[2];     // 17             Character     altLoc       Alternate location.
  char resName[3]; // 18 19 20        Residue name  resName      Residue name.
  char chainID[2];    // 22             Character     chainID      Chain identifier.
  char resSeq[4];  // 23 24 25 26        Integer       resSeq       Residue sequence.
  char iCode[2];      // 27             AChar         iCode       Code 4 insertion of res
  char X[7];  // 31 32 33 34 35 36 37 38        Real(8.3)     x     X in Angstroms
  char Y[7];  // 39 40 41 42 43 44 45 46        Real(8.3)     y     Y in Angstroms
  char Z[7];  // 47 48 49 50 51 52 53 54        Real(8.3)     z     Z in Angstroms
  char occupancy[6];  // 55 56 57 58 59 60        Real(6.2)    occupancy     Occupancy.
  char tempFactor[6]; // 61 62 63 64 65 66        Real(6.2)    tempFactor   B-factor.
  char element[2];    // 77 78   LString(2)  element    Element symbol, right-justified.
  char charge[2];     // 79 80        LString(2)   charge        Charge on the atom.


  dumpfile<<"REMARK    GENERATED FROM STRUCTURE3D"<<endl;
  dumpfile<<"TITLE"<<_title;
  dumpfile<<"MODEL"<<endl;

  for(int iline=0; iline<_beads.size(); ++iline){

    if(_beads[iline].getHet()) sprintf(record,"HETATM");
    else                       sprintf(record,"ATOM  ");

    sprintf(serial,"%5d",iline+1);

    sprintf(name," %s    ",_beads.at(iline).getAtomType().c_str());
    name[4]='\0';

    sprintf(altLoc," \0");
    
    //    sprintf(resName,"  G\0");
    sprintf(resName,"  %s   ",_beads.at(iline).getResName().c_str());
    resName[3]='\0';
    
    sprintf(chainID,"A\0");
    
    sprintf(resSeq,"%4d",_beads.at(iline).getResNum());

    sprintf(iCode," \0");

    Vector3d coords=_beads.at(iline).getCoordinates();
    sprintf(X,"%8.3lf",coords.X);
    sprintf(Y,"%8.3lf",coords.Y);
    sprintf(Z,"%8.3lf",coords.Z);


    dumpfile<<record
	    <<serial<<" "
	    <<name
	    <<altLoc
	    <<resName<<" "
	    <<chainID
	    <<resSeq
	    <<iCode<<"   "
	    <<X
	    <<Y
	    <<Z
	    <<endl;
    
  }//endloop
   
  dumpfile<<"TER"<<endl;
  dumpfile<<"ENDMDL"<<endl;

  return -1;
}


double Structure3d::fit(Structure3d *ref_struc, double *weight){
  double **str1, **str2;



   if ((*ref_struc).getSize()!=_size){
     cerr<<endl<<"Structure3d.fit:: FATAL ERROR -> the two structures to align doesn't have the same size!"<<endl;
     exit(1);   }

   str1=d2t(3,_size);
   str2=d2t(3,_size);

   for(int i=0;i<_size;++i){
     str1[0][i]=(*ref_struc).getBead(i).getCoordinates().X;
     str1[1][i]=(*ref_struc).getBead(i).getCoordinates().Y;
     str1[2][i]=(*ref_struc).getBead(i).getCoordinates().Z;    

     str2[0][i]=_beads[i].getCoordinates().X;
     str2[1][i]=_beads[i].getCoordinates().Y;
     str2[2][i]=_beads[i].getCoordinates().Z;
   }
   
   /* ______________ QUI FACCIO L'ALLINEAMENTO ____________ */
   //   double *weight=d1t(_size);
   //   for(int i=0; i<_size;++i) weight[i]=1.0;
   double rmsd = QCPAlignment(str1,str2,_size,weight);
   //   double rmsd=0;
   /* _____________________________________________________ */


   for(int i=0;i<_size;++i){
     double x=str2[0][i];
     double y=str2[1][i];
     double z=str2[2][i];
     _beads[i].setCoordinates(Vector3d(x,y,z));
   }
   
   free(str1);
   free(str2);
   
   return rmsd;

}
