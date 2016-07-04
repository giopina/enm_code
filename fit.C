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
#include "Matrix.h" //definition of the class CMatrix
//#include "PrincipalComp.h"
#include "Structure3d.h"

using namespace std;

void help_display(){
  cout<<"Help:"<<endl<<
    "-f --> input file  (required) format .pdb"<<endl<<
    "-o --> output file (required) format .pdb"<<endl<<
    "-s --> reference structure file format .pdb (default reference=first frame of input trajectory)"<<endl<<
    "-beads --> beads to consider (default = ALL )"<<endl<<
    "-lasthr --> fit only on the first/last 3 beads"<<endl<<
    "-par --> parameters file (DEFAULT = none)"<<endl<<
     "-h help"<<endl;;
  return;
}

int main(int argc, char **argv){//PROGRAM MAIN
  char file_name[200], outf_name[200], fitf_name[200];
  char par_name[200];
  bool par_has_name=false;
  int n_comp=0;
  int ntop=-1;
  bool LASTHR=false;
  bool FITONFIRST=true;
  //  bool PROJECT=false;

  vector<string> bead_list;
  // list.push_back("C1'");
  // list.push_back("C2");
  // list.push_back("P");


  cout<<" --- Program FIT --- "<<endl<<endl;
  int i=1;
  int nopar=0;
  while(i<argc){
    if(!strncmp(argv[i],"-f",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(file_name,argv[i]);
      nopar+=1;
      cout<<"will analize the file "<<file_name<<endl;
    }
    else if(!strncmp(argv[i],"-o",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(outf_name,argv[i]);
      nopar+=1;
      cout<<"will write the trajectory on the file "<<outf_name<<endl;
    }
    else if(!strncmp(argv[i],"-s",2)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(fitf_name,argv[i]);
      nopar+=1;
      FITONFIRST=false;
      cout<<"will fit on the structure on the file "<<fitf_name<<endl;
    }
    else if(!strncmp(argv[i],"-lasthr",2)){
      //      i=i+1;
      //      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      //      sprintf(outf_name,argv[i]);
      //      nopar+=1;
      LASTHR=true;
      cout<<"will fit on first+last 3 residues"<<endl;
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
//    else if(!strncmp(argv[i],"-beads",6)){
//      bead_list.clear();
//      i=i+1;
//      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; help_display(); return 0;}
//      while(strncmp(argv[i],"[",1)!=0){cout<<"ERROR: in parameter -beads"<<endl; help_display(); return 0;};
//      i=i+1;
//      while(strncmp(argv[i],"]",1)!=0){
//	char name[4];
//	sprintf(name,argv[i]);
//	bead_list.push_back(name);
//	i=i+1;
//      }
//      if(bead_list.size()==0) {
//	cout<<"WARNING: no beads declared!"<<endl;;
//      }
//    }    
    else if(!strncmp(argv[i],"-par",4)){
      i=i+1;
      if(i>argc-1){cout<<"ERROR: empty parameter"<<endl; return 0;}
      sprintf(par_name,argv[i]);
      cout<<"will read parameters from the file "<<par_name<<endl;
      par_has_name=true;
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
  if(nopar<2){cout<<"ERROR: insufficient number of parameters"<<endl;
    help_display;
    return 0;
  }

  
  Structure3d structure(0);

  cout<<"Beads to consider are: ";
  //  cout<<bead_list.size();
  for(int i=0; i<bead_list.size();++i){
    cout<<bead_list.at(i)<<" ";}
  cout<<endl<<endl;



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
  BeadList lista(bead_list); // setto i numeri dei residui di output
  lista.setResNum(res_list);
  structure.setBeadList(lista);

  // +++ OPEN the INPUTfile +++
  FILE* fp;
  if ((fp = fopen (file_name, "r")) == NULL) {
    fprintf (stderr,"Could not open file %s.\n.Exiting.\n", file_name);
    exit (1);    }


  // +++ OPEN the OUTPUTfile +++
  ofstream dumpfile;
  dumpfile.open(outf_name);

  int error=0;
  int Nsteps=0;

  Structure3d refStructure(0);
  refStructure.setBeadList(lista);
  
  error=structure.readFromPDBFile(fp);
  int N=structure.getSize();
  //  int size=3*_N;
  cout<<"#taglia : "<<N<<endl;
  cout<<"#error : "<<error<<endl;

  //putting the correct weights
  double* weight=d1t(N);
  if(LASTHR){  
  weight[0]=1.0;
  weight[1]=1.0;
  weight[2]=1.0;
  //  weight[3]=1.0;
  weight[N-1]=1.0;
  weight[N-2]=1.0;
  weight[N-3]=1.0;
  }
  else{ for(int i=0; i<N;++i) weight[i]=1.0; }

  if(FITONFIRST)  refStructure=structure;
  else refStructure.readFromPDBFile(fitf_name);
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  while(error==0){                                     // * * *   
    ++Nsteps;					       // * * * 
    if(Nsteps%1000==0) cout<<"Step "<<Nsteps<<endl;    // * * * 
    //cout<<"Step "<<Nsteps<<"\t";                       // * * * 
    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    double rmsd=structure.fit(&refStructure,weight);
    //    cout<<"RMSD="<<rmsd<<"\t";
    //    cout<<"size="<<structure.getSize()<<"\t";
    //    cout<<"TITLE: "<<structure.getTitle()<<endl;
    structure.dumpPDBFile(dumpfile);

    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    //    cout<<"read"<<endl;
    error=structure.readFromPDBFile(fp);  	       // * * *
    //    cout<<error<<endl;
  }//end while  			 	       // * * *
  //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  fclose(fp);
  dumpfile.close();
  free_d1t(weight);

  return 0;

}//END PROGRAM MAIN
