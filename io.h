#ifndef IO_H_
#define IO_H_

#include <stdio.h>
#include <stdlib.h>

/**********/
static FILE *open_file_r (char *name){
  FILE *fp;
  if ((fp = fopen (name, "r")) == NULL)
    {
      printf ("Could not open file %s.\n.Exiting.\n", name);
      exit (1);
    }
  return (fp);
}
/**********/
static FILE *open_file_w (char *name){
  FILE *fp;
  if ((fp = fopen (name, "w")) == NULL)
    {
      printf ("Could not open file %s.\n.Exiting.\n", name);
      exit (1);
    }
  return (fp);
}
/**********/
static FILE *open_file_a (char *name){
  FILE *fp;
  if ((fp = fopen (name, "a")) == NULL)
    {
      printf ("Could not open file %s.\n.Exiting.\n", name);
      exit (1);
    }
  return (fp);
}
/**********/
static void close_file (FILE * fp){
  fclose (fp);
}

#endif /* IO_H_ */


/**********/

//// ############ functions to read input from pdb ############
//int get_mol_len_xyz(ifstream &ifile){
//   int mol_len=0;
//   char line[200];
//  //find out the length of the molecule
//   ifile>>mol_len;   
//   ifile.seekg(0, ios::beg);
//   //  ifile.getline(line,200);
//   //   cout<<line<<endl;
//   return 3*mol_len;
//} //END FUNCTION GET_MOL_LEN_XYZ
//
//int read_coords_xyz(ifstream &ifile, double *coordinates, int size){
//  double cacca;
//  int idum;
//  double dum;
//  char chdum;
//  char line[200];
//  //    cout<<"yeaahh"<<endl;
//  ifile>>idum;
//  //cout<<idum<<endl;
//  ifile.getline(line,200);
//  ifile.getline(line,200); 
//  //cout<<line<<endl;
//  // cout<<merda<<" ";
//  for(int i=0; i<size/3;i++){
//    ifile>>chdum>>coordinates[3*i]>>coordinates[3*i+1]>>coordinates[3*i+2];
//    //   cout<<chdum<<" "<<coordinates[3*i]<<endl;
//  }//ENDFOR
//
//  if(ifile.eof()) return 1;
//  else return 0;
//
//}//END FUNCTION READ_COORDS_XYZ
//#####################################################
