#include <vector>

double correlation(std::vector<double> data1, std::vector<double> data2){
  double corr=0.0;
  int ndata1=data1.size();
  if(!(ndata1==data2.size())){std::cerr<<"ERROR: correlation->lunghezza dati diverse. Non posso correlarle. "<<ndata1<<" "<<data2.size()<<std::endl;}
  else{
    double sum_X=0.0;
    double sum_Y=0.0;
    double sum_X2=0.0;
    double sum_Y2=0.0;
    double sum_XY=0.0;      
    for(int i=0;i<ndata1;++i){
      sum_X+=data1[i];
      sum_X2+=data1[i]*data1[i];
      sum_Y+=data2[i];
      sum_Y2+=data2[i]*data2[i];
      sum_XY+=data1[i]*data2[i];
    }//enddo i
    corr=(sum_XY-(sum_X*sum_Y)/ndata1)/(sqrt((sum_X2-(sum_X*sum_X)/ndata1)*(sum_Y2-(sum_Y*sum_Y)/ndata1)));
  }
  return corr;
}

double kendall_tau(std::vector<double> data1, std::vector<double> data2){
  int ndata=data1.size();
  double kend=0.0;
  if(!(ndata==data2.size())){cout<<"ERROR: lunghezza dati data1 diversa da fluttuazioni. Non posso correlarle. "<<ndata<<" "<<data2.size()<<endl;}
  else{
    double C=0;
    double D=0;
    double C_n_2=ndata*(ndata-1)/2.;      
    for(int i=0;i<ndata;++i){
      for(int j=i+1;j<ndata;++j){
	if((data2[i]-data2[j])*(data1[i]-data1[j])>0.0)
	  C+=1;
	else
	  if((data2[i]-data2[j])*(data1[i]-data1[j])<0.0)
	    D+=1;
      }//enddo j
    }//enddo i    
    kend=(C-D)/C_n_2;
  }
  return kend;
}
  
