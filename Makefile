LAPACK=-I${MKL_HOME}/include/ -L${MKL_LIB} -liomp5 -lmkl_intel_lp64 -lpthread  -Wl,--start-group  -lmkl_intel_thread -lmkl_sequential -lmkl_core -Wl,--end-group -lm -L${MKLPATH} ${MKLPATH}/libmkl_solver_lp64.a

all: ribogm pca compdist genconf distC2 ang_C1p_P enm_new countn fast_rwsip fast_bhatt fit entrcomp

#program to do the Elastic Network Model with different beads, different cutoffs, different models
ribogm: RiboGM.o Matrix.o Structure3d.o ElasticNet.o CovMatrix1d.o
	g++ -O3 ${LAPACK} RiboGM.o Matrix.o Structure3d.o ElasticNet.o CovMatrix1d.o -o ribogm.x
RiboGM.o: RiboGM.C
	g++ -O3 -c RiboGM.C

# program for principal component analysis
pca: pca.o Matrix.o PrincipalComp.o Structure3d.o CovMatrix1d.o
	g++ -O3 pca.o Matrix.o PrincipalComp.o Structure3d.o CovMatrix1d.o ${LAPACK} -o pca.x
pca.o: pca.C
	g++ -O3 -c pca.C

#compute the fluctuations of the distances in a trajectory (generated from genconf) (probably obsolete)
compdist: compdist.o Structure3d.o
	g++ -O3 compdist.o Structure3d.o -o compdist.x
compdist.o: compdist.C
	g++ -O3 -c compdist.C

#generate gaussian trajectories from one ENM
genconf: genConf.o Structure3d.o Matrix.o PrincipalComp.o CovMatrix1d.o
	g++ -O3 genConf.o Structure3d.o PrincipalComp.o Matrix.o CovMatrix1d.o ${LAPACK} -o genConf.x
genConf.o: genConf.C
	g++ -O3 -c genConf.C

# compute the fluctuation of distances between neighbouring C2
distC2:  distC2.o Structure3d.o Matrix.o PrincipalComp.o ElasticNet.o CovMatrix1d.o
	g++ -O3 -fopenmp distC2.o Structure3d.o PrincipalComp.o Matrix.o ElasticNet.o CovMatrix1d.o ${LAPACK} -o distC2.x
distC2.o: distC2.C		
	g++ -O3 -fopenmp -c distC2.C

# compute the fluctuation of angles between neighbouring C2
ang_C1p_P:  ang_C1p_P.o Structure3d.o Matrix.o PrincipalComp.o ElasticNet.o CovMatrix1d.o
	g++ -O3 -fopenmp ang_C1p_P.o Structure3d.o PrincipalComp.o Matrix.o ElasticNet.o CovMatrix1d.o ${LAPACK} -o ang_C1p_P.x
ang_C1p_P.o: ang_C1p_P.C		
	g++ -O3 -fopenmp -c ang_C1p_P.C

# elastic network model program: fast, safe, beautiful
enm_new:  enm_new.o Structure3d.o Matrix.o PrincipalComp.o ElasticNet.o CovMatrix1d.o
	g++ -O3 -fopenmp enm_new.o Structure3d.o PrincipalComp.o Matrix.o ElasticNet.o CovMatrix1d.o ${LAPACK} -o enm_new.x
enm_new.o: enm_new.C		
	g++ -O3 -fopenmp -c enm_new.C

#count the number of neighbours of each bead, for a given cutoff
countn:  countn.o Structure3d.o Matrix.o PrincipalComp.o ElasticNet.o CovMatrix1d.o
	g++ -O3 countn.o Structure3d.o PrincipalComp.o Matrix.o ElasticNet.o CovMatrix1d.o ${LAPACK} -o countn.x
countn.o: countn.C		
	g++ -O3 -c countn.C

# compute (fast and sound) the RWSIP with a given PCA output (hopefully)
fast_rwsip:  fast_rwsip.o Structure3d.o Matrix.o PrincipalComp.o ElasticNet.o CovMatrix1d.o
	g++ -O3 fast_rwsip.o Structure3d.o PrincipalComp.o Matrix.o ElasticNet.o CovMatrix1d.o ${LAPACK} -o fast_rwsip.x
fast_rwsip.o: fast_rwsip.C		
	g++ -O3 -c fast_rwsip.C

# compute (fast and sound) the RWSIP with a given PCA output (hopefully)
fast_bhatt:  fast_Bhatt.o Structure3d.o Matrix.o PrincipalComp.o ElasticNet.o CovMatrix1d.o
	g++ -O3 fast_Bhatt.o Structure3d.o PrincipalComp.o Matrix.o ElasticNet.o CovMatrix1d.o ${LAPACK} -o fast_bhatt.x
fast_bhatt.o: fast_Bhatt.C		
	g++ -O3 -c fast_Bhatt.C

#align different structures (.pdb trajectory)
fit: fit.o Structure3d.o
	g++ -O3 fit.o Structure3d.o -o fit.x
fit.o: fit.C qcprot.cpp
	g++ -O3 -c fit.C 

# compute the entropy
entrcomp: EntrComp.o ElasticNet.o Structure3d.o Matrix.o CovMatrix1d.o
	g++ -O3 ${LAPACK} EntrComp.o Matrix.o Structure3d.o ElasticNet.o CovMatrix1d.o -o entrcomp.x
EntrComp.o: EntrComp.C
	g++ -O3 -c EntrComp.C

## some classes 
ElasticNet.o: ElasticNet.cc lapack_matrix_routines_wrapper_v3.cc
	g++ -O3 -fopenmp -c ElasticNet.cc

Structure3d.o: Structure3d.cpp qcprot.cpp
	g++ -O3 -c Structure3d.cpp

Matrix.o: Matrix.cc lapack_matrix_routines_wrapper_v3.cc
	g++ -O3 -fopenmp -c Matrix.cc

CovMatrix1d.o: CovMatrix1d.cc lapack_matrix_routines_wrapper_v3.cc
	g++ -O3 -fopenmp -c CovMatrix1d.cc

PrincipalComp.o: PrincipalComp.cc
	g++ -O3 -fopenmp -c PrincipalComp.cc

NormalModes.o: NormalModes.cpp
	g++ -O3 -c NormalModes.cpp

clean: 
	rm -rf fit.o ElasticNet.o Structure3d.o Matrix.o CovMatrix1d.o PrincipalComp.o generateConfigs.o dist.o RiboGM.o EntrComp.o pca.o distC2.o Entrcomp.o fast_rwsip.o fast_bhatt.o enm_new.o countn.o compdist.o genConf.o ang_C1p_P.o