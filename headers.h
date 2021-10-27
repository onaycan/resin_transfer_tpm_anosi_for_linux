#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<iomanip>
#include<time.h>
#include<string>
#include <odb_API.h>
#include "csr_Matrix.h"

using namespace std;

//global definitions of the global functions
int node_counter(ifstream &_nodes_infile, int &_numberofnodes); 
int element_counter(ifstream &_nodes_infile, int &_numberofnodes); 
void transpose(double tens[3][3], double tens_T[3][3]); 
void matrix_product(double tens1[3][3],double tens2[3][3],double result[3][3]);
double trace(double tens[3][3]);
double cronecker(int i,int j);
double symmetry_tensor(int i,int j, int k, int l);
void du_copy(double du[],double du_prev[],ofstream& debugfile,int numberofnodes);
double eulernorm_residual(double fi[],int numberofnodes);
double error_indu(double du[],double du_prev[],ofstream& debugfile,int numberofnodes);
double det(double a[][3]);
void inverse(double tens[3][3], double inv_tens[3][3]);
void SetArrayZero(double _v[], int _m);


//functions of linear algebra 




class Node
{
private:
int nodelabel;
double *X_lamb; //undeformed coordinates
double *u_lamb; //total deformation
double *du_lamb; //total deformation
int *natural_bc;
double *natural_val; //total deformation
int *nonzero_natural_bc;
int *neumann_bc;
double *neumann_val;
double ns; //solid volume fraction in the deformed configuration 
double nS; //solid volume fraction in the reference configuration 
vector<int> neighbour_nodes;
double inpermeability; 
double layup;

public:
Node(){}  //default constructer
friend class Element; 
//friend functions 
void Push_Neighbour_In(int _n);
friend void K_allocater(Node _nodes[], csrM & _K, int** _Kloc, int** _eft, int _numberofnodes, int _numberofelements);
friend void apply_natural_bc(Node _nodes[],ifstream& _natural_bc,int _numberofnodes, double _u[], ofstream & _debugfile);
friend void apply_neumann_bc(double _force[], Node _nodes[], ifstream& _force_bc, int numberofnodes);
friend void apply_inpermeability_bc(Node _nodes[], ifstream& _inpermeability_bc, int _numberofnodes);
friend void set_layup_integration(Node _nodes[], Element _elements[], int _numberofnodes, int _numberofelements);
friend void node_allocater(Node _nodes[], ifstream &_nodes_infile, int &_numberofnodes);
friend void element_allocater(Node _nodes[], Element _elements[], ifstream &_elements_infile, int &_numberofnodes, int &_numberofelements);
friend void Kglobal_reducer(Node _nodes[], double** _Kglobal, double _fi[], ofstream& _debugfile, int _numberofnodes);
friend void Kglobal_ff_reducer(Node _nodes[], double** _Kglobal, double _fi[], ofstream& _debugfile, int _numberofnodes);
friend void csrKglobal_reducer(Node _nodes[], csrM & _K, double _fi[], ofstream& _debugfile, int _numberofnodes);
friend void csrKglobal_ff_reducer(Node _nodes[], csrM & _K, double _fi[], ofstream& _debugfile, int _numberofnodes);
friend void Nodes_neumann_update(Node nodes[], double du[], int numberofdatapoints);
	//gauss points friend functions
friend void Grad_at_Gauss(Node _nodes[],
	                  Element _elements[],
	                  double _ksi, double _eta, double _zeta,
                          double &_detj,
                          int _currentelement,
	                  double B[][32], //Strain displacement matrix 
                          double B_zero[][32],
		          double F[], //deformation gradient 
                          ofstream &_debugfile);
friend void Elasticity_at_Gauss(double D[][6],  
           		        double F[], //deformation gradient 
				double S[],
                                ofstream &_debugfile,
				double lambda);
friend void K_gauss_geometrical(Node _nodes[],
	                        Element _elements[],
	                        double _ksi, double _eta, double _zeta,
		        	double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
			        double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
				double F[],
				double S[],
			        double _Kgeo[][32],//the geometrical stiffness matrix to be passed 
			        double& detj,
                                int _currentelement,
                                ofstream &_debugfile);
friend void K_element(Node _nodes[],
	              Element _elements[],
                      int _currentelement,
                      ofstream &_debugfile,
                      double t);
friend void K_global(Node _nodes[],
	             Element _elements[],
	             int _numberofelements,
                     ofstream &_debugfile,
		     double t);
friend void gausselim(double rhs[],double** K,int numberofdatapoints,ofstream& outfile,Node nodes[], double du[]);
friend void incremental_force_initial(Node nodes[], Element elements[], double fi[], double numberofloadsteps, double l, ofstream& debugfile,int numberofnodes,int numberofelements, double incremental_force[]);
friend void update_fi(Node nodes[], Element elements[], double fi[], double numberofloadsteps, double l, ofstream& debugfile,int numberofnodes,int numberofelements, double incremental_force[]);
friend void result_writer(Node _nodes[], Element _elements[],
                          int _numberofnodes, int _numberofelements);
friend void odb_model_writer(Node _nodes[], Element _elements[],
                              int _numberofnodes, int _numberofelements,char *argv[],odb_Odb& myodb);

friend void odb_result_writer(Node _nodes[], Element _elements[],
                              int _numberofnodes, int _numberofelements,char *argv[],odb_Odb& odb, double l, double numberofloadsteps);
friend void fullincrementalforcevector(Node nodes[],double incrementalforcevector[], double du[],
				double** Kglobal, ofstream& debugfile, int numberofnodes, double fi[]);
friend void csrfullincrementalforcevector(Node nodes[],double incrementalforcevector[], double du[],
		                	  csrM K, ofstream& debugfile, int numberofnodes, double fi[]);

friend void initialize(Node _nodes[], double _numberofloadsteps, ofstream& _debugfile, int _numberofnodes, double _du[], double _u[], double _incremental_force[], double _force[], double _l);
friend void Nodes_disp_update(Node _nodes[], int _numberofnodes);
friend void Nodes_disp_downdate(Node _nodes[], int _numberofnodes);

// TPM STIFFNESS FUNCTIONS 

friend void calc_K_Wlambdau_at_gauss(Node _nodes[],
    	                             Element _elements[],
	                             double _ksi, double _eta, double _zeta,
         	                     double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                             double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                     double F[],
	                             double S[],
	                             double _K_Wlambdau[][32],//the geometrical stiffness matrix to be passed 
	                             double& detj,
                                     int _currentelement,
                                     ofstream &_debugfile,
	 		             double _t,
                                     double _coeff,
                                     double current_lambda);


friend void calc_K_Wlambda_at_gauss(Node _nodes[],
	                            Element _elements[],
	                            double _ksi, double _eta, double _zeta,
         	                    double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                            double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                    double F[],
	                            double S[],
        	                    double _K_Wlambda[][32],//the geometrical stiffness matrix to be passed 
	                            double& detj,
                                    int _currentelement,
                                    ofstream &_debugfile,
	 		            double _t,
                                    double _coeff);


friend void calc_K_M1uJ_at_gauss(Node _nodes[],
	                         Element _elements[],
	                         double _ksi, double _eta, double _zeta,
                	         double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                         double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                 double F[],
	                         double S[],
	                         double _K_M1uJ[][32],//the geometrical stiffness matrix to be passed 
	                         double& detj,
                                 int _currentelement,
                                 ofstream &_debugfile,
			         double _t,
                                 double _coeff);

friend void calc_K_M1uGeo_at_gauss(Node _nodes[],
	                           Element _elements[],
	                           double _ksi, double _eta, double _zeta,
         	                   double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                           double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                   double F[],
	                           double S[],
	                           double _K_M1uGeo[][32],//the geometrical stiffness matrix to be passed 
	                           double& detj,
                                   int _currentelement,
                                   ofstream &_debugfile,
	    		           double _t,
				   double _coeff);

friend void Pore_Elasticity_at_Gauss(double D_pore[][6],  
	    	                     double F[], //deformation gradient 
			             ofstream &_debugfile,
			             double lambda);

friend void calc_K_M1uMat_at_gauss(Node _nodes[],
	                           Element _elements[],
	                           double _ksi, double _eta, double _zeta,
         	                   double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                           double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                   double F[],
	                           double S[],
	                           double _K_M1uMat[][32],//the geometrical stiffness matrix to be passed 
	                           double& detj,
                                   int _currentelement,
                                   ofstream &_debugfile,
	    		           double _t,
				   double _coeff);

friend void calc_K_M1uN_at_gauss(Node _nodes[],
	                         Element _elements[],
	                         double _ksi, double _eta, double _zeta,
         	                 double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                         double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                 double F[],
	                         double S[],
	                         double _K_M1uN[][32],//the geometrical stiffness matrix to be passed 
	                         double& detj,
                                 int _currentelement,
                                 ofstream &_debugfile,
	    		         double _t,
                                 double _coeff);


friend void calc_K_M2uJ_at_gauss(Node _nodes[],
	                         Element _elements[],
	                         double _ksi, double _eta, double _zeta,
         	                 double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                         double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                 double F[],
	                         double S[],
	                         double _K_M2uJ[][32],//the geometrical stiffness matrix to be passed 
	                         double& detj,
                                 int _currentelement,
                                 ofstream &_debugfile,
     	    	                 double _t,
			         double _coeff,
				 double _current_layup);

friend void calc_K_M2uMat_at_gauss(Node _nodes[],
	                           Element _elements[],
	                           double _ksi, double _eta, double _zeta,
         	                   double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                           double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                   double F[],
	                           double S[],
	                           double _K_M2uMat[][32],//the geometrical stiffness matrix to be passed 
	                           double& detj,
                                   int _currentelement,
                                   ofstream &_debugfile,
     		                   double _t,
				   double _current_layup);


friend void calc_K_M2lambda_at_gauss(Node _nodes[],
	                             Element _elements[],
	                             double _ksi, double _eta, double _zeta,
         	                     double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                             double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                     double F[],
	                             double S[],
	                             double _K_M2lambda[][32],//the geometrical stiffness matrix to be passed 
	                             double& detj,
                                     int _currentelement,
                                     ofstream &_debugfile,
     		                     double _t,
				     double _current_layup);


//Internal forces caused by the conservation of mass equations 
friend void calc_internal_M1(Node _nodes[],
	                     Element _elements[],
	                     double _ksi, double _eta, double _zeta,
               	             double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                     double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                             double F[],
	                     double S[],
	                     double _internal_M1[32],//internalforce of the conservation of mass 1 
	                     double& detj,
                             int _currentelement,
                             ofstream &_debugfile,
	                     double _t,
			     double _coeff);

friend void calc_internal_M2(Node _nodes[],
	                     Element _elements[],
	                     double _ksi, double _eta, double _zeta,
               	             double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                     double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                             double F[],
	                     double S[],
	                     double _internal_M1[32],//internalforce of the conservation of mass 1 
	                     double& detj,
                             int _currentelement,
                             ofstream &_debugfile,
	                     double _t,
		             double _current_layup);



};

class Element
{
private:
int elementlabel;
int* conname;
int *connid;
Node* elementnodes;
double Kelem[32][32];
double preS[6];
double ns[8];
double v_fluid[24];
double v_fs[24];
double layup[8];

public:
Element(){};
//friend functions 
friend class Node; 
	//gauss points friend functions 
friend void node_allocater(Node _nodes[], ifstream &_nodes_infile, int &_numberofnodes);
friend void element_allocater(Node _nodes[], Element _elements[], ifstream &_elements_infile, int &_numberofnodes, int &_numberofelements);
friend void set_layup_integration(Node _nodes[], Element _elements[], int _numberofnodes, int _numberofelements);
friend void eftable(Element _elements[], int** _eft,ofstream& _debugfile,int _numberofelements);
friend void Grad_at_Gauss(Node _nodes[],
	                  Element _elements[],
	                  double _ksi, double _eta, double _zeta,
                          double &_detj,
                          int _currentelement,
	                  double B[][32], //Strain displacement matrix 
                          double B_zero[][32],
		          double F[], //deformation gradient 
                          ofstream &_debugfile);
friend void Elasticity_at_Gauss(double D[][6],  
           		        double F[], //deformation gradient 
				double S[],
                                ofstream &_debugfile,
				double lambda);
friend void K_gauss_geometrical(Node _nodes[],
	                        Element _elements[],
	                        double _ksi, double _eta, double _zeta,
		        	double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
			        double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
				double F[],
				double S[],
			        double _Kgeo[][32],//the geometrical stiffness matrix to be passed 
			        double& detj,
                                int _currentelement,
                                ofstream &_debugfile);
friend void K_element(Node _nodes[],
	              Element _elements[],
                      int _currentelement,
                      ofstream &_debugfile,
                      double t);
friend void K_global(Node _nodes[],
	             Element _elements[],
	             int _numberofelements,
                     ofstream &_debugfile,
                     double t);
friend void K_assembly(Element _elements[], 
                       double** _Kglobal, 
                       int** _eft,ofstream& _debugfile,
                       int _numberofnodes, 
                       int _numberofelements);
friend void csrK_assembly(Element _elements[], csrM & _K,  int** _eft, int _numberofnodes, int _numberofelements);
friend void incremental_force_initial(Node nodes[], Element elements[], double fi[], double numberofloadsteps, double l, ofstream& debugfile,int numberofnodes,int numberofelements, double incremental_force[]);
friend void update_fi(Node nodes[], Element elements[], double fi[], double numberofloadsteps, double l, ofstream& debugfile,int numberofnodes,int numberofelements, double incremental_force[]);
friend void result_writer(Node _nodes[], Element _elements[],
                          int _numberofnodes, int _numberofelements);
friend void odb_model_writer(Node _nodes[], Element _elements[],
                              int _numberofnodes, int _numberofelements,char *argv[],odb_Odb& myodb);
friend void odb_result_writer(Node _nodes[], Element _elements[],
                              int _numberofnodes, int _numberofelements,char *argv[],odb_Odb& odb, double l, double numberofloadsteps);

// TPM STIFFNESS FUNCTIONS 

friend void calc_K_Wlambdau_at_gauss(Node _nodes[],
    	                             Element _elements[],
	                             double _ksi, double _eta, double _zeta,
         	                     double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                             double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                     double F[],
	                             double S[],
	                             double _K_Wlambdau[][32],//the geometrical stiffness matrix to be passed 
	                             double& detj,
                                     int _currentelement,
                                     ofstream &_debugfile,
	 		             double _t,
                                     double _coeff,
                                     double current_lambda);

friend void calc_K_Wlambda_at_gauss(Node _nodes[],
	                            Element _elements[],
	                            double _ksi, double _eta, double _zeta,
         	                    double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                            double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                    double F[],
	                            double S[],
        	                    double _K_Wlambda[][32],//the geometrical stiffness matrix to be passed 
	                            double& detj,
                                    int _currentelement,
                                    ofstream &_debugfile,
	 		            double _t,
                                     double _coeff);



friend void calc_K_M1uJ_at_gauss(Node _nodes[],
	                         Element _elements[],
	                         double _ksi, double _eta, double _zeta,
                	         double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                         double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                 double F[],
	                         double S[],
	                         double _K_M1uJ[][32],//the geometrical stiffness matrix to be passed 
	                         double& detj,
                                 int _currentelement,
                                 ofstream &_debugfile,
			         double _t,
                                 double _coeff);

friend void calc_K_M1uGeo_at_gauss(Node _nodes[],
	                           Element _elements[],
	                           double _ksi, double _eta, double _zeta,
         	                   double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                           double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                   double F[],
	                           double S[],
	                           double _K_M1uGeo[][32],//the geometrical stiffness matrix to be passed 
	                           double& detj,
                                   int _currentelement,
                                   ofstream &_debugfile,
	    		           double _t,
				   double _coeff);


friend void Pore_Elasticity_at_Gauss(double D_pore[][6],  
	    	                     double F[], //deformation gradient 
			             ofstream &_debugfile,
			             double lambda);

friend void calc_K_M1uMat_at_gauss(Node _nodes[],
	                           Element _elements[],
	                           double _ksi, double _eta, double _zeta,
         	                   double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                           double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                   double F[],
	                           double S[],
	                           double _K_M1uMat[][32],//the geometrical stiffness matrix to be passed 
	                           double& detj,
                                   int _currentelement,
                                   ofstream &_debugfile,
	    		           double _t,
                                   double _coeff);

friend void calc_K_M1uN_at_gauss(Node _nodes[],
	                         Element _elements[],
	                         double _ksi, double _eta, double _zeta,
         	                 double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                         double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                 double F[],
	                         double S[],
	                         double _K_M1uN[][32],//the geometrical stiffness matrix to be passed 
	                         double& detj,
                                 int _currentelement,
                                 ofstream &_debugfile,
	    		         double _t,
                                 double _coeff);


friend void calc_K_M2uJ_at_gauss(Node _nodes[],
	                         Element _elements[],
	                         double _ksi, double _eta, double _zeta,
         	                 double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                         double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                 double F[],
	                         double S[],
	                         double _K_M2uJ[][32],//the geometrical stiffness matrix to be passed 
	                         double& detj,
                                 int _currentelement,
                                 ofstream &_debugfile,
     	    	                 double _t,
			         double _coeff,
				 double _current_layup);

friend void calc_K_M2uMat_at_gauss(Node _nodes[],
	                           Element _elements[],
	                           double _ksi, double _eta, double _zeta,
         	                   double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                           double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                   double F[],
	                           double S[],
	                           double _K_M2uMat[][32],//the geometrical stiffness matrix to be passed 
	                           double& detj,
                                   int _currentelement,
                                   ofstream &_debugfile,
     		                   double _t,
				   double _current_layup);

friend void calc_K_M2lambda_at_gauss(Node _nodes[],
	                             Element _elements[],
	                             double _ksi, double _eta, double _zeta,
         	                     double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                             double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                                     double F[],
	                             double S[],
	                             double _K_M2lambda[][32],//the geometrical stiffness matrix to be passed 
	                             double& detj,
                                     int _currentelement,
                                     ofstream &_debugfile,
     		                     double _t,
				     double _current_layup);

//Internal forces caused by the conservation of mass equations 
friend void calc_internal_M1(Node _nodes[],
	                     Element _elements[],
	                     double _ksi, double _eta, double _zeta,
               	             double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                     double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                             double F[],
	                     double S[],
	                     double _internal_M1[32],//internalforce of the conservation of mass 1 
	                     double& detj,
                             int _currentelement,
                             ofstream &_debugfile,
	                     double _t,
			     double _coeff);

friend void calc_internal_M2(Node _nodes[],
	                     Element _elements[],
	                     double _ksi, double _eta, double _zeta,
               	             double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	                     double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                             double F[],
	                     double S[],
	                     double _internal_M1[32],//internalforce of the conservation of mass 1 
	                     double& detj,
                             int _currentelement,
                             ofstream &_debugfile,
	                     double _t,
			     double _current_layup);


};

