#include "headers.h"
#include <odb_API.h>
#include <vector>
#include <ctime>
//Mkl libraries
#include "mkl.h" 
#include "mkl_lapack.h" 
#include "mkl_blas.h" 
#include "mkl_service.h"
#include "mkl_spblas.h"
#include "mkl_pardiso.h"
#include "mkl_dss.h"

using namespace std;

int ABQmain(int argc,char *argv[])
{


//values to be changed as input 
double numberofloadsteps=5.0;
double ltwocurrentfi=10000.0;
double ltwoinitialfi=10000.0;

ifstream nodes_infile(argv[1]);
ifstream nodes_infile_alloc(argv[1]);
int numberofnodes;
Node *nodes; 

ifstream elements_infile(argv[1]);
ifstream elements_infile_alloc(argv[1]);
int numberofelements;
Element* elements;

//files of boundary conditions 
ifstream natural_bc("natural_bc.inp");
ifstream neumann_bc("neumann_bc.inp");
ifstream inpermeability_bc("inpermeability_bc.inp");
//debugfile 
ofstream debugfile("debugfile.deb");


//Global declerations, thought initially as dense matrices 
int **eft;
int** Kloc;
double* force;
double* force_constant;
double* df;
double* incremental_force;
double *u;
double* du;
double* du_prev;
double *fi;

//counting 
node_counter(nodes_infile,numberofnodes);
element_counter(elements_infile,numberofelements);


//allocation and initialization of the parameters 
nodes=new Node[numberofnodes];
elements=new Element[numberofelements];

node_allocater(nodes,nodes_infile_alloc,numberofnodes);
element_allocater(nodes, elements, elements_infile_alloc, numberofnodes, numberofelements);

//allocation of the global vectors 
Kloc=new int*[4*numberofnodes];
eft=new int*[32];
force=new double[4*numberofnodes];
df=new double[4*numberofnodes];
u=new double[4*numberofnodes];
du=new double[4*numberofnodes];
du_prev=new double[4*numberofnodes];
fi=new double[4*numberofnodes];
incremental_force=new double[4*numberofnodes];
for(int i=0;i<32;i++)
eft[i]=new int[numberofelements];

for(int i=0;i<4*numberofnodes;i++)
{

Kloc[i]=new int[4*numberofnodes];
for(int j=0;j<4*numberofnodes;j++)
Kloc[i][j]=0;
force[i]=0.0;
df[i]=0.0;
incremental_force[i]=0.0;
du_prev[i]=0.0;
u[i]=0.0;
du[i]=0.0;
}






// Generating the odb 
std::string odb_name;
odb_name.assign(argv[1]);

odb_name.erase(odb_name.end()-3,odb_name.end());
odb_name.push_back('o');
odb_name.push_back('d');
odb_name.push_back('b');


char* odb_name_inchar=new char[odb_name.length()];

    odb_String name("simpleModel");
    odb_String analysisTitle("ODB created with C++ ODB API");
    odb_String description("example illustrating C++ ODB API");
    for(int i=0;i<odb_name.length();i++)
    odb_name_inchar[i]=odb_name[i];

    odb_String path(odb_name_inchar);
    
    odb_Odb &myodb = Odb(name,
                       analysisTitle,
                       description,
                       path);


odb_model_writer(nodes,elements,numberofnodes,numberofelements,argv,myodb);


apply_natural_bc(nodes,natural_bc,numberofnodes,u,debugfile);
apply_neumann_bc(force,nodes,neumann_bc,numberofnodes);
apply_inpermeability_bc(nodes, inpermeability_bc, numberofnodes);
set_layup_integration(nodes, elements, numberofnodes, numberofelements);


eftable(elements,eft,debugfile,numberofelements);







unsigned begin;
unsigned end;
double elapsed_secs;



cout<<"Start Allocation..."<<endl;
begin=clock();
csrM K; 
//K.CreateOfHollow(4*numberofnodes);
K_allocater(nodes,K,Kloc,eft,numberofnodes,numberofelements);
end=clock();
cout<<"Allocation time:"<<(end-begin)<<endl;


for(double l=0;l<numberofloadsteps;l++)
{

double sometime=1; //steady state inposition 

//double sometime=0.0;
int numberofinneriterations=0;
ltwocurrentfi=100000.0;
ltwoinitialfi=100000.0;
double errormeasure;
initialize(nodes,numberofloadsteps,debugfile,numberofnodes,du,u,incremental_force,force,sometime);

cout<<'\t'<<"frame "<<l<<endl;
debugfile<<'\t'<<"frame "<<l<<endl;

int domore=1;
Nodes_disp_update(nodes,numberofnodes);
//incremental_force_initial(nodes,elements,fi,numberofloadsteps,sometime,debugfile,numberofnodes,numberofelements,incremental_force);
do
{

Nodes_disp_downdate(nodes,numberofnodes);

K_global(nodes,elements,numberofelements,debugfile,sometime);


//cout<<"Start Assembly..."<<endl;
//begin=clock();
csrK_assembly(elements,K,eft,numberofnodes,numberofelements);
//end=clock();
//cout<<"Assembly time:"<<(end-begin)<<endl;



Nodes_disp_update(nodes,numberofnodes);


update_fi(nodes,elements,fi,numberofloadsteps,sometime,debugfile,numberofnodes,numberofelements,incremental_force);


//cout<<"Start First reduction..."<<endl;
//begin=clock();
csrKglobal_reducer(nodes,K,fi,debugfile,numberofnodes);
//end=clock();
//cout<<"First reduction time:"<<end-begin<<endl;


/*
if(numberofinneriterations==0 && l==0)
csrfullincrementalforcevector(nodes,incremental_force,du,K,debugfile,numberofnodes,fi);
*/


//if(numberofinneriterations==0.0)
//du_copy(du,du_prev,debugfile,numberofnodes);


//cout<<"Start second Reduction..."<<endl;
//begin=clock();
csrKglobal_ff_reducer(nodes,K,fi,debugfile,numberofnodes);
//end=clock();
//cout<<"Second reduction time:"<<end-begin<<endl;


//cout<<"Start Solving..."<<endl;
//begin=clock();
Pardiso_UnsymSolver(K,du,fi,4*numberofnodes,'Y');
//end=clock();
//cout<<"Solver time:"<<end-begin<<endl;

if(numberofinneriterations==0.0)
du_copy(du,du_prev,debugfile,numberofnodes);

Nodes_neumann_update(nodes, du, 4*numberofnodes); 

Nodes_disp_update(nodes,numberofnodes);


//if(numberofinneriterations>0.0)
errormeasure=error_indu(du,du_prev,debugfile,numberofnodes);
//else 
//errormeasure=1.0;

debugfile<<'\t'<<"iteration"<<numberofinneriterations<<'\t'<<errormeasure<<endl; 
cout<<'\t'<<"iteration"<<numberofinneriterations<<'\t'<<errormeasure<<endl; 
if(errormeasure<0.01 && numberofinneriterations>=5.0)
domore=0;
if(errormeasure>0.01 && numberofinneriterations>=25.0)
domore=0;
numberofinneriterations++;
K.SetZero();

}
while(domore==1);

if(l==numberofloadsteps-1)
Nodes_disp_update(nodes,numberofnodes);


odb_result_writer(nodes,elements,numberofnodes,numberofelements,argv,myodb,l,numberofloadsteps);


numberofinneriterations=0;
}


//result_writer(nodes,elements,numberofnodes,numberofelements);
//odb_result_writer(nodes,elements,numberofnodes,numberofelements,argv,myodb);

myodb.save();
myodb.close();

//cout<<"The file name you have given: "<<argv[1]<<endl;
return 0;
};
