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


void csrM::CsrAddCsr(csrM L, csrM R, double alpha)
{

iM.push_back(0); 
//can: loop over the rows of the left matrix (in the same way around the rows of the second matrix)
    for(int l=1;l<L.iM.size();l++)
    {
        //can:first insert the (l-1)th rowth elements of the left matrix into the possible solution 
     vM.insert(vM.end(),L.vM.begin()+L.iM[l-1],L.vM.begin()+L.iM[l]); 
     jM.insert(jM.end(),L.jM.begin()+L.iM[l-1],L.jM.begin()+L.iM[l]);
     int search=jM.size()-L.iM[l]+L.iM[l-1];
     for(int r=R.iM[l-1];r<R.iM[l];r++)//can:elements in the (l-1)'th row of the second matrix 
     {
         int search_end=jM.size();
         int found=0;
         while(search<search_end && found==0)
         {
         if(R.jM[r]<jM[search])
         {
             vM.insert(vM.begin()+search,alpha*R.vM[r]); 
             jM.insert(jM.begin()+search,R.jM[r]);
             found=1;
         }
         else if(R.jM[r]==jM[search])
         {
             vM[search]=vM[search]+alpha*R.vM[r];
             found=1; 
         }
         else
             search++; 
         };
         if(found==0)
         {
             vM.push_back(alpha*R.vM[r]);
             jM.push_back(R.jM[r]); 
         }

     }

     iM.push_back(jM.size());
    }

}


void csrM::CsrMultCsr(csrM L, csrM R, double _sparseprec)
{



iM.push_back(0); 
//can: loop over the rows of the left matrix
    for(int rL=1;rL<L.iM.size();rL++)
    {   
        //cout<<"pr0cessing row:"<<rL-1<<endl;
        int row_counter=0;
        //can: loop over the number of entries in row rL-1 
        for(int r=L.iM[rL-1]; r<L.iM[rL]; r++)
        {
            for(int c=R.iM[L.jM[r]];c<R.iM[L.jM[r]+1];c++)//can: each c is a row(loc) of the second matrix
            {
                double currentval=L.vM[r]*R.vM[c];
                if(fabs(currentval)>_sparseprec)
                {
                if(row_counter==0)
                {
                vM.push_back(currentval);
                jM.push_back(R.jM[c]);
                }
                else
                {
                    int ord=1;
                    int current_location;
                    int inserted=0;
                    int added=0; 
                    while(ord<=row_counter && added==0)
                    {
                        if(R.jM[c]<jM[jM.size()-ord])
                        {
                        current_location=jM.size()-ord;
                        inserted=1;
                        added=0; 
                        }
                        if (R.jM[c]==jM[jM.size()-ord])
                        {
                        current_location=jM.size()-ord;
                        inserted=0;
                        added=1; 
                        }
                    ord++;
                    };
                    if(inserted==1)
                        {
                        vM.insert(vM.begin()+current_location,currentval);
                        jM.insert(jM.begin()+current_location,R.jM[c]);
                        }
                    else if(added==1)   
                        {
                        vM[current_location]=vM[current_location]+currentval;
                        row_counter--;
                        }
                    if(added==0 && inserted==0)
                        {
                        vM.push_back(currentval);
                        jM.push_back(R.jM[c]);
                        }
                }
                row_counter++;
                }
            }
        }
        
        iM.push_back(jM.size());
    }

}


void csrM::writeCsrM()
{

for(int i=0;i<vM.size();i++)
cout<<vM[i]<<'\t';
cout<<endl; 

for(int i=0;i<jM.size();i++)
cout<<jM[i]<<'\t';
cout<<endl; 

for(int i=0;i<iM.size();i++)
cout<<iM[i]<<'\t';
cout<<endl; 

}


void csrM::writeCsrM(ofstream& outfile)
{

for(int i=0;i<vM.size();i++)
outfile<<vM[i]<<'\t';
outfile<<endl; 

for(int i=0;i<jM.size();i++)
outfile<<jM[i]<<'\t';
outfile<<endl; 

for(int i=0;i<iM.size();i++)
outfile<<iM[i]<<'\t';
outfile<<endl; 

}


void csrM::writeCsrM_arrayform(ofstream& outfile)
{
for(int i=0;i<iM.size()-1;i++)
{
double row_sum=0.0;
for(int j=0;j<iM.size()-1;j++)
{
int zero=0;
for(int k=iM[i];k<iM[i+1];k++)
if(j==jM[k])
{
//row_sum=row_sum+vM[k];
outfile<<vM[k]<<setw(20);
zero=1;
}
if(zero==0)
outfile<<0.0<<setw(20);
}
//outfile<<row_sum<<endl;
outfile<<endl;
}

}

void csrM::writeCsrM_arrayform()
{
for(int i=0;i<iM.size()-1;i++)
{
double row_sum=0.0;
for(int j=0;j<iM.size()-1;j++)
{
int zero=0;
for(int k=iM[i];k<iM[i+1];k++)
if(j==jM[k])
{
//row_sum=row_sum+vM[k];
cout<<vM[k]<<'\t';
zero=1;
}
if(zero==0)
cout<<0.0<<'\t';
}
//outfile<<row_sum<<endl;
cout<<endl;
}

}




void csrM::Normalize_to_stochastic(ofstream& outfile)
{
for(int i=0;i<iM.size()-1;i++)
{
double row_sum=0.0;
for(int j=0;j<iM.size()-1;j++)
for(int k=iM[i];k<iM[i+1];k++)
if(j==jM[k])
row_sum=row_sum+vM[k];


if(row_sum==0.0)
cout<<"WARNING: the row sum is found to be zero in normalize to stochastic"<<endl;


for(int j=0;j<iM.size()-1;j++)
for(int k=iM[i];k<iM[i+1];k++)
if(j==jM[k])
vM[k]=vM[k]/row_sum;

}

}






void csrM::writeCsrM_row(int row, ofstream& outfile)
{

    for(int i=iM[row];i<iM[row+1];i++)
        outfile<<jM[i]<<endl;


}

void csrM::ConvertToIdentity(int _m)
{

KillCsrM();
vM.resize(_m);
jM.resize(_m); 
iM.resize(_m+1);

for(int i=0;i<_m;i++)
{
vM[i]=1.0; 
jM[i]=i;
iM[i]=i;
}
iM[_m]=_m;

}


void csrM::CreateOfHollow(int _m)
{
KillCsrM();
vM.resize(_m);
jM.resize(_m); 
iM.resize(_m+1);

for(int i=0;i<_m;i++)
{
vM[i]=0.0; 
jM[i]=i;
iM[i]=i;
}
iM[_m]=_m;
}

void csrM::SortCsrM(double _sparseprec)
{
csrM I; 
I.ConvertToIdentity(iM.size()-1); 
csrM sorted; 
sorted.CsrMultCsr(*this,I,_sparseprec);
I.KillCsrM();
this->KillCsrM();
this->AsignCsrM(sorted);
sorted.KillCsrM();
}


void csrM::AsignCsrM(csrM right)
{

for(int i=0;i<right.vM.size();i++)
vM.push_back(right.vM[i]); 

for(int i=0;i<right.jM.size();i++)
jM.push_back(right.jM[i]); 

for(int i=0;i<right.iM.size();i++)
iM.push_back(right.iM[i]); 

}


void csrM::KillCsrM()
{
iM.clear();
vM.clear();
jM.clear();
}

void csrM::MinusCsrM()
{
for(int i=0;i<vM.size();i++)
vM[i]=-vM[i];
}

void csrM::MaxNormCsrM(double &_norm)
{
_norm=0.0; 
for(int i=0;i<vM.size();i++)
if(fabs(vM[i])>_norm)
_norm=fabs(vM[i]);
}



void csrM::MaxNormCsrM_local(double &_norm, int& vMloc, int iter)
{
_norm=0.0; 
if(iter==0)
{
for(int i=0;i<vM.size();i++)
if(fabs(vM[i])>_norm)
_norm=fabs(vM[i]);
}
else
_norm=fabs(vM[vMloc]);
}


void csrM::NeumannInverter(csrM _M, double _tol, char complement, char killhim,
                           double &_sparsity_loss, double _sparseprec)
{
csrM I; 
csrM ThisPast; 
csrM This_min_ThisPast; 
int nrows=_M.iM.size()-1;
I.ConvertToIdentity(nrows);
this->ConvertToIdentity(nrows);

//can: if the complement of the matrix i.e. (I-_M) is 
//can: given as input, then the process below is unnecessary...
if(complement=='N')
{
csrM Imin_M; 
//_M.MinusCsrM();
Imin_M.CsrAddCsr(I,_M,-1.0);
_M.KillCsrM();
_M.AsignCsrM(Imin_M);
Imin_M.KillCsrM();
}


double norm_initial=0.1;
double norm_next=0.1; 
double norm=0.1; 
int itr=0;
int vMloc;
while(norm>_tol && itr<100)
{
//norm_past=norm_next; 
csrM MMpast;
MMpast.CsrMultCsr(_M,*this,_sparseprec); 
this->KillCsrM();
this->CsrAddCsr(MMpast,I,1.0); 
MMpast.KillCsrM();
if(itr>=1)
{
    This_min_ThisPast.CsrAddCsr(*this,ThisPast,-1.0); 
    ThisPast.KillCsrM();
    This_min_ThisPast.MaxNormCsrM(norm);
    This_min_ThisPast.KillCsrM();
}
ThisPast.AsignCsrM(*this);
if(itr%10==0)
{
cout<<"iteration: "<<itr<<endl; 
cout<<"maximum norm of the last power: "<<norm<<endl;
}
itr++; 
};
ThisPast.KillCsrM();
I.KillCsrM();
_sparsity_loss=(iM[iM.size()-1])/(_M.iM[_M.iM.size()-1]);
if(killhim=='Y')
_M.KillCsrM();
}


void csrM::SparseSolverInverter(csrM _M, double _tol, char complement, char killhim,
                                double &_sparsity_loss, double _sparseprec)
{
int m=_M.iM.size()-1;
double* current_rhs=new double[m];
double* current_lhs=new double[m];
this->CreateOfHollow(m);

for(int i=0;i<m;i++)
{
SetArrayZero(current_rhs,m);
SetArrayZero(current_lhs,m);
current_rhs[i]=1;
    Pardiso_UnsymSolver(_M,current_lhs,current_rhs,m,'N'); 
for(int j=0;j<m;j++)
if(fabs(current_lhs[j])>_sparseprec)
this->Push_Value_incsrM(current_lhs[j],j,i); 
double percent_finished=(double(i)+1)*100.0/double(m);
cout<<"%"<<percent_finished<<": percent is inverted"<<'\r';
}
_sparsity_loss=(iM[iM.size()-1])/(_M.iM[_M.iM.size()-1]);
if(killhim=='Y')
_M.KillCsrM();
}

void csrM::Pardiso_ExpInverter(csrM _K, double _tol, char complement, char killhim,
                               double &_sparsity_loss, double _sparseprec)
{
char comment='N';
int m=_K.iM.size()-1;
double* _rhs=new double[m];
double* _lhs=new double[m];
csrM thisbutcsc; 


//can: if the complement of the matrix i.e. (I-_M) is 
//can: given as input, then the process below is unnecessary...
//can: In fact, for pardiso this process is a must: 
if(complement=='N')
{
csrM I; 
I.ConvertToIdentity(m);
this->CsrAddCsr(I,_K,-1.0);
_K.KillCsrM();
_K.AsignCsrM(*this);
this->KillCsrM();
}

int n=_K.iM.size()-1; 
//cout<<"n:"<<n<<endl;
double* a=new double[_K.vM.size()];
//cout<<"a:"<<endl;
for(int i=0;i<_K.vM.size();i++)
{
a[i]=_K.vM[i];
//cout<<a[i]<<'\t';
}
//cout<<endl;
//cout<<"ia:"<<endl;
int*ia=new int[_K.iM.size()];
for(int i=0;i<_K.iM.size();i++)
{
ia[i]=_K.iM[i]+1;
//cout<<ia[i]<<'\t';
}
//cout<<endl;
//cout<<"ja:"<<endl;
int* ja=new int[_K.jM.size()];
for(int i=0;i<_K.jM.size();i++)
{
ja[i]=_K.jM[i]+1;
//cout<<ja[i]<<'\t';
}
//cout<<endl;


int mtype = 11; /* Real unsymmetric matrix */
/* RHS and solution vectors.*/
int nrhs = 1; /* Number of right hand sides. */
/* Internal solver memory pointer pt, */
/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
/* or void *pt[64] should be OK on both architectures */
void *pt[64];
/* Pardiso control parameters.*/
int iparm[64];
int maxfct, mnum, phase, error, msglvl;
/* Auxiliary variables. */
int i;
double ddum; /* Double dummy*/
int idum; /* Integer dummy.*/
/* --------------------------------------------------------------------*/
/* .. Setup Pardiso control parameters.*/
/* --------------------------------------------------------------------*/
for (i = 0; i < 64; i++) {
iparm[i] = 0;
}
iparm[0] = 1; /* No solver default*/
iparm[1] = 2; /* Fill-in reordering from METIS */
/* Numbers of processors, value of MKL_NUM_THREADS */
iparm[2] = mkl_get_max_threads();
iparm[3] = 0; /* No iterative-direct algorithm */
iparm[4] = 0; /* No user fill-in reducing permutation */
iparm[5] = 0; /* Write solution into x */
iparm[6] = 16; /* Default logical fortran unit number for output */
iparm[7] = 2; /* Max numbers of iterative refinement steps */
iparm[8] = 0; /* Not in use*/
iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
iparm[11] = 0; /* Not in use*/
iparm[12] = 0; /* Not in use*/
iparm[13] = 0; /* Output: Number of perturbed pivots */
iparm[14] = 0; /* Not in use*/
iparm[15] = 0; /* Not in use*/
iparm[16] = 0; /* Not in use*/
iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
iparm[18] = -1; /* Output: Mflops for LU factorization */
iparm[19] = 0; /* Output: Numbers of CG Iterations */
maxfct = 1; /* Maximum number of numerical factorizations. */
mnum = 1; /* Which factorization to use. */
msglvl = 0; /* Don't print statistical information in file */
error = 0; /* Initialize error flag */
/* --------------------------------------------------------------------*/
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* --------------------------------------------------------------------*/
for (i = 0; i < 64; i++) {
pt[i] = 0;
}
/* --------------------------------------------------------------------*/
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* --------------------------------------------------------------------*/
phase = 11;
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, a, ia, ja, &idum, &nrhs,
iparm, &msglvl, &ddum, &ddum, &error);
if (error != 0) {
printf("\nERROR during symbolic factorization: %d", error);
exit(1);
}
if(comment=='Y')
{
printf("\nReordering completed ... ");
printf("\nNumber of nonzeros in factors = %d", iparm[17]);
printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
}
/* --------------------------------------------------------------------*/
/* .. Numerical factorization.*/
/* --------------------------------------------------------------------*/
phase = 22;
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, a, ia, ja, &idum, &nrhs,
iparm, &msglvl, &ddum, &ddum, &error);
if (error != 0) {
printf("\nERROR during numerical factorization: %d", error);
exit(2);
}
if(comment=='Y')
printf("\nFactorization completed ... ");
/* --------------------------------------------------------------------*/
/* .. Back substitution and iterative refinement. */
/* --------------------------------------------------------------------*/
phase = 33;
iparm[7] = 2; /* Max numbers of iterative refinement steps. */
/* Set right hand side to one.*/

/*
for (i = 0; i < n; i++) {
b[i] = 1;
}
*/ 
this->iM.resize(m+1);
for(int i=0;i<m+1;i++)
this->iM[i]=0.0;

for(int i=0;i<m;i++)
{
SetArrayZero(_rhs,m);
SetArrayZero(_lhs,m);
_rhs[i]=1;
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, a, ia, ja, &idum, &nrhs,
iparm, &msglvl, _rhs, _lhs, &error);
if (error != 0) {
printf("\nERROR during solution: %d", error);
exit(3);
}
if(comment=='Y')
printf("\nSolve completed ... "); 

//can: everyting is added if it is more the sparse precision
//-can: in compresses sparse column format into the *this
int visiting_count=0; 
for(int j=0;j<m;j++)
if(fabs(_lhs[j])>_sparseprec)
{
visiting_count++;
thisbutcsc.vM.push_back(_lhs[j]);
thisbutcsc.jM.push_back(j);
if(visiting_count==1)
thisbutcsc.iM.push_back(thisbutcsc.jM.size()-1);
this->iM[j+1]=this->iM[j+1]+1;
}
if(visiting_count==0)
thisbutcsc.iM.push_back(thisbutcsc.jM.size());



double percent_finished=(double(i)+1)*100.0/double(m);
cout<<"%"<<percent_finished<<": percent is inverted "<<'\r';
}

/* --------------------------------------------------------------------*/
/* .. Termination and release of memory. */
/* --------------------------------------------------------------------*/
phase = -1; /* Release internal memory. */
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, &ddum, ia, ja, &idum, &nrhs,
iparm, &msglvl, &ddum, &ddum, &error);
delete []a;
delete []ia;
delete []ja;
cout<<endl;

thisbutcsc.iM.push_back(thisbutcsc.jM.size());
//this->AsignCsrM(thisbutcsc); 
//this->ConvertFromCSC(thisbutcsc); 
cout<<"start of conversion of csc to csr"<<endl; 
//can: this method is a specific way of doing it 
for(int i=0;i<m;i++)
this->iM[i+1]=this->iM[i+1]+this->iM[i];
this->vM.resize(thisbutcsc.vM.size()); 
for(int i=0;i<thisbutcsc.vM.size();i++)
this->vM[i]=0.0; 
this->jM.resize(thisbutcsc.jM.size()); 
cout<<"start of conversion of csc to csr resized"<<endl; 
for(int i=0;i<m;i++)
for(int j=thisbutcsc.iM[i];j<thisbutcsc.iM[i+1];j++)
{
//row of this: thisbutcsc.jM[j]
//column of this: i 
//value of this: thisbutcsc.vM[j]
if(fabs(thisbutcsc.vM[j])>_sparseprec)
{
int searchofempty=0;
int search=0; 
while(searchofempty==0)
{
if(vM[iM[thisbutcsc.jM[j]]+search]==0.0)
{
vM[iM[thisbutcsc.jM[j]]+search]=thisbutcsc.vM[j];
jM[iM[thisbutcsc.jM[j]]+search]=i;
searchofempty=1; 
}
else 
search++;
};
}
}
cout<<"end of conversion of csc to csr"<<endl; 
thisbutcsc.KillCsrM();
_sparsity_loss=(iM[iM.size()-1])/(_K.iM[_K.iM.size()-1]);
if(killhim=='Y')
_K.KillCsrM();





}



void csrM::ConvertFromCSC(csrM thisbutcsc)
{
cout<<"conversion from csc to csr is (succesfully) started"<<endl; 
int m=thisbutcsc.iM.size()-1; 
    this->CreateOfHollow(m); 
    for(int i=0;i<m;i++)
        for(int j=thisbutcsc.iM[i];j<thisbutcsc.iM[i+1];j++)
            this->Push_Value_incsrM(thisbutcsc.vM[j],thisbutcsc.jM[j],i);
cout<<"conversion from csc to csr is (succesfully) completed"<<endl; 
/*
int info; 
int job[7]; 
job[0]=0; //can: from csr to csc 
job[1]=0; //can: zero based csr 
job[2]=0; //can: zero based csc 
job[5]=1;
job[6]=1; //can: all the output files are filled; 
int m=thisbutcsc.iM.size()-1; 
double* acsr=new double[thisbutcsc.vM.size()];; 
int* ja=new int[thisbutcsc.jM.size()];; 
int* ia=new int[thisbutcsc.iM.size()];; 
double* acsc=new double[thisbutcsc.vM.size()];
int* ja1=new int[thisbutcsc.jM.size()];
int* ia1=new int[thisbutcsc.iM.size()];

for(int i=0;i<thisbutcsc.vM.size();i++)
acsr[i]=thisbutcsc.vM[i]; 
for(int i=0;i<thisbutcsc.jM.size();i++)
ja[i]=thisbutcsc.jM[i]; 
for(int i=0;i<thisbutcsc.iM.size();i++)
ia[i]=thisbutcsc.iM[i]; 

cout<<"conversion from csc to csr is (succesfully) started"<<endl; 
mkl_dcsrcsc(job, &m, acsr, ja, ia, acsc, ja1, ia1, &info);
cout<<"conversion from csc to csr is (succesfully) completed"<<endl; 


delete []acsr;
delete []ja;
delete []ia;

this->vM.resize(thisbutcsc.vM.size());
this->jM.resize(thisbutcsc.jM.size());
this->iM.resize(thisbutcsc.iM.size());

for(int i=0;i<thisbutcsc.vM.size();i++)
this->vM[i]=acsc[i]; 
for(int i=0;i<thisbutcsc.jM.size();i++)
this->jM[i]=ja1[i]; 
for(int i=0;i<thisbutcsc.iM.size();i++)
this->iM[i]=ia1[i]; 

delete []acsc;
delete []ja1;
delete []ia1;
*/ 
}


void csrM::Push_Value_incsrM(double val, int row, int column)
{
int i=iM[row];
int inserted=0;
int added=0;


while(i<iM[row+1] && inserted==0 && added==0)
{
if(column<jM[i])
{
jM.insert(jM.begin()+i,column);
vM.insert(vM.begin()+i,val);
inserted=1;
//cout<<"inserted, row, column: "<<val<<endl;
//cout<<row<<endl; 
//cout<<column<<endl;
}
else if(column==jM[i])
{
vM[i]=vM[i]+val;
added=1;
//cout<<"added, row, column: "<<val<<endl;
//cout<<row<<endl; 
//cout<<column<<endl;
}
i++;
};

if(inserted==0 && added==0)
{
jM.insert(jM.begin()+i,column);
vM.insert(vM.begin()+i,val);
//cout<<"pushed, row, column: "<<val<<endl;
//cout<<row<<endl; 
//cout<<column<<endl;
}

if(added==0)
for(int j=row+1;j<iM.size();j++)
iM[j]=iM[j]+1;

//cout<<"added:"<<added<<endl;
//writeCsrM();
}


void csrM::Insert_Value_incsrM(double val, int row, int column)
{
int i=iM[row];
int inserted=0;


while(i<iM[row+1] && inserted==0)
{
if(jM[i]==column)
{
vM[i]=vM[i]+val;
inserted=1;
}
else
i++; 
};

if(inserted==0)
cout<<"Error in inversion in csrM!"<<endl;

}



/*
void csrM::Insert_Matrix_incsrM(ublas::c_matrix<double,3,3> val, int row, int column)
{
int i=iM[row];
int inserted=0;


while(i<iM[row+1] && inserted==0)
{
if(jM[i]==column)
{
vM[i]=vM[i]+val(0,0);
vM[i+1]=vM[i+1]+val(0,1);
vM[i+2]=vM[i+2]+val(0,2);

int second_row=iM[row+1]-iM[row]+i;
vM[second_row]=vM[second_row]+val(1,0);
vM[second_row+1]=vM[second_row+1]+val(1,1);
vM[second_row+2]=vM[second_row+2]+val(1,2);

int third_row=iM[row+2]-iM[row]+i;
vM[third_row]=vM[third_row]+val(2,0);
vM[third_row+1]=vM[third_row+1]+val(2,1);
vM[third_row+2]=vM[third_row+2]+val(2,2);

inserted=1;
}
else
i=i+3; 
};

if(inserted==0)
cout<<"Error in inversion in csrM!"<<endl;

}
*/

void csrM::Insert_Matrix_incsrM(double Kelem_node[4][4], int row, int column)
{
int i=iM[row];
int inserted=0;


while(i<iM[row+1] && inserted==0)
{
if(jM[i]==column)
{

vM[i]=vM[i]+Kelem_node[0][0];
vM[i+1]=vM[i+1]+Kelem_node[0][1];
vM[i+2]=vM[i+2]+Kelem_node[0][2];
vM[i+3]=vM[i+3]+Kelem_node[0][3];


int second_row=iM[row+1]-iM[row]+i;
vM[second_row]=vM[second_row]+Kelem_node[1][0];
vM[second_row+1]=vM[second_row+1]+Kelem_node[1][1];
vM[second_row+2]=vM[second_row+2]+Kelem_node[1][2];
vM[second_row+3]=vM[second_row+3]+Kelem_node[1][3];


int third_row=iM[row+2]-iM[row]+i;
vM[third_row]=vM[third_row]+Kelem_node[2][0];
vM[third_row+1]=vM[third_row+1]+Kelem_node[2][1];
vM[third_row+2]=vM[third_row+2]+Kelem_node[2][2];
vM[third_row+3]=vM[third_row+3]+Kelem_node[2][3];


int fourth_row=iM[row+3]-iM[row]+i;
vM[fourth_row]=vM[fourth_row]+Kelem_node[3][0];
vM[fourth_row+1]=vM[fourth_row+1]+Kelem_node[3][1];
vM[fourth_row+2]=vM[fourth_row+2]+Kelem_node[3][2];
vM[fourth_row+3]=vM[fourth_row+3]+Kelem_node[3][3];


inserted=1;
}
else
i=i+4; 
};

if(inserted==0)
cout<<"Error in insersion of Matrix in csrM!"<<endl;

};

void csrM::Shrink_Column(int column)
{
csrM Ipart; 
Ipart.ConvertToIdentity(iM.size()-1); 
Ipart.vM.erase(Ipart.vM.begin()+column);
Ipart.jM.erase(Ipart.jM.begin()+column);
for(int c=column+1;c<Ipart.iM.size();c++)
Ipart.iM[c]=c-1; 
//cout<<"identity part:"<<endl;
//Ipart.writeCsrM_arrayform();

double _sparseprec=0.0;
csrM red; 
red.CsrMultCsr(*this,Ipart,_sparseprec); 
red.Push_Value_incsrM(1.0,column,column); 
KillCsrM(); 
AsignCsrM(red);
red.KillCsrM();
}


void csrM::Shrink_Columns(vector<int> columns)
{
for(int i=0;i<columns.size();i++)
for(int c=0;c<iM.size();c++)
if(columns[i]==iM[c])
vM[c]=0.0;
/*
for(int i=0;i<columns.size();i++)
this->Push_Value_incsrM(1.0,columns[i],columns[i]); 
*/
/*
csrM Ipart; 
Ipart.ConvertToIdentity(iM.size()-1); 
for(int i=0;i<columns.size();i++)
{
Ipart.vM[columns[i]]=0.0;
//Ipart.vM.erase(Ipart.vM.begin()+columns[i]);
//Ipart.jM.erase(Ipart.jM.begin()+columns[i]);
//for(int c=columns[i]+1;c<Ipart.iM.size();c++)
//Ipart.iM[c]=c-1; 
}
//cout<<"identity part:"<<endl;
//Ipart.writeCsrM_arrayform();

double _sparseprec=-1.0;
csrM red; 
red.CsrMultCsr(*this,Ipart,_sparseprec); 

for(int i=0;i<columns.size();i++)
red.Push_Value_incsrM(1.0,columns[i],columns[i]); 
KillCsrM(); 
AsignCsrM(red);
red.KillCsrM();
*/
}

void csrM::Shrink_Row(int row)
{
csrM Ipart; 
Ipart.ConvertToIdentity(iM.size()-1); 
Ipart.vM.erase(Ipart.vM.begin()+row);
Ipart.jM.erase(Ipart.jM.begin()+row);
for(int r=row+1;r<Ipart.iM.size();r++)
Ipart.iM[r]=r-1; 
//cout<<"identity part:"<<endl;
//Ipart.writeCsrM_arrayform();

double _sparseprec=0.0;
csrM red; 
red.CsrMultCsr(Ipart,*this,_sparseprec); 
red.Push_Value_incsrM(1.0,row,row); 
KillCsrM(); 
AsignCsrM(red);
red.KillCsrM();
}

void csrM::Shrink_Rows(vector<int> rows)
{
for(int i=0;i<rows.size();i++)
for(int j=iM[rows[i]];j<iM[rows[i]+1];j++)
vM[j]=0.0;
for(int i=0;i<rows.size();i++)
this->Insert_Value_incsrM(1.0,rows[i],rows[i]); 

/*
csrM Ipart; 
Ipart.ConvertToIdentity(iM.size()-1); 

for(int i=0;i<rows.size();i++)
{
Ipart.vM[rows[i]]=0.0;
//Ipart.vM.erase(Ipart.vM.begin()+rows[i]);
//Ipart.jM.erase(Ipart.jM.begin()+rows[i]);
//for(int r=rows[i]+1;r<Ipart.iM.size();r++)
//Ipart.iM[r]=r-1; 
}

//cout<<"identity part:"<<endl;
//Ipart.writeCsrM_arrayform();

double _sparseprec=-1.0;
csrM red; 
red.CsrMultCsr(Ipart,*this,_sparseprec); 
for(int i=0;i<rows.size();i++)
red.Push_Value_incsrM(1.0,rows[i],rows[i]); 
KillCsrM(); 
AsignCsrM(red);
red.KillCsrM();
*/
}




void csrM::CsrMultVec(double _v[], double _result[], int _size)
{
for(int i=0;i<_size;i++)
{
double value=0.0;
for(int j=iM[i];j<iM[i+1];j++)
value=value+vM[j]*_v[jM[j]];
_result[i]=value;
}
}


void csrM::SetZero()
{
for(int i=0;i<vM.size();i++)
vM[i]=0.0;
}







void csrM::ScalarToTensorial(csrM Mscalar)
{
int size_scalar=Mscalar.iM.size()-1;
int size_tensorial=size_scalar*3;
double zero=0.0;

iM.resize(size_tensorial+1);
//jM.resize(Mscalar.jM.size()*9); 
//vM.resize(Mscalar.vM.size()*9);

for(int i=0;i<size_scalar;i++)
for(int rtenso=0;rtenso<3;rtenso++)
{
for(int j=Mscalar.iM[i];j<Mscalar.iM[i+1];j++)
{
    jM.push_back(Mscalar.jM[j]*3);
    jM.push_back(Mscalar.jM[j]*3+1);
    jM.push_back(Mscalar.jM[j]*3+2);
    vM.push_back(zero);
    vM.push_back(zero);
    vM.push_back(zero);
}
iM[3*i+rtenso+1]=jM.size();
}
}



int Pardiso_UnsymSolver(csrM _K, double _lhs[], double _rhs[], int _size, char comment)
{

int n=_K.iM.size()-1; 
//cout<<"n:"<<n<<endl;
double* a=new double[_K.vM.size()];
//cout<<"a:"<<endl;
for(int i=0;i<_K.vM.size();i++)
{
a[i]=_K.vM[i];
//cout<<a[i]<<'\t';
}
//cout<<endl;
//cout<<"ia:"<<endl;
int*ia=new int[_K.iM.size()];
for(int i=0;i<_K.iM.size();i++)
{
ia[i]=_K.iM[i]+1;
//cout<<ia[i]<<'\t';
}
//cout<<endl;
//cout<<"ja:"<<endl;
int* ja=new int[_K.jM.size()];
for(int i=0;i<_K.jM.size();i++)
{
ja[i]=_K.jM[i]+1;
//cout<<ja[i]<<'\t';
}
//cout<<endl;

/*
int n = 8;
int ia[ 9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
int ja[18] = { 1, 3, 6, 7,
2, 3, 5,
3, 8,
4, 7,
5, 6, 7,
6, 8,
7,
8 };
double a[18] = { 7.0, 1.0, 2.0, 7.0,
-4.0, 8.0, 2.0,
1.0, 5.0,
7.0, 9.0,
5.0, 1.0, 5.0,
-1.0, 5.0,
11.0,
5.0 };
double  b[8], x[8];
*/
int mtype = 11; /* Real unsymmetric matrix */
/* RHS and solution vectors.*/
int nrhs = 1; /* Number of right hand sides. */
/* Internal solver memory pointer pt, */
/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
/* or void *pt[64] should be OK on both architectures */
void *pt[64];
/* Pardiso control parameters.*/
int iparm[64];
int maxfct, mnum, phase, error, msglvl;
/* Auxiliary variables. */
int i;
double ddum; /* Double dummy*/
int idum; /* Integer dummy.*/
/* --------------------------------------------------------------------*/
/* .. Setup Pardiso control parameters.*/
/* --------------------------------------------------------------------*/
for (i = 0; i < 64; i++) {
iparm[i] = 0;
}
iparm[0] = 1; /* No solver default*/
iparm[1] = 2; /* Fill-in reordering from METIS */
/* Numbers of processors, value of MKL_NUM_THREADS */
iparm[2] = mkl_get_max_threads();
iparm[3] = 0; /* No iterative-direct algorithm */
iparm[4] = 0; /* No user fill-in reducing permutation */
iparm[5] = 0; /* Write solution into x */
iparm[6] = 16; /* Default logical fortran unit number for output */
iparm[7] = 2; /* Max numbers of iterative refinement steps */
iparm[8] = 0; /* Not in use*/
iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
iparm[11] = 0; /* Not in use*/
iparm[12] = 0; /* Not in use*/
iparm[13] = 0; /* Output: Number of perturbed pivots */
iparm[14] = 0; /* Not in use*/
iparm[15] = 0; /* Not in use*/
iparm[16] = 0; /* Not in use*/
iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
iparm[18] = -1; /* Output: Mflops for LU factorization */
iparm[19] = 0; /* Output: Numbers of CG Iterations */
maxfct = 1; /* Maximum number of numerical factorizations. */
mnum = 1; /* Which factorization to use. */
msglvl = 0; /* Don't print statistical information in file */
error = 0; /* Initialize error flag */
/* --------------------------------------------------------------------*/
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* --------------------------------------------------------------------*/
for (i = 0; i < 64; i++) {
pt[i] = 0;
}
/* --------------------------------------------------------------------*/
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* --------------------------------------------------------------------*/
phase = 11;
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, a, ia, ja, &idum, &nrhs,
iparm, &msglvl, &ddum, &ddum, &error);
if (error != 0) {
printf("\nERROR during symbolic factorization: %d", error);
exit(1);
}
if(comment=='Y')
{
//printf("\nReordering completed ... ");
//printf("\nNumber of nonzeros in factors = %d", iparm[17]);
//printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
}
/* --------------------------------------------------------------------*/
/* .. Numerical factorization.*/
/* --------------------------------------------------------------------*/
phase = 22;
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, a, ia, ja, &idum, &nrhs,
iparm, &msglvl, &ddum, &ddum, &error);
if (error != 0) {
printf("\nERROR during numerical factorization: %d", error);
exit(2);
}
if(comment=='Y')
//printf("\nFactorization completed ... ");
/* --------------------------------------------------------------------*/
/* .. Back substitution and iterative refinement. */
/* --------------------------------------------------------------------*/
phase = 33;
iparm[7] = 2; /* Max numbers of iterative refinement steps. */
/* Set right hand side to one.*/

/*
for (i = 0; i < n; i++) {
b[i] = 1;
}
*/ 
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, a, ia, ja, &idum, &nrhs,
iparm, &msglvl, _rhs, _lhs, &error);
if (error != 0) {
printf("\nERROR during solution: %d", error);
exit(3);
}
if(comment=='Y')
//printf("\nSolve completed ... ");


/*
printf("\nThe solution of the system is: ");
for (i = 0; i < n; i++) {
printf("\n x [%d] = % f", i, _lhs[i] );
}
printf ("\n");
*/ 
/* --------------------------------------------------------------------*/
/* .. Termination and release of memory. */
/* --------------------------------------------------------------------*/
phase = -1; /* Release internal memory. */
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, &ddum, ia, ja, &idum, &nrhs,
iparm, &msglvl, &ddum, &ddum, &error);
delete []a;
delete []ia;
delete []ja;
return 0;
cout<<endl;
}



void node_allocater(Node _nodes[], ifstream &_nodes_infile, int &_numberofnodes)
{

//Start allocating 
for(int n=0;n<_numberofnodes;n++)
{
_nodes[n].X_lamb=new double[4];
_nodes[n].u_lamb=new double[4];
_nodes[n].du_lamb=new double[4];
_nodes[n].natural_bc=new int[4];
_nodes[n].natural_val=new double[4];
_nodes[n].nonzero_natural_bc=new int[4];
_nodes[n].neumann_bc=new int[4];
_nodes[n].neumann_val=new double[4];
}
//End Allocating 


char current_line[256];
_nodes_infile.getline(current_line,256);
int condition=0;
while(condition==0)
{
_nodes_infile.getline(current_line,256);
if(current_line[1]=='N' && current_line[2]=='O' && current_line[3]=='D' && current_line[4]=='E')
condition=1; 
}

for(int n=0;n<_numberofnodes;n++)
{
char comma;
_nodes_infile>>_nodes[n].nodelabel; 
_nodes_infile>>comma;
for(int i=0;i<2;i++)
{
_nodes_infile>>_nodes[n].X_lamb[i]; 
_nodes_infile>>comma;
}
_nodes_infile>>_nodes[n].X_lamb[2];
_nodes[n].X_lamb[3]=0.0; 
for(int i=0;i<4;i++)
{
_nodes[n].u_lamb[i]=0.0;
_nodes[n].du_lamb[i]=0.0;
_nodes[n].natural_bc[i]=0;
_nodes[n].natural_val[i]=0.0;
_nodes[n].nonzero_natural_bc[i]=0;
_nodes[n].neumann_bc[i]=0;
_nodes[n].neumann_val[i]=0.0;
}
_nodes[n].nS=0.8; //solid volume fraction in the reference coordinate system 
_nodes[n].ns=0.0; // '' in the current coordinate system 
_nodes[n].inpermeability=0.08;
_nodes[n].layup=0.0;
_nodes[n].neighbour_nodes.push_back(n);
}

/*
for(int n=0;n<_numberofnodes;n++)
{
cout<<_nodes[n].nodelabel<<'\t';
cout<<_nodes[n].X[0]<<'\t';
cout<<_nodes[n].X[1]<<'\t';
cout<<_nodes[n].X[2]<<endl;
}
*/
}

void element_allocater(Node _nodes[], Element _elements[], ifstream &_elements_infile, int &_numberofnodes, int &_numberofelements)
{

//Start allocating 
for(int e=0;e<_numberofelements;e++)
{
_elements[e].conname=new int[8];
_elements[e].connid=new int[8];
_elements[e].elementnodes=new Node[8];
}
//End Allocating 


char current_line[256];
_elements_infile.getline(current_line,256);
int condition=0;
while(condition==0)
{
_elements_infile.getline(current_line,256);
//cout<<current_line<<endl;
if(current_line[1]=='E' && current_line[2]=='L' && current_line[3]=='E' && current_line[4]=='M' 
&& current_line[5]=='E' && current_line[6]=='N' && current_line[7]=='T')
condition=1; 
}

for(int e=0;e<_numberofelements;e++)
{
char comma;
_elements_infile>>_elements[e].elementlabel; 
_elements_infile>>comma;
for(int i=0;i<7;i++)
{
_elements_infile>>_elements[e].conname[i]; 
_elements_infile>>comma;
}
_elements_infile>>_elements[e].conname[7]; 

//assignment of the connectivity id 
for(int n=0;n<_numberofnodes;n++)
for(int i=0;i<8;i++)
if(_nodes[n].nodelabel==_elements[e].conname[i])
_elements[e].connid[i]=n;

for(int i=0;i<6;i++)
_elements[e].preS[i]=0.0;

_elements[e].preS[0]=-10.0;
_elements[e].preS[1]=-10.0;
_elements[e].preS[2]=-10.0;

for(int i=0;i<8;i++)
for(int j=0;j<8;j++)
_nodes[_elements[e].connid[i]].Push_Neighbour_In(_elements[e].connid[j]);



}

/*
for(int e=0;e<_numberofelements;e++)
{
cout<<_elements[e].elementlabel<<'\t';
for(int i=0;i<7;i++)
cout<<_elements[e].conname[i]<<'\t';
cout<<_elements[e].conname[8]<<endl;
}
*/
}


void Node::Push_Neighbour_In(int _n)
{
int already_inside=0;
for(int i=0;i<neighbour_nodes.size();i++)
if(_n==neighbour_nodes[i])
already_inside=1;

if(already_inside==0)
neighbour_nodes.push_back(_n);


};



void K_allocater(Node _nodes[], csrM & _K, int** _Kloc,  int** _eft, int _numberofnodes, int _numberofelements)
{
//Allocation through nodes into csr 
/*
for(int n=0;n<_numberofnodes;n++)
for(int s=0;s<_nodes[n].neighbour_nodes.size();s++)
for(int i=0;i<4;i++)
for(int j=0;j<4;j++)
{
_K.Push_Value_incsrM(0.0,n*4+i,_nodes[n].neighbour_nodes[s]*4+j);
double percent_finished=(double(n)+1)*100.0/double(_numberofnodes);
cout<<"%"<<percent_finished<<" compressed sparse row stiffness is allocated              "<<'\r';
}
*/

//Allocation through elements into csr 
/*
for(int e=0;e<_numberofelements;e++)
for(int m=0;m<32;m++)
for(int n=0;n<32;n++)
{
_K.Push_Value_incsrM(0.0,_eft[m][e],_eft[n][e]);
double percent_finished=(double(e)+1)*100.0/double(_numberofelements);
cout<<"%"<<percent_finished<<" compressed sparse row stiffness is allocated"<<'\r'; 
}
*/



for(int e=0;e<_numberofelements;e++)

for(int m=0;m<32;m++)
for(int n=0;n<32;n++)
_Kloc[_eft[m][e]][_eft[n][e]]=1; 

int counter=0;
for(int i=0;i<4*_numberofnodes;i++)
{
int newrow=1;
for(int j=0;j<4*_numberofnodes;j++)
if(_Kloc[i][j]==1)
{
_K.vM.push_back(0.0);
_K.jM.push_back(j);
if(newrow==1)
_K.iM.push_back(counter);
counter++;
newrow=0;
}
double percent_finished=(double(i)+1)*100.0/double(_numberofnodes*4);
cout<<"%"<<percent_finished<<" compressed sparse row stiffness is allocated"<<'\r'; 
}


_K.iM.push_back(_K.vM.size());
cout<<_K.vM.size()<<endl;

};






void K_element(Node _nodes[],
	       Element _elements[],
               int _currentelement,
               ofstream &_debugfile,
	       double t)
{


double beta_2=0.25;
double beta_1=0.5;

double Kgeo[32][32];
for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_elements[_currentelement].Kelem[i][j]=0.0; 

int gauss_counter=0;
//call for the standard gauss points 
for(double g3=-1.0;g3<2.0;g3=g3+2.0)
for(double g2=-1.0;g2<2.0;g2=g2+2.0)
for(double g1=-1.0;g1<2.0;g1=g1+2.0)
{



double current_layup; 
current_layup=_elements[_currentelement].layup[gauss_counter];
gauss_counter++;

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
Kgeo[i][j]=0.0;

//local declerations 
double detj=0.0;
double B[6][32];//strain displacement matrix at the gauss point
double B_rate[6][32];
double B_zero[6][32];//strain displacement matrix at the gauss point
double F[9];//deformation gradient at the gauss point 
double D[6][6];//elasticity matrix at the gauss point
double Bt_dot_D[32][6];
double Bt[32][6];
double S[6];
double N[8]; //Shape Functions at the current stage, necessary for the interpolation of the lambda at the gauss location 
double current_lambda=0.0;
double ksi=g1*1.0/sqrt(3);
double eta=g2*1.0/sqrt(3);
double zeta=g3*1.0/sqrt(3);


//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N[shape_counter]=0.125*(1+sign_ksi*ksi)*(1+sign_eta*eta)*(1+sign_zeta*zeta);
shape_counter++;
}
double temp=N[2];
N[2]=N[3];
N[3]=temp;
temp=N[6];
N[6]=N[7];
N[7]=temp;
//end definition shape value vector

//interpolation of the current lambda. It is assummed that the reference pore pressure is to be taken here 
//this is a quite strange problem 
for(int s=0;s<8;s++)
current_lambda=current_lambda+N[s]*_nodes[_elements[_currentelement].connid[s]].u_lamb[3];

for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
{
B[j][i]=0.0;
B_zero[j][i]=0.0;
}



//computation of the noonlinear strain-displacement matrix and the deformation gradient at the gauss points
Grad_at_Gauss(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),detj,_currentelement,B,B_zero,F,_debugfile);

/*
_debugfile<<"B_zero"<<endl;
for(int i=0;i<6;i++)
for(int j=0;j<32;j++)
_debugfile<<i<<'\t'<<j<<'\t'<<B_zero[i][j]<<endl;
*/


//the material coefficients are defined inside of the function beneath
Elasticity_at_Gauss(D,F,S,_debugfile,current_lambda); //for the time being this function is independent of the space
//Contribution of the geometric stiffness 




K_gauss_geometrical(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,F,S,Kgeo,detj,_currentelement,_debugfile);//TO BE FINISHED 




//start tpm stiffness parts. 
double K_WMat[32][32];
double K_Wlambda[32][32];
double K_Wlambdau[32][32];
double K_M1uJ[32][32];
double K_M1uGeo[32][32];
double K_M1uMat[32][32];
double K_M1uN[32][32];
double K_M2uJ[32][32];
double K_M2uMat[32][32];
double K_M2lambda[32][32];


for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
{
K_WMat[i][j]=0.0;
K_Wlambda[i][j]=0.0;
K_Wlambdau[i][j]=0.0;
K_M1uJ[i][j]=0.0;
K_M1uGeo[i][j]=0.0;
K_M1uMat[i][j]=0.0;
K_M1uN[i][j]=0.0;
K_M2uJ[i][j]=0.0;
K_M2uMat[i][j]=0.0;
K_M2lambda[i][j]=0.0;		
}

double coeff;
if(t>0.0)
coeff=1.0*beta_2/beta_1/t;
else if(t==0.0)
coeff=1.0;

//for(int i=0;i<32;i++)
//for(int j=0;j<6;j++)
//B_rate[i][j]=B[j][i]-B_zero[j][i];


for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
B_rate[j][i]=B[j][i];


//this was open previously 
//if(t==0.0)
for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
Bt[i][j]=B_rate[j][i];





//remember the B matrix in the internal forces are the rate one, (the variation of the rate of the Cauchy-Green Strain tensor )
calc_K_Wlambdau_at_gauss(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B_rate,D,F,S,K_Wlambdau,detj,_currentelement,_debugfile,t,coeff,current_lambda);
calc_K_Wlambda_at_gauss(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B_rate,D,F,S,K_Wlambda,detj,_currentelement,_debugfile,t,coeff);

calc_K_M1uJ_at_gauss(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,F,S,K_M1uJ,detj,_currentelement,_debugfile,t,coeff);
calc_K_M1uGeo_at_gauss(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B_rate,D,F,S,K_M1uGeo,detj,_currentelement,_debugfile,t,coeff);
//In the case below the matrix is given as full
calc_K_M1uMat_at_gauss(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,F,S,K_M1uMat,detj,_currentelement,_debugfile,t,coeff);
//calc_K_M1uN_at_gauss(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,F,S,K_M1uN,detj,_currentelement,_debugfile,t,coeff);

calc_K_M2uJ_at_gauss(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,F,S,K_M2uJ,detj,_currentelement,_debugfile,t,coeff,current_layup);
calc_K_M2uMat_at_gauss(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,F,S,K_M2uMat,detj,_currentelement,_debugfile,t,current_layup);
calc_K_M2lambda_at_gauss(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,F,S,K_M2lambda,detj,_currentelement,_debugfile,t,current_layup);

//end tpm stiffness parts 



for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
Bt_dot_D[i][j]=0.0; 

for(int i=0;i<32;i++)
for(int k=0;k<6;k++)
for(int j=0;j<6;j++)
Bt_dot_D[i][k]=Bt_dot_D[i][k]+Bt[i][j]*D[j][k]; 

for(int i=0;i<32;i++)
for(int k=0;k<32;k++)
for(int j=0;j<6;j++)
K_WMat[i][k]=K_WMat[i][k]+Bt_dot_D[i][j]*B[j][k];



for(int i=0;i<32;i++)
for(int k=0;k<32;k++)
_elements[_currentelement].Kelem[i][k]=_elements[_currentelement].Kelem[i][k]+K_WMat[i][k]*coeff;

/*
if(_currentelement==5)
{
_debugfile<<"K_Mat"<<endl;
_debugfile<<K_WMat[0][0]<<'\t'<<K_WMat[0][1]<<'\t'<<K_WMat[0][2]<<'\t'<<K_WMat[0][3]<<endl;
_debugfile<<K_WMat[1][0]<<'\t'<<K_WMat[1][1]<<'\t'<<K_WMat[1][2]<<'\t'<<K_WMat[1][3]<<endl;
_debugfile<<K_WMat[2][0]<<'\t'<<K_WMat[2][1]<<'\t'<<K_WMat[2][2]<<'\t'<<K_WMat[2][3]<<endl;
_debugfile<<K_WMat[3][0]<<'\t'<<K_WMat[3][1]<<'\t'<<K_WMat[3][2]<<'\t'<<K_WMat[3][3]<<endl;
}
*/

for(int i=0;i<32;i++)
for(int k=0;k<32;k++)
_elements[_currentelement].Kelem[i][k]=_elements[_currentelement].Kelem[i][k]+(Kgeo[i][k])*coeff;
//_elements[_currentelement].Kelem[i][k]=_elements[_currentelement].Kelem[i][k]+(Kgeo[i][k]+K_Wlambdau[i][k])*coeff;
//_elements[_currentelement].Kelem[i][k]=_elements[_currentelement].Kelem[i][k]+(Kgeo[i][k]+K_Wlambdau[i][k])*1.0*beta_2/beta_1/t;



for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
//if(_nodes[_elements[_currentelement].connid[i/4]].neumann_bc[i%4]==0 && _nodes[_elements[_currentelement].connid[j/4]].neumann_bc[j%4]==0)
{
_elements[_currentelement].Kelem[i][j]=_elements[_currentelement].Kelem[i][j]
+K_M1uJ[i][j]
+K_M1uGeo[i][j]
+K_M1uMat[i][j]
+K_M1uN[i][j]
+K_Wlambda[i][j]
+K_Wlambdau[i][j]
+K_M2lambda[i][j];
//if(_nodes[_elements[_currentelement].connid[i/4]].neumann_bc[i%4]==0 && _nodes[_elements[_currentelement].connid[j/4]].neumann_bc[j%4]==0)
_elements[_currentelement].Kelem[i][j]=_elements[_currentelement].Kelem[i][j]+
+K_M2uJ[i][j]
+K_M2uMat[i][j];
}


/*
if(_currentelement==5)
{

_debugfile<<"KGeo"<<endl;
_debugfile<<Kgeo[0][0]<<'\t'<<Kgeo[0][1]<<'\t'<<Kgeo[0][2]<<'\t'<<Kgeo[0][3]<<endl;
_debugfile<<Kgeo[1][0]<<'\t'<<Kgeo[1][1]<<'\t'<<Kgeo[1][2]<<'\t'<<Kgeo[1][3]<<endl;
_debugfile<<Kgeo[2][0]<<'\t'<<Kgeo[2][1]<<'\t'<<Kgeo[2][2]<<'\t'<<Kgeo[2][3]<<endl;
_debugfile<<Kgeo[3][0]<<'\t'<<Kgeo[3][1]<<'\t'<<Kgeo[3][2]<<'\t'<<Kgeo[3][3]<<endl;
//_debugfile<<"S"<<endl;
//_debugfile<<F[0]<<'\t'<<F[1]<<'\t'<<F[2]<<'\t'<<F[3]<<'\t'<<F[4]<<'\t'<<F[5]<<endl;
_debugfile<<"K_Wlambda"<<endl;
_debugfile<<K_Wlambda[0][0]<<'\t'<<K_Wlambda[0][1]<<'\t'<<K_Wlambda[0][2]<<'\t'<<K_Wlambda[0][3]<<endl;
_debugfile<<K_Wlambda[1][0]<<'\t'<<K_Wlambda[1][1]<<'\t'<<K_Wlambda[1][2]<<'\t'<<K_Wlambda[1][3]<<endl;
_debugfile<<K_Wlambda[2][0]<<'\t'<<K_Wlambda[2][1]<<'\t'<<K_Wlambda[2][2]<<'\t'<<K_Wlambda[2][3]<<endl;
_debugfile<<K_Wlambda[3][0]<<'\t'<<K_Wlambda[3][1]<<'\t'<<K_Wlambda[3][2]<<'\t'<<K_Wlambda[3][3]<<endl;

_debugfile<<"K_Wlambdau"<<endl;
_debugfile<<K_Wlambdau[0][0]<<'\t'<<K_Wlambdau[0][1]<<'\t'<<K_Wlambdau[0][2]<<'\t'<<K_Wlambdau[0][3]<<endl;
_debugfile<<K_Wlambdau[1][0]<<'\t'<<K_Wlambdau[1][1]<<'\t'<<K_Wlambdau[1][2]<<'\t'<<K_Wlambdau[1][3]<<endl;
_debugfile<<K_Wlambdau[2][0]<<'\t'<<K_Wlambdau[2][1]<<'\t'<<K_Wlambdau[2][2]<<'\t'<<K_Wlambdau[2][3]<<endl;
_debugfile<<K_Wlambdau[3][0]<<'\t'<<K_Wlambdau[3][1]<<'\t'<<K_Wlambdau[3][2]<<'\t'<<K_Wlambdau[3][3]<<endl;

_debugfile<<"K_M1uj"<<endl;
_debugfile<<K_M1uJ[0][0]<<'\t'<<K_M1uJ[0][1]<<'\t'<<K_M1uJ[0][2]<<'\t'<<K_M1uJ[0][3]<<endl;
_debugfile<<K_M1uJ[1][0]<<'\t'<<K_M1uJ[1][1]<<'\t'<<K_M1uJ[1][2]<<'\t'<<K_M1uJ[1][3]<<endl;
_debugfile<<K_M1uJ[2][0]<<'\t'<<K_M1uJ[2][1]<<'\t'<<K_M1uJ[2][2]<<'\t'<<K_M1uJ[2][3]<<endl;
_debugfile<<K_M1uJ[3][0]<<'\t'<<K_M1uJ[3][1]<<'\t'<<K_M1uJ[3][2]<<'\t'<<K_M1uJ[3][3]<<endl;
_debugfile<<"K_M1uGeo"<<endl;
_debugfile<<K_M1uGeo[0][0]<<'\t'<<K_M1uGeo[0][1]<<'\t'<<K_M1uGeo[0][2]<<'\t'<<K_M1uGeo[0][3]<<endl;
_debugfile<<K_M1uGeo[1][0]<<'\t'<<K_M1uGeo[1][1]<<'\t'<<K_M1uGeo[1][2]<<'\t'<<K_M1uGeo[1][3]<<endl;
_debugfile<<K_M1uGeo[2][0]<<'\t'<<K_M1uGeo[2][1]<<'\t'<<K_M1uGeo[2][2]<<'\t'<<K_M1uGeo[2][3]<<endl;
_debugfile<<K_M1uGeo[3][0]<<'\t'<<K_M1uGeo[3][1]<<'\t'<<K_M1uGeo[3][2]<<'\t'<<K_M1uGeo[3][3]<<endl;
_debugfile<<"K_M1uMat"<<endl;
_debugfile<<K_M1uMat[0][0]<<'\t'<<K_M1uMat[0][1]<<'\t'<<K_M1uMat[0][2]<<'\t'<<K_M1uMat[0][3]<<endl;
_debugfile<<K_M1uMat[1][0]<<'\t'<<K_M1uMat[1][1]<<'\t'<<K_M1uMat[1][2]<<'\t'<<K_M1uMat[1][3]<<endl;
_debugfile<<K_M1uMat[2][0]<<'\t'<<K_M1uMat[2][1]<<'\t'<<K_M1uMat[2][2]<<'\t'<<K_M1uMat[2][3]<<endl;
_debugfile<<K_M1uMat[3][0]<<'\t'<<K_M1uMat[3][1]<<'\t'<<K_M1uMat[3][2]<<'\t'<<K_M1uMat[3][3]<<endl;
_debugfile<<"K_M1uN"<<endl;
_debugfile<<K_M1uN[0][0]<<'\t'<<K_M1uN[0][1]<<'\t'<<K_M1uN[0][2]<<'\t'<<K_M1uN[0][3]<<endl;
_debugfile<<K_M1uN[1][0]<<'\t'<<K_M1uN[1][1]<<'\t'<<K_M1uN[1][2]<<'\t'<<K_M1uN[1][3]<<endl;
_debugfile<<K_M1uN[2][0]<<'\t'<<K_M1uN[2][1]<<'\t'<<K_M1uN[2][2]<<'\t'<<K_M1uN[2][3]<<endl;
_debugfile<<K_M1uN[3][0]<<'\t'<<K_M1uN[3][1]<<'\t'<<K_M1uN[3][2]<<'\t'<<K_M1uN[3][3]<<endl;


_debugfile<<"K_M2uJ"<<endl;
_debugfile<<K_M2uJ[0][0]<<'\t'<<K_M2uJ[0][1]<<'\t'<<K_M2uJ[0][2]<<'\t'<<K_M2uJ[0][3]<<endl;
_debugfile<<K_M2uJ[1][0]<<'\t'<<K_M2uJ[1][1]<<'\t'<<K_M2uJ[1][2]<<'\t'<<K_M2uJ[1][3]<<endl;
_debugfile<<K_M2uJ[2][0]<<'\t'<<K_M2uJ[2][1]<<'\t'<<K_M2uJ[2][2]<<'\t'<<K_M2uJ[2][3]<<endl;
_debugfile<<K_M2uJ[3][0]<<'\t'<<K_M2uJ[3][1]<<'\t'<<K_M2uJ[3][2]<<'\t'<<K_M2uJ[3][3]<<endl;
_debugfile<<"K_M2uMat"<<endl;
_debugfile<<K_M2uMat[0][0]<<'\t'<<K_M2uMat[0][1]<<'\t'<<K_M2uMat[0][2]<<'\t'<<K_M2uMat[0][3]<<endl;
_debugfile<<K_M2uMat[1][0]<<'\t'<<K_M2uMat[1][1]<<'\t'<<K_M2uMat[1][2]<<'\t'<<K_M2uMat[1][3]<<endl;
_debugfile<<K_M2uMat[2][0]<<'\t'<<K_M2uMat[2][1]<<'\t'<<K_M2uMat[2][2]<<'\t'<<K_M2uMat[2][3]<<endl;
_debugfile<<K_M2uMat[3][0]<<'\t'<<K_M2uMat[3][1]<<'\t'<<K_M2uMat[3][2]<<'\t'<<K_M2uMat[3][3]<<endl;
_debugfile<<"K_M2lambda"<<endl;
_debugfile<<K_M2lambda[0][0]<<'\t'<<K_M2lambda[0][1]<<'\t'<<K_M2lambda[0][2]<<'\t'<<K_M2lambda[0][3]<<endl;
_debugfile<<K_M2lambda[1][0]<<'\t'<<K_M2lambda[1][1]<<'\t'<<K_M2lambda[1][2]<<'\t'<<K_M2lambda[1][3]<<endl;
_debugfile<<K_M2lambda[2][0]<<'\t'<<K_M2lambda[2][1]<<'\t'<<K_M2lambda[2][2]<<'\t'<<K_M2lambda[2][3]<<endl;
_debugfile<<K_M2lambda[3][0]<<'\t'<<K_M2lambda[3][1]<<'\t'<<K_M2lambda[3][2]<<'\t'<<K_M2lambda[3][3]<<endl;

}

*/

}





};




//Strain-Displacement Matrix at a specific gauss point 
void Grad_at_Gauss(Node _nodes[],
	           Element _elements[],
	           double _ksi, double _eta, double _zeta,
                   double &_detj,
                   int _currentelement,
	           double B[][32], //Strain displacement matrix 
		   double B_zero[][32],
		   double F[], //deformation gradient 
                   ofstream &_debugfile
		   )
{
double jacob[3][3];
double detj; 
double shapederivative[9][24];
double modifiedinversejacob_B[6][9];
double modifiedinversejacob_F[9][9];
double Finterpolator[9][24];
double Gradu[9];
double GraduT_tensor[3][3];
double BLin[6][24];
double BNLin[6][24];
double B_zero_mech[6][24];

for(int i=0;i<6;i++)
for(int j=0;j<24;j++)
{
//With BLin it is not necessarily meant the Linear Strain displacement matrix 
//It represents only the previous multiplication of the Nonlinear strain displacement matrix 
BLin[i][j]=0.0; 
BNLin[i][j]=0.0; 
B_zero_mech[i][j]=0.0;
}


for(int i=0;i<6;i++)
for(int j=0;j<32;j++)
{
B[i][j]=0.0;
B_zero[i][j]=0.0;
}
for(int i=0;i<9;i++)
for(int j=0;j<24;j++)
{
shapederivative[i][j]=0.0;
Finterpolator[i][j]=0.0;
Gradu[i]=0.0;
}
for(int i=0;i<6;i++)
for(int j=0;j<9;j++)
modifiedinversejacob_B[i][j]=0.0;


for(int i=0;i<9;i++)
for(int j=0;j<9;j++)
modifiedinversejacob_F[i][j]=0.0;


shapederivative[0][0]=shapederivative[3][1]=shapederivative[6][2]=-0.125*(1.0-_eta)*(1.0-_zeta);
shapederivative[1][0]=shapederivative[4][1]=shapederivative[7][2]=-0.125*(1.0-_ksi)*(1.0-_zeta);
shapederivative[2][0]=shapederivative[5][1]=shapederivative[8][2]=-0.125*(1.0-_ksi)*(1.0-_eta);
//2 three columns 
shapederivative[0][3]=shapederivative[3][4]=shapederivative[6][5]=0.125*(1.0-_eta)*(1.0-_zeta);
shapederivative[1][3]=shapederivative[4][4]=shapederivative[7][5]=-0.125*(1.0+_ksi)*(1.0-_zeta);
shapederivative[2][3]=shapederivative[5][4]=shapederivative[8][5]=-0.125*(1.0+_ksi)*(1.0-_eta);
//3 three columns 
shapederivative[0][6]=shapederivative[3][7]=shapederivative[6][8]=0.125*(1.0+_eta)*(1.0-_zeta);
shapederivative[1][6]=shapederivative[4][7]=shapederivative[7][8]=0.125*(1.0+_ksi)*(1.0-_zeta);
shapederivative[2][6]=shapederivative[5][7]=shapederivative[8][8]=-0.125*(1.0+_ksi)*(1.0+_eta);
//4 three columns 
shapederivative[0][9]=shapederivative[3][10]=shapederivative[6][11]=-0.125*(1.0+_eta)*(1.0-_zeta);
shapederivative[1][9]=shapederivative[4][10]=shapederivative[7][11]=0.125*(1.0-_ksi)*(1.0-_zeta);
shapederivative[2][9]=shapederivative[5][10]=shapederivative[8][11]=-0.125*(1.0-_ksi)*(1.0+_eta);
//5 three columns 
shapederivative[0][12]=shapederivative[3][13]=shapederivative[6][14]=-0.125*(1.0-_eta)*(1.0+_zeta);
shapederivative[1][12]=shapederivative[4][13]=shapederivative[7][14]=-0.125*(1.0-_ksi)*(1.0+_zeta);
shapederivative[2][12]=shapederivative[5][13]=shapederivative[8][14]=0.125*(1.0-_ksi)*(1.0-_eta);
//6 three columns 
shapederivative[0][15]=shapederivative[3][16]=shapederivative[6][17]=0.125*(1.0-_eta)*(1.0+_zeta);
shapederivative[1][15]=shapederivative[4][16]=shapederivative[7][17]=-0.125*(1.0+_ksi)*(1.0+_zeta);
shapederivative[2][15]=shapederivative[5][16]=shapederivative[8][17]=0.125*(1.0+_ksi)*(1.0-_eta);
//7 three columns 
shapederivative[0][18]=shapederivative[3][19]=shapederivative[6][20]=0.125*(1.0+_eta)*(1.0+_zeta);
shapederivative[1][18]=shapederivative[4][19]=shapederivative[7][20]=0.125*(1.0+_ksi)*(1.0+_zeta);
shapederivative[2][18]=shapederivative[5][19]=shapederivative[8][20]=0.125*(1.0+_ksi)*(1.0+_eta);
//8 three columns 
shapederivative[0][21]=shapederivative[3][22]=shapederivative[6][23]=-0.125*(1.0+_eta)*(1.0+_zeta);
shapederivative[1][21]=shapederivative[4][22]=shapederivative[7][23]=0.125*(1.0-_ksi)*(1.0+_zeta);
shapederivative[2][21]=shapederivative[5][22]=shapederivative[8][23]=0.125*(1.0-_ksi)*(1.0+_eta);

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0; 

//CAUTION! MOST PROBABLY THE TOTAL DEFORMATION AND THUS THE CURRENT COORDINATES HAVE TO BE USED HERE ! 
for(int i=0;i<3;i++)
for(int j=0;j<8;j++)
{
jacob[0][i]=jacob[0][i]+shapederivative[i][j*3]*_nodes[_elements[_currentelement].connid[j]].X_lamb[0];
jacob[1][i]=jacob[1][i]+shapederivative[i][j*3]*_nodes[_elements[_currentelement].connid[j]].X_lamb[1];
jacob[2][i]=jacob[2][i]+shapederivative[i][j*3]*_nodes[_elements[_currentelement].connid[j]].X_lamb[2];
}

detj=det(jacob);
_detj=detj;
modifiedinversejacob_B[0][0]=modifiedinversejacob_B[3][3]=modifiedinversejacob_B[5][6]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
modifiedinversejacob_B[0][1]=modifiedinversejacob_B[3][4]=modifiedinversejacob_B[5][7]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
modifiedinversejacob_B[0][2]=modifiedinversejacob_B[3][5]=modifiedinversejacob_B[5][8]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;

modifiedinversejacob_B[1][3]=modifiedinversejacob_B[3][0]=modifiedinversejacob_B[4][6]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
modifiedinversejacob_B[1][4]=modifiedinversejacob_B[3][1]=modifiedinversejacob_B[4][7]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
modifiedinversejacob_B[1][5]=modifiedinversejacob_B[3][2]=modifiedinversejacob_B[4][8]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;

modifiedinversejacob_B[2][6]=modifiedinversejacob_B[4][3]=modifiedinversejacob_B[5][0]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj;
modifiedinversejacob_B[2][7]=modifiedinversejacob_B[4][4]=modifiedinversejacob_B[5][1]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;
modifiedinversejacob_B[2][8]=modifiedinversejacob_B[4][5]=modifiedinversejacob_B[5][2]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;







modifiedinversejacob_F[0][0]=modifiedinversejacob_F[3][3]=modifiedinversejacob_F[6][6]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
modifiedinversejacob_F[0][1]=modifiedinversejacob_F[3][4]=modifiedinversejacob_F[6][7]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
modifiedinversejacob_F[0][2]=modifiedinversejacob_F[3][5]=modifiedinversejacob_F[6][8]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;

modifiedinversejacob_F[1][0]=modifiedinversejacob_F[4][3]=modifiedinversejacob_F[7][6]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
modifiedinversejacob_F[1][1]=modifiedinversejacob_F[4][4]=modifiedinversejacob_F[7][7]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
modifiedinversejacob_F[1][2]=modifiedinversejacob_F[4][5]=modifiedinversejacob_F[7][8]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;


modifiedinversejacob_F[2][0]=modifiedinversejacob_F[5][3]=modifiedinversejacob_F[8][6]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj;
modifiedinversejacob_F[2][1]=modifiedinversejacob_F[5][4]=modifiedinversejacob_F[8][7]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;
modifiedinversejacob_F[2][2]=modifiedinversejacob_F[5][5]=modifiedinversejacob_F[8][8]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;

/*
if(_currentelement==0)
_debugfile<<"det J: "<<detj<<endl;


if(_currentelement==0)
{
_debugfile<<"Something "<<_currentelement<<endl;
for(int i=0;i<3;i++)
for(int k=0;k<3;k++)
_debugfile<<jacob[i][k]<<setw(10)<<endl;
}
*/

for(int i=0;i<9;i++)
for(int k=0;k<24;k++)
for(int j=0;j<9;j++)
Finterpolator[i][k]=Finterpolator[i][k]+modifiedinversejacob_F[i][j]*shapederivative[j][k];





for(int i=0;i<9;i++)
for(int k=0;k<24;k++)
Gradu[i]=Gradu[i]+Finterpolator[i][k]*_nodes[_elements[_currentelement].connid[k/3]].u_lamb[k%3];



for(int i=0;i<9;i++)
{
F[i]=Gradu[i];
if(i==0 || i==4 || i==8)
F[i]=F[i]+1.0;  
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
GraduT_tensor[j][i]=F[i*3+j];


for(int i=0;i<6;i++)
for(int k=0;k<24;k++)
for(int j=0;j<9;j++)
BLin[i][k]=BLin[i][k]+modifiedinversejacob_B[i][j]*shapederivative[j][k];






for(int i=0;i<6;i++)
for(int k=0;k<24;k++)
for(int j=0;j<24;j++)
BNLin[i][k]=BNLin[i][k]+BLin[i][j]*GraduT_tensor[j%3][k%3]*cronecker(k/3,j/3);



//instead of deformation gradient, the zero deformation gradient is taken here 

for(int i=0;i<6;i++)
for(int k=0;k<24;k++)
for(int j=0;j<24;j++)
B_zero_mech[i][k]=B_zero_mech[i][k]+BLin[i][j]*cronecker(j%3,k%3)*cronecker(k/3,j/3);



for(int i=0;i<6;i++)
for(int k=0;k<24;k++)
B[i][k+k/3]=BNLin[i][k];

for(int i=0;i<6;i++)
for(int k=0;k<24;k++)
B_zero[i][k+k/3]=B_zero_mech[i][k];



/*
if(_currentelement==5)
{
_debugfile<<"B "<<_currentelement<<endl;
_debugfile<<B[0][0]<<'\t'<<B[0][1]<<'\t'<<B[0][2]<<'\t'<<B[0][3]<<endl;
_debugfile<<B[1][0]<<'\t'<<B[1][1]<<'\t'<<B[1][2]<<'\t'<<B[1][3]<<endl;
_debugfile<<B[2][0]<<'\t'<<B[2][1]<<'\t'<<B[2][2]<<'\t'<<B[2][3]<<endl;
_debugfile<<B[3][0]<<'\t'<<B[3][1]<<'\t'<<B[3][2]<<'\t'<<B[3][3]<<endl;
_debugfile<<B[4][0]<<'\t'<<B[4][1]<<'\t'<<B[4][2]<<'\t'<<B[4][3]<<endl;
_debugfile<<B[5][0]<<'\t'<<B[5][1]<<'\t'<<B[5][2]<<'\t'<<B[5][3]<<endl;
}
*/

/*
if(_currentelement==0)
{
_debugfile<<"Something "<<_currentelement<<endl;
for(int i=0;i<24;i++)
{
for(int k=0;k<24;k++)
_debugfile<<GraduT_tensor[i%3][k%3]*cronecker(i/3,k/3)<<setw(10);
_debugfile<<endl;
}
}
*/
};


void Elasticity_at_Gauss(double D[][6],  
		         double F[], //deformation gradient 
			 double S[],
                         ofstream &_debugfile,
			 double lambda)
{
double Ftens[3][3];
double FtensT[3][3];
double FtensInv[3][3];
double Cauchy[3][3]; //Cauchy-Green Strain Tensor 
double CauchyInv[3][3]; //Cauchy-Green Strain Tensor 
double gammablk=100; //gamma coefficient for the bulk energy 
double beta=4.5; //beta coefficient for something else 
double c[3][3][3][3];  //constitutive tensor in tensor notation 
double C[3][3][3][3]; //bulk constitutive tensor 
double E=5000.0;
double v=0.3;
double lam=E*v/(1.0+v)/(1.0-2.0*v);
double Mu=E*(1.0-v)/(1.0+v)/(1.0-2.0*v);

Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];

transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchy); 
inverse(Cauchy,CauchyInv);
double Ithree=det(Cauchy); 
double Ione=trace(Cauchy); 
inverse(Ftens,FtensInv);
double detF=det(Ftens);


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
for(int k=0;k<3;k++)
for(int l=0;l<3;l++)
{
c[i][j][k][l]=(lam*((2.0*detF-1.0)*cronecker(i,j)*cronecker(k,l)-2.0*symmetry_tensor(i,j,k,l)*(detF-1.0))+2.0*Mu/detF*symmetry_tensor(i,j,k,l));
//coefficient 4 is for the cauchy to green-lagrange strain tensor conversion 
}

for(int I=0;I<3;I++)
for(int J=0;J<3;J++)
for(int K=0;K<3;K++)
for(int L=0;L<3;L++)
C[I][J][K][L]=0.0;

for(int I=0;I<3;I++)
for(int J=0;J<3;J++)
for(int K=0;K<3;K++)
for(int L=0;L<3;L++)
for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
for(int k=0;k<3;k++)
for(int l=0;l<3;l++)
{
C[I][J][K][L]=C[I][J][K][L]+FtensInv[I][i]*FtensInv[J][j]*c[i][j][k][l]*FtensInv[K][k]*FtensInv[L][l]*detF;
}


//effect of fluid pore pressure
for(int I=0;I<3;I++)
for(int J=0;J<3;J++)
for(int K=0;K<3;K++)
for(int L=0;L<3;L++)
C[I][J][K][L]=C[I][J][K][L]+1.0*lambda*(CauchyInv[I][K]*CauchyInv[J][L]+CauchyInv[I][L]*CauchyInv[J][K]);




/*
for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
for(int k=0;k<3;k++)
for(int l=0;l<3;l++)
{
cbulk[i][j][k][l]=4*gammablk*beta*pow(Ithree,-1.0*beta)*cronecker(i,j)*cronecker(k,l)
                 +4*gammablk*pow(Ithree,-1.0*beta)*symmetry_tensor(i,j,k,l);
c[i][j][k][l]=+cbulk[i][j][k][l]+c[i][j][k][l]; 
}
*/  
D[0][0]=C[0][0][0][0];
D[0][1]=C[0][0][1][1];
D[0][2]=C[0][0][2][2];
D[0][3]=C[0][0][0][1];
D[0][4]=C[0][0][1][2];
D[0][5]=C[0][0][2][0];
 
D[1][0]=C[1][1][0][0];
D[1][1]=C[1][1][1][1];
D[1][2]=C[1][1][2][2];
D[1][3]=C[1][1][0][1];
D[1][4]=C[1][1][1][2];
D[1][5]=C[1][1][2][0];

D[2][0]=C[2][2][0][0];
D[2][1]=C[2][2][1][1];
D[2][2]=C[2][2][2][2];
D[2][3]=C[2][2][0][1];
D[2][4]=C[2][2][1][2];
D[2][5]=C[2][2][2][0];

D[3][0]=C[0][1][0][0];
D[3][1]=C[0][1][1][1];
D[3][2]=C[0][1][2][2];
D[3][3]=C[0][1][0][1];
D[3][4]=C[0][1][1][2];
D[3][5]=C[0][1][2][0];

D[4][0]=C[1][2][0][0];
D[4][1]=C[1][2][1][1];
D[4][2]=C[1][2][2][2];
D[4][3]=C[1][2][0][1];
D[4][4]=C[1][2][1][2];
D[4][5]=C[1][2][2][0];

D[5][0]=C[2][0][0][0];
D[5][1]=C[2][0][1][1];
D[5][2]=C[2][0][2][2];
D[5][3]=C[2][0][0][1];
D[5][4]=C[2][0][1][2];
D[5][5]=C[2][0][2][0];

/*
_debugfile<<"Something "<<endl;
for(int i=0;i<6;i++)
for(int j=0;j<6;j++)
_debugfile<<D[i][j]<<setw(10);
_debugfile<<endl;
*/

//computation of the stress
double cauchy_stress[3][3];
double finger[3][3];
double Stens[3][3];

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
{
cauchy_stress[i][j]=0.0;
Stens[i][j]=0.0;
}
matrix_product(Ftens,FtensT,finger);

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
cauchy_stress[i][j]=(cauchy_stress[i][j]+lam*(detF-1.0)*cronecker(i,j)+Mu/detF*(finger[i][j]-cronecker(i,j)));
//coefficient 2 is for the cauchy to green lagrange strain tensor conversion 
for(int I=0;I<3;I++)
for(int J=0;J<3;J++)
for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
Stens[I][J]=Stens[I][J]+detF*FtensInv[i][I]*cauchy_stress[i][j]*FtensInv[j][J];

//effect of fluis pore pressure on Second Pioloa Kirchoff stress tensor 
for(int I=0;I<3;I++)
for(int J=0;J<3;J++)
Stens[I][J]=Stens[I][J]-lambda*CauchyInv[I][J];


S[0]=Stens[0][0];
S[1]=Stens[1][1];
S[2]=Stens[2][2];
S[3]=Stens[0][1];
S[4]=Stens[1][2];
S[5]=Stens[2][0];


//thought for the pore pressure 
//for(int i=0;i<3;i++)
//S[i]=S[i]-10000.0;


};







void K_gauss_geometrical(Node _nodes[],
	                 Element _elements[],
	                 double _ksi, double _eta, double _zeta,
			 double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
			 double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
			 double F[],
			 double S[],
			 double _Kgeo[][32],//the geometrical stiffness matrix to be passed 
			 double& detj,
                         int _currentelement,
                         ofstream &_debugfile)
{
double Kgeo_temp[24][24];
double dNdksi[24]; 
double inversejacob[24][24]; 
double dNdX[24]; 
double jacob[3][3]; 
//double E[6]; //Green Lagrange Strain Tensor in Vector format 
//double S[6]; //Second Piola-Krichoff Material in Vector format 
double Ftens[3][3];
double FtensT[3][3];
double Etens[3][3];

/*
for(int i=0;i<6;i++)
{
E[i]=0.0;
S[i]=0.0;
}
*/
dNdksi[0]=-0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[1]=-0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[2]=-0.125*(1.0-_ksi)*(1.0-_eta);


dNdksi[3]=0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[4]=-0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[5]=-0.125*(1.0+_ksi)*(1.0-_eta);


dNdksi[6]=0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[7]=0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[8]=-0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[9]=-0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[10]=0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[11]=-0.125*(1.0-_ksi)*(1.0+_eta);

dNdksi[12]=-0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[13]=-0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[14]=0.125*(1.0-_ksi)*(1.0-_eta);

dNdksi[15]=0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[16]=-0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[17]=0.125*(1.0+_ksi)*(1.0-_eta);

dNdksi[18]=0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[19]=0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[20]=0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[21]=-0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[22]=0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[23]=0.125*(1.0-_ksi)*(1.0+_eta);




for(int i=0;i<24;i++)
{
for(int j=0;j<24;j++)
inversejacob[i][j]=0.0; 
dNdX[i]=0.0; 
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0;


for(int i=0;i<3;i++) // row label 
for(int j=0;j<8;j++) 
{
jacob[i][0]=jacob[i][0]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[0];
jacob[i][1]=jacob[i][1]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[1];
jacob[i][2]=jacob[i][2]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[2];
}

detj=jacob[0][0]*(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])-
     jacob[0][1]*(jacob[1][0]*jacob[2][2]-jacob[2][0]*jacob[1][2])+
     jacob[0][2]*(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1]);


inversejacob[0][0]=inversejacob[3][3]=inversejacob[6][6]=inversejacob[9][9]=inversejacob[12][12]=inversejacob[15][15]=inversejacob[18][18]=inversejacob[21][21]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
inversejacob[0][1]=inversejacob[3][4]=inversejacob[6][7]=inversejacob[9][10]=inversejacob[12][13]=inversejacob[15][16]=inversejacob[18][19]=inversejacob[21][22]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
inversejacob[0][2]=inversejacob[3][5]=inversejacob[6][8]=inversejacob[9][11]=inversejacob[12][14]=inversejacob[15][17]=inversejacob[18][20]=inversejacob[21][23]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj; 

inversejacob[1][0]=inversejacob[4][3]=inversejacob[7][6]=inversejacob[10][9]=inversejacob[13][12]=inversejacob[16][15]=inversejacob[19][18]=inversejacob[22][21]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
inversejacob[1][1]=inversejacob[4][4]=inversejacob[7][7]=inversejacob[10][10]=inversejacob[13][13]=inversejacob[16][16]=inversejacob[19][19]=inversejacob[22][22]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
inversejacob[1][2]=inversejacob[4][5]=inversejacob[7][8]=inversejacob[10][11]=inversejacob[13][14]=inversejacob[16][17]=inversejacob[19][20]=inversejacob[22][23]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;

inversejacob[2][0]=inversejacob[5][3]=inversejacob[8][6]=inversejacob[11][9]=inversejacob[14][12]=inversejacob[17][15]=inversejacob[20][18]=inversejacob[23][21]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;
inversejacob[2][1]=inversejacob[5][4]=inversejacob[8][7]=inversejacob[11][10]=inversejacob[14][13]=inversejacob[17][16]=inversejacob[20][19]=inversejacob[23][22]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;
inversejacob[2][2]=inversejacob[5][5]=inversejacob[8][8]=inversejacob[11][11]=inversejacob[14][14]=inversejacob[17][17]=inversejacob[20][20]=inversejacob[23][23]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;


for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
dNdX[i]=dNdX[i]+inversejacob[i][j]*dNdksi[j];


/*
for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
FtensT[i][j]=F[i*3+j];

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
Ftens[i][j]=FtensT[j][i];

matrix_product(FtensT,Ftens,Etens);


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
Etens[i][j]=Etens[i][j]-cronecker(i,j);

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
Etens[i][j]=0.5*Etens[i][j];

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
E[i+j*3]=Etens[i][j];


for(int i=0;i<6;i++)
for(int j=0;j<6;j++)
S[i]=S[i]+D[i][j]*E[j];
*/


//CAUTION: IT HAS TO BE CHECKED IF THE GEOMETRICAL STIFFNESS MATRIX HAS TO BE 
//INITIALIZED AS ZERO BEFORE THE CALL OF THIS FUNCTION
//,x: i or j 
//,y: i+1 or j+1 
//,z: i+2 or j+2 
for(int i=0;i<24;i++)
for(int k=0;k<24;k++)
Kgeo_temp[i][k]=0.0;

for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
{
if(i%3==0 && j%3==0)
Kgeo_temp[i][j]=Kgeo_temp[i][j]+
(dNdX[i]*S[0]*dNdX[j]
+dNdX[i+1]*S[1]*dNdX[j+1]
+dNdX[i+2]*S[2]*dNdX[j+2]
+dNdX[i]*S[3]*dNdX[j+1]
+dNdX[i+1]*S[3]*dNdX[j]
+dNdX[i+1]*S[4]*dNdX[j+2]
+dNdX[i+2]*S[4]*dNdX[j+1]
+dNdX[i]*S[5]*dNdX[j+2]
+dNdX[i+2]*S[5]*dNdX[j]);//*detj; 
else if(i%3!=0 && j>0)
Kgeo_temp[i][j]=Kgeo_temp[i-1][j-1]; 
}



for(int i=0;i<24;i++)
for(int k=0;k<24;k++)
_Kgeo[i+i/3][k+k/3]=_Kgeo[i+i/3][k+k/3]+Kgeo_temp[i][k];

}


void incremental_force_initial(Node nodes[], Element elements[], double fi[], double numberofloadsteps, double l, ofstream& debugfile,int numberofnodes,int numberofelements, double incremental_force[])
{





double beta_2=0.25;
double beta_1=0.5;



double coeff;
if(l>0.0)
coeff=1.0*beta_2/beta_1/l;
else if(l==0.0)
coeff=1.0;






double internal_M1[32];
double internal_M2[32];




for(int e=0;e<numberofelements;e++)
{
for(int i=0;i<32;i++)
{
internal_M1[i]=0.0;
internal_M2[i]=0.0;
}

int gauss_counter=0;


for(double g3=-1.0;g3<2.0;g3=g3+2.0)
for(double g2=-1.0;g2<2.0;g2=g2+2.0)
for(double g1=-1.0;g1<2.0;g1=g1+2.0)
{



double current_layup; 
current_layup=elements[e].layup[gauss_counter];
gauss_counter++;

//local declerations 
double detj=0.0;
double B[6][32];//strain displacement matrix at the gauss point
double B_zero[6][32];//initial strain displacement matrix at the gauss point
double F[9];//deformation gradient at the gauss point 
double D[6][6];//elasticity matrix at the gauss point
double Bt[32][6];
double S[6];
double N[8]; //Shape Functions at the current stage, necessary for the interpolation of the lambda at the gauss location 
double current_lambda=0.0;
double ksi=g1*1.0/sqrt(3);
double eta=g2*1.0/sqrt(3);
double zeta=g3*1.0/sqrt(3);


//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N[shape_counter]=0.125*(1+sign_ksi*ksi)*(1+sign_eta*eta)*(1+sign_zeta*zeta);
shape_counter++;
}
double temp=N[2];
N[2]=N[3];
N[3]=temp;
temp=N[6];
N[6]=N[7];
N[7]=temp;
//end definition shape value vector

//interpolation of the current lambda. It is assummed that the reference pore pressure is to be taken here 
for(int s=0;s<8;s++)
current_lambda=current_lambda+N[s]*nodes[elements[e].connid[s]].u_lamb[3];



//computation of the noonlinear strain-displacement matrix and the deformation gradient at the gauss points
Grad_at_Gauss(nodes,elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),detj,e,B,B_zero,F,debugfile);
//the material coefficients are defined inside of the function beneath
Elasticity_at_Gauss(D,F,S,debugfile,current_lambda); //for the time being this function is independent of the space
//Contribution of the geometric stiffness 
//K_gauss_geometrical(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,Kgeo,detj,_currentelement,_debugfile);//TO BE FINISHED 





calc_internal_M1(nodes,elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,F,S,internal_M1,detj,e,debugfile,l,coeff);
calc_internal_M2(nodes,elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,F,S,internal_M2,detj,e,debugfile,l,current_layup);
}



//m1 
for(int c=0;c<8;c++)
{
if(nodes[elements[e].connid[c]].neumann_bc[0]==1)
incremental_force[elements[e].connid[c]*4]=incremental_force[elements[e].connid[c]*4]+internal_M2[c*4]; 
if(nodes[elements[e].connid[c]].neumann_bc[1]==1)
incremental_force[elements[e].connid[c]*4+1]=incremental_force[elements[e].connid[c]*4+1]+internal_M2[c*4+1];
if(nodes[elements[e].connid[c]].neumann_bc[2]==1)
incremental_force[elements[e].connid[c]*4+2]=incremental_force[elements[e].connid[c]*4+2]+internal_M2[c*4+2]; 
if(nodes[elements[e].connid[c]].neumann_bc[3]==1)
incremental_force[elements[e].connid[c]*4+3]=incremental_force[elements[e].connid[c]*4+3]+internal_M2[c*4+3];
}




/*
if(e==0)
{

debugfile<<"M1internal"<<endl;
for(int c=0;c<8;c++)
debugfile<<internal_M1[c*4+0]<<'\t'<<internal_M1[c*4+1]<<'\t'<<internal_M1[c*4+2]<<'\t'<<internal_M1[c*4+3]<<endl;
debugfile<<"M2internal"<<endl;
for(int c=0;c<8;c++)
debugfile<<internal_M2[c*4+0]<<'\t'<<internal_M2[c*4+1]<<'\t'<<internal_M2[c*4+2]<<'\t'<<internal_M2[c*4+3]<<endl;
debugfile<<"Winternal"<<endl;
for(int c=0;c<8;c++)
debugfile<<Bt_dot_S[c*4+0]<<'\t'<<Bt_dot_S[c*4+1]<<'\t'<<Bt_dot_S[c*4+2]<<'\t'<<Bt_dot_S[c*4+3]<<endl;
}
*/
}


//for(int i=0;i<4*numberofnodes;i++)
//{
//if(nodes[i/4].neumann_bc[i%4]==0)
//fi[i]=-internalforcevector[i];
//else
//fi[i]=incremental_force[i]-internalforcevector[i];
//fi[i]=incremental_force[i];
//}
//fi[i]=incremental_force[i]-internalforcevector[i];
/*
debugfile<<"fi"<<endl;
for(int i=0;i<4*numberofnodes;i++)
if(nodes[i/4].neumann_bc[i%4]==1)
debugfile<<fi[i]<<'\t'; 
debugfile<<endl; 
*/
/*
debugfile<<"du"<<endl;
for(int i=0;i<4*numberofnodes;i++)
debugfile<<nodes[i/4].u_lamb[i%4]<<'\t'; 
debugfile<<endl; 
*/
/*
debugfile<<"internalforce"<<endl;
for(int i=0;i<3*numberofnodes;i++)
//if(nodes[i/3].neumann_bc[2]==1 && i%3==2)
debugfile<<internalforcevector[i]<<'\t'; 
debugfile<<endl; 
*/
/*
debugfile<<"externalforce"<<endl;
for(int i=0;i<3*numberofnodes;i++)
if(nodes[i/3].neumann_bc[2]==1 && i%3==2)
debugfile<<incrementalforcevector[i]<<'\t'; 
debugfile<<endl; 

debugfile<<"fi"<<endl;
for(int i=0;i<3*numberofnodes;i++)
if(nodes[i/3].neumann_bc[2]==1 && i%3==2)
debugfile<<fi[i]<<'\t'; 
debugfile<<endl; 
*/
}; 







void update_fi(Node nodes[], Element elements[], double fi[], double numberofloadsteps, double l, ofstream& debugfile,int numberofnodes,int numberofelements, double incremental_force[])
{





double beta_2=0.25;
double beta_1=0.5;



double coeff;
if(l>0.0)
coeff=1.0*beta_2/beta_1/l;
else if(l==0.0)
coeff=1.0;




double* internalforcevector=new double[4*numberofnodes]; 
for(int i=0;i<4*numberofnodes;i++)
internalforcevector[i]=0.0; 


double Bt_dot_S[32];
double internal_M1[32];
double internal_M2[32];



for(int e=0;e<numberofelements;e++)
{
for(int i=0;i<32;i++)
{
Bt_dot_S[i]=0.0; 
internal_M1[i]=0.0;
internal_M2[i]=0.0;
}

int gauss_counter=0;


for(double g3=-1.0;g3<2.0;g3=g3+2.0)
for(double g2=-1.0;g2<2.0;g2=g2+2.0)
for(double g1=-1.0;g1<2.0;g1=g1+2.0)
{



double current_layup; 
current_layup=elements[e].layup[gauss_counter];
gauss_counter++;

//local declerations 
double detj=0.0;
double B[6][32];//strain displacement matrix at the gauss point
double B_zero[6][32];//initial strain displacement matrix at the gauss point
double F[9];//deformation gradient at the gauss point 
double D[6][6];//elasticity matrix at the gauss point
double Bt[32][6];
double S[6];
double N[8]; //Shape Functions at the current stage, necessary for the interpolation of the lambda at the gauss location 
double current_lambda=0.0;
double ksi=g1*1.0/sqrt(3);
double eta=g2*1.0/sqrt(3);
double zeta=g3*1.0/sqrt(3);


//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N[shape_counter]=0.125*(1+sign_ksi*ksi)*(1+sign_eta*eta)*(1+sign_zeta*zeta);
shape_counter++;
}
double temp=N[2];
N[2]=N[3];
N[3]=temp;
temp=N[6];
N[6]=N[7];
N[7]=temp;
//end definition shape value vector

//interpolation of the current lambda. It is assummed that the reference pore pressure is to be taken here 
for(int s=0;s<8;s++)
current_lambda=current_lambda+N[s]*nodes[elements[e].connid[s]].u_lamb[3];



//computation of the noonlinear strain-displacement matrix and the deformation gradient at the gauss points
Grad_at_Gauss(nodes,elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),detj,e,B,B_zero,F,debugfile);
//the material coefficients are defined inside of the function beneath
Elasticity_at_Gauss(D,F,S,debugfile,current_lambda); //for the time being this function is independent of the space
//Contribution of the geometric stiffness 
//K_gauss_geometrical(_nodes,_elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,Kgeo,detj,_currentelement,_debugfile);//TO BE FINISHED 




for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
Bt[i][j]=B[j][i]-B_zero[j][i];

//this was open previously
//if(l==0.0)
for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
Bt[i][j]=B[j][i];



for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
Bt_dot_S[i]=Bt_dot_S[i]+Bt[i][j]*S[j];//*detj; 
/*
if(e==0)
{
debugfile<<"Something "<<e<<endl;
for(int j=0;j<24;j++)
debugfile<<Bt_dot_S[j]<<setw(10)<<endl;
}
*/


calc_internal_M1(nodes,elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,F,S,internal_M1,detj,e,debugfile,l,coeff);
calc_internal_M2(nodes,elements,g1*1.0/sqrt(3),g2*1.0/sqrt(3),g3*1.0/sqrt(3),B,D,F,S,internal_M2,detj,e,debugfile,l,current_layup);

}




for(int c=0;c<8;c++)
{
//if(nodes[elements[e].connid[c]].neumann_bc[0]==0)
internalforcevector[elements[e].connid[c]*4]=internalforcevector[elements[e].connid[c]*4]+Bt_dot_S[c*4]*coeff; 
//if(nodes[elements[e].connid[c]].neumann_bc[1]==0)
internalforcevector[elements[e].connid[c]*4+1]=internalforcevector[elements[e].connid[c]*4+1]+Bt_dot_S[c*4+1]*coeff;
//if(nodes[elements[e].connid[c]].neumann_bc[2]==0)
internalforcevector[elements[e].connid[c]*4+2]=internalforcevector[elements[e].connid[c]*4+2]+Bt_dot_S[c*4+2]*coeff; 
//if(nodes[elements[e].connid[c]].neumann_bc[3]==0)
internalforcevector[elements[e].connid[c]*4+3]=internalforcevector[elements[e].connid[c]*4+3]+Bt_dot_S[c*4+3]*coeff;
}


//m1 
for(int c=0;c<8;c++)
{
//if(nodes[elements[e].connid[c]].neumann_bc[0]==0)
internalforcevector[elements[e].connid[c]*4]=internalforcevector[elements[e].connid[c]*4]+internal_M1[c*4]; 
//if(nodes[elements[e].connid[c]].neumann_bc[1]==0)
internalforcevector[elements[e].connid[c]*4+1]=internalforcevector[elements[e].connid[c]*4+1]+internal_M1[c*4+1];
//if(nodes[elements[e].connid[c]].neumann_bc[2]==0)
internalforcevector[elements[e].connid[c]*4+2]=internalforcevector[elements[e].connid[c]*4+2]+internal_M1[c*4+2]; 
//if(nodes[elements[e].connid[c]].neumann_bc[3]==0)
internalforcevector[elements[e].connid[c]*4+3]=internalforcevector[elements[e].connid[c]*4+3]+internal_M1[c*4+3];
}

//m2 
for(int c=0;c<8;c++)
{
if(nodes[elements[e].connid[c]].neumann_bc[0]==0)
internalforcevector[elements[e].connid[c]*4]=internalforcevector[elements[e].connid[c]*4]+internal_M2[c*4]; 
if(nodes[elements[e].connid[c]].neumann_bc[1]==0)
internalforcevector[elements[e].connid[c]*4+1]=internalforcevector[elements[e].connid[c]*4+1]+internal_M2[c*4+1];
if(nodes[elements[e].connid[c]].neumann_bc[2]==0)
internalforcevector[elements[e].connid[c]*4+2]=internalforcevector[elements[e].connid[c]*4+2]+internal_M2[c*4+2]; 
if(nodes[elements[e].connid[c]].neumann_bc[3]==0)
internalforcevector[elements[e].connid[c]*4+3]=internalforcevector[elements[e].connid[c]*4+3]+internal_M2[c*4+3];
}


/*
if(e==0)
{

debugfile<<"M1internal"<<endl;
for(int c=0;c<8;c++)
debugfile<<internal_M1[c*4+0]<<'\t'<<internal_M1[c*4+1]<<'\t'<<internal_M1[c*4+2]<<'\t'<<internal_M1[c*4+3]<<endl;
debugfile<<"M2internal"<<endl;
for(int c=0;c<8;c++)
debugfile<<internal_M2[c*4+0]<<'\t'<<internal_M2[c*4+1]<<'\t'<<internal_M2[c*4+2]<<'\t'<<internal_M2[c*4+3]<<endl;
debugfile<<"Winternal"<<endl;
for(int c=0;c<8;c++)
debugfile<<Bt_dot_S[c*4+0]<<'\t'<<Bt_dot_S[c*4+1]<<'\t'<<Bt_dot_S[c*4+2]<<'\t'<<Bt_dot_S[c*4+3]<<endl;
}
*/
}


for(int i=0;i<4*numberofnodes;i++)
{
//if(nodes[i/4].neumann_bc[i%4]==0)
fi[i]=-internalforcevector[i];
//else
//fi[i]=incremental_force[i]-internalforcevector[i];
//fi[i]=incremental_force[i];
}
//fi[i]=incremental_force[i]-internalforcevector[i];

debugfile<<"fi"<<endl;
for(int i=0;i<4*numberofnodes;i++)
if(nodes[i/4].neumann_bc[i%4]==1)
debugfile<<fi[i]<<'\t'; 
debugfile<<endl; 

/*
debugfile<<"du"<<endl;
for(int i=0;i<4*numberofnodes;i++)
debugfile<<nodes[i/4].u_lamb[i%4]<<'\t'; 
debugfile<<endl; 
*/
/*
debugfile<<"internalforce"<<endl;
for(int i=0;i<3*numberofnodes;i++)
//if(nodes[i/3].neumann_bc[2]==1 && i%3==2)
debugfile<<internalforcevector[i]<<'\t'; 
debugfile<<endl; 
*/
/*
debugfile<<"externalforce"<<endl;
for(int i=0;i<3*numberofnodes;i++)
if(nodes[i/3].neumann_bc[2]==1 && i%3==2)
debugfile<<incrementalforcevector[i]<<'\t'; 
debugfile<<endl; 

debugfile<<"fi"<<endl;
for(int i=0;i<3*numberofnodes;i++)
if(nodes[i/3].neumann_bc[2]==1 && i%3==2)
debugfile<<fi[i]<<'\t'; 
debugfile<<endl; 
*/
}; 




void calc_K_Wlambdau_at_gauss(Node _nodes[],
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
                              double current_lambda)
{




double beta_2=0.25;
double beta_1=0.5;
double Ftens[3][3];
double FtensT[3][3];
double Cauchytens[3][3];
double CauchytensInv[3][3];
double Cauchy_rate[6];
double CauchyInv[6];
double dNdX[24]; 
double detF;
double nodal_u_only[32];
double N_lambda[32];

//can: startof:=these values are for the geometric stiffnes similar part and its convergence 
double K_M1uJ_temp[24][24];
double K_M1uJ_temp_large[32][32];
double dNdksi[24]; 
double inversejacob[24][24]; 
double jacob[3][3]; 
//can: endof:=these values are for the geometric stiffnes similar part and its convergence 




Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];

transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchytens); 
inverse(Cauchytens,CauchytensInv);
detF=det(Ftens);




CauchyInv[0]=CauchytensInv[0][0];
CauchyInv[1]=CauchytensInv[1][1];
CauchyInv[2]=CauchytensInv[2][2];
CauchyInv[3]=CauchytensInv[0][1];
CauchyInv[4]=CauchytensInv[1][2];
CauchyInv[5]=CauchytensInv[2][0];

/*
Cauchy_rate[0]=1.0*_coeff*(Cauchytens[0][0]-1.0);
Cauchy_rate[1]=1.0*_coeff*(Cauchytens[1][1]-1.0);
Cauchy_rate[2]=1.0*_coeff*(Cauchytens[2][2]-1.0);
Cauchy_rate[3]=1.0*2.0*_coeff*Cauchytens[0][1];
Cauchy_rate[4]=1.0*2.0*_coeff*Cauchytens[1][2];
Cauchy_rate[5]=1.0*2.0*_coeff*Cauchytens[2][0];
*/

Cauchy_rate[0]=1.0*_coeff*(Cauchytens[0][0]-1.0);
Cauchy_rate[1]=1.0*_coeff*(Cauchytens[1][1]-1.0);
Cauchy_rate[2]=1.0*_coeff*(Cauchytens[2][2]-1.0);
Cauchy_rate[3]=1.0*_coeff*Cauchytens[0][1];
Cauchy_rate[4]=1.0*_coeff*Cauchytens[1][2];
Cauchy_rate[5]=1.0*_coeff*Cauchytens[2][0];


for(int n=0;n<32;n++)
nodal_u_only[n]=0.0;

//there might be a problem here !!!
for(int n=0;n<24;n++)
nodal_u_only[n+n/3]=_nodes[_elements[_currentelement].connid[n/3]].u_lamb[n%4];


//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N_lambda[shape_counter+0]=0.0;
N_lambda[shape_counter+1]=0.0;
N_lambda[shape_counter+2]=0.0;
N_lambda[shape_counter+3]=0.125*(1+sign_ksi*_ksi)*(1+sign_eta*_eta)*(1+sign_zeta*_zeta);
shape_counter=shape_counter+4;
}
/*
for(int n=0;n<4;n++)
{
double temp=N_lambda[2*4+n];
N_lambda[2*4+n]=N_lambda[3*4+n];
N_lambda[3*4+n]=temp;
temp=N_lambda[6*4+n];
N_lambda[6*4+n]=N_lambda[7*4+n];
N_lambda[7*4+n]=temp;
}
*/
double temp=N_lambda[11];
N_lambda[11]=N_lambda[15];
N_lambda[15]=temp;
temp=N_lambda[27];
N_lambda[27]=N_lambda[31];
N_lambda[31]=temp;
//end definition shape value vector






double somecoeff=0.0;
for(int n=0;n<6;n++)
somecoeff=somecoeff+Cauchy_rate[n]*CauchyInv[n];




//computation of the third part 
dNdksi[0]=-0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[1]=-0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[2]=-0.125*(1.0-_ksi)*(1.0-_eta);


dNdksi[3]=0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[4]=-0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[5]=-0.125*(1.0+_ksi)*(1.0-_eta);


dNdksi[6]=0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[7]=0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[8]=-0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[9]=-0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[10]=0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[11]=-0.125*(1.0-_ksi)*(1.0+_eta);

dNdksi[12]=-0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[13]=-0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[14]=0.125*(1.0-_ksi)*(1.0-_eta);

dNdksi[15]=0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[16]=-0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[17]=0.125*(1.0+_ksi)*(1.0-_eta);

dNdksi[18]=0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[19]=0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[20]=0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[21]=-0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[22]=0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[23]=0.125*(1.0-_ksi)*(1.0+_eta);




for(int i=0;i<24;i++)
{
for(int j=0;j<24;j++)
inversejacob[i][j]=0.0; 
dNdX[i]=0.0; 
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0;


for(int i=0;i<3;i++) // row label 
for(int j=0;j<8;j++) 
{
jacob[i][0]=jacob[i][0]+dNdksi[j*3+i]*(_nodes[_elements[_currentelement].connid[j]].X_lamb[0]+_nodes[_elements[_currentelement].connid[j]].u_lamb[0]);
jacob[i][1]=jacob[i][1]+dNdksi[j*3+i]*(_nodes[_elements[_currentelement].connid[j]].X_lamb[1]+_nodes[_elements[_currentelement].connid[j]].u_lamb[1]);
jacob[i][2]=jacob[i][2]+dNdksi[j*3+i]*(_nodes[_elements[_currentelement].connid[j]].X_lamb[2]+_nodes[_elements[_currentelement].connid[j]].u_lamb[2]);
}

detj=jacob[0][0]*(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])-
     jacob[0][1]*(jacob[1][0]*jacob[2][2]-jacob[2][0]*jacob[1][2])+
     jacob[0][2]*(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1]);


inversejacob[0][0]=inversejacob[3][3]=inversejacob[6][6]=inversejacob[9][9]=inversejacob[12][12]=inversejacob[15][15]=inversejacob[18][18]=inversejacob[21][21]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
inversejacob[0][1]=inversejacob[3][4]=inversejacob[6][7]=inversejacob[9][10]=inversejacob[12][13]=inversejacob[15][16]=inversejacob[18][19]=inversejacob[21][22]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
inversejacob[0][2]=inversejacob[3][5]=inversejacob[6][8]=inversejacob[9][11]=inversejacob[12][14]=inversejacob[15][17]=inversejacob[18][20]=inversejacob[21][23]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj; 

inversejacob[1][0]=inversejacob[4][3]=inversejacob[7][6]=inversejacob[10][9]=inversejacob[13][12]=inversejacob[16][15]=inversejacob[19][18]=inversejacob[22][21]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
inversejacob[1][1]=inversejacob[4][4]=inversejacob[7][7]=inversejacob[10][10]=inversejacob[13][13]=inversejacob[16][16]=inversejacob[19][19]=inversejacob[22][22]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
inversejacob[1][2]=inversejacob[4][5]=inversejacob[7][8]=inversejacob[10][11]=inversejacob[13][14]=inversejacob[16][17]=inversejacob[19][20]=inversejacob[22][23]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;

inversejacob[2][0]=inversejacob[5][3]=inversejacob[8][6]=inversejacob[11][9]=inversejacob[14][12]=inversejacob[17][15]=inversejacob[20][18]=inversejacob[23][21]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;
inversejacob[2][1]=inversejacob[5][4]=inversejacob[8][7]=inversejacob[11][10]=inversejacob[14][13]=inversejacob[17][16]=inversejacob[20][19]=inversejacob[23][22]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;
inversejacob[2][2]=inversejacob[5][5]=inversejacob[8][8]=inversejacob[11][11]=inversejacob[14][14]=inversejacob[17][17]=inversejacob[20][20]=inversejacob[23][23]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;


for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
dNdX[i]=dNdX[i]+inversejacob[i][j]*dNdksi[j];





double dlambda_du[32];

for(int i=0;i<32;i++)
dlambda_du[i]=0.0;




for(int i=0;i<24;i++) 
dlambda_du[i+i/3]=dNdX[i]*current_lambda;






double first_v[32]; //left node driven 




//they are all supposed to be initialized as zero 
for(int n=0;n<32;n++)
first_v[n]=0.0;



for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
first_v[i]=first_v[i]+CauchyInv[j]*B[j][i];





for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
{
_K_Wlambdau[j][i]=_K_Wlambdau[j][i]-first_v[i]*dlambda_du[j]*_coeff;
//_K_Wlambdau[i][j]=_K_Wlambdau[j][i];
}



/*
double trans[32][32];
double self[32][32];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
{
trans[j][i]=_K_Wlambdau[i][j];
self[i][j]=_K_Wlambdau[i][j];
}

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_Wlambdau[i][j]=0.5*(self[i][j]+trans[i][j]);
*/
/*
double Ktrans[32][32];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
Ktrans[j][i]=_K_Wlambdau[i][j];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_Wlambdau[i][j]=0.5*(_K_Wlambdau[i][j]+Ktrans[i][j]);
*/
//integration over the left and right nodes 

/*
//loop over the right node 
for(int L=0;L<8;L++)
{
//for the sake of easiness, there are in total 3 vectors splitted 
double first_v[32]; //left node driven 
double second_v[32]; //right node driven 
double third_v[32]; //right node driven 


//they are all supposed to be initialized as zero 
for(int n=0;n<32;n++)
first_v[n]=0.0;



for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
if(i/4==L) //only the corresponding nodes are taken 
first_v[i]=first_v[i]+CauchyInv[j]*B[j][i];






//loop over the right node
for(int M=0;M<8;M++)
{
for(int n=0;n<32;n++)
second_v[n]=0.0;

for(int n=0;n<4;n++)
second_v[M*4+n]=-0.5*N_lambda[M*4+n];




for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_Wlambda[i][j]=_K_Wlambda[i][j]+first_v[i]*second_v[j];

}
}
*/
};













void calc_K_Wlambda_at_gauss(Node _nodes[],
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
                             double _coeff)
{




double beta_2=0.25;
double beta_1=0.5;
double Ftens[3][3];
double FtensT[3][3];
double Cauchytens[3][3];
double CauchytensInv[3][3];
double Cauchy_rate[6];
double CauchyInv[6];
double dNdX[24]; 
double detF;
double nodal_u_only[32];
double N_lambda[32];

//can: startof:=these values are for the geometric stiffnes similar part and its convergence 
double K_M1uJ_temp[24][24];
double K_M1uJ_temp_large[32][32];
double dNdksi[24]; 
double inversejacob[24][24]; 
double jacob[3][3]; 
//can: endof:=these values are for the geometric stiffnes similar part and its convergence 




Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];

transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchytens); 
inverse(Cauchytens,CauchytensInv);
detF=det(Ftens);




CauchyInv[0]=CauchytensInv[0][0];
CauchyInv[1]=CauchytensInv[1][1];
CauchyInv[2]=CauchytensInv[2][2];
CauchyInv[3]=CauchytensInv[0][1];
CauchyInv[4]=CauchytensInv[1][2];
CauchyInv[5]=CauchytensInv[2][0];


Cauchy_rate[0]=1.0*_coeff*(Cauchytens[0][0]-1.0);
Cauchy_rate[1]=1.0*_coeff*(Cauchytens[1][1]-1.0);
Cauchy_rate[2]=1.0*_coeff*(Cauchytens[2][2]-1.0);
Cauchy_rate[3]=1.0*_coeff*Cauchytens[0][1];
Cauchy_rate[4]=1.0*_coeff*Cauchytens[1][2];
Cauchy_rate[5]=1.0*_coeff*Cauchytens[2][0];




for(int n=0;n<32;n++)
nodal_u_only[n]=0.0;

//there might be a problem here !!!
for(int n=0;n<24;n++)
nodal_u_only[n+n/3]=_nodes[_elements[_currentelement].connid[n/3]].u_lamb[n%4];


//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N_lambda[shape_counter+0]=0.0;
N_lambda[shape_counter+1]=0.0;
N_lambda[shape_counter+2]=0.0;
N_lambda[shape_counter+3]=0.125*(1+sign_ksi*_ksi)*(1+sign_eta*_eta)*(1+sign_zeta*_zeta);
shape_counter=shape_counter+4;
}
/*
for(int n=0;n<4;n++)
{
double temp=N_lambda[2*4+n];
N_lambda[2*4+n]=N_lambda[3*4+n];
N_lambda[3*4+n]=temp;
temp=N_lambda[6*4+n];
N_lambda[6*4+n]=N_lambda[7*4+n];
N_lambda[7*4+n]=temp;
}
*/
double temp=N_lambda[11];
N_lambda[11]=N_lambda[15];
N_lambda[15]=temp;
temp=N_lambda[27];
N_lambda[27]=N_lambda[31];
N_lambda[31]=temp;
//end definition shape value vector



double somecoeff=0.0;
for(int n=0;n<6;n++)
somecoeff=somecoeff+Cauchy_rate[n]*CauchyInv[n];




//computation of the third part 
dNdksi[0]=-0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[1]=-0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[2]=-0.125*(1.0-_ksi)*(1.0-_eta);


dNdksi[3]=0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[4]=-0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[5]=-0.125*(1.0+_ksi)*(1.0-_eta);


dNdksi[6]=0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[7]=0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[8]=-0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[9]=-0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[10]=0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[11]=-0.125*(1.0-_ksi)*(1.0+_eta);

dNdksi[12]=-0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[13]=-0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[14]=0.125*(1.0-_ksi)*(1.0-_eta);

dNdksi[15]=0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[16]=-0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[17]=0.125*(1.0+_ksi)*(1.0-_eta);

dNdksi[18]=0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[19]=0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[20]=0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[21]=-0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[22]=0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[23]=0.125*(1.0-_ksi)*(1.0+_eta);




for(int i=0;i<24;i++)
{
for(int j=0;j<24;j++)
inversejacob[i][j]=0.0; 
dNdX[i]=0.0; 
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0;


for(int i=0;i<3;i++) // row label 
for(int j=0;j<8;j++) 
{
jacob[i][0]=jacob[i][0]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[0];
jacob[i][1]=jacob[i][1]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[1];
jacob[i][2]=jacob[i][2]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[2];
}

detj=jacob[0][0]*(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])-
     jacob[0][1]*(jacob[1][0]*jacob[2][2]-jacob[2][0]*jacob[1][2])+
     jacob[0][2]*(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1]);


inversejacob[0][0]=inversejacob[3][3]=inversejacob[6][6]=inversejacob[9][9]=inversejacob[12][12]=inversejacob[15][15]=inversejacob[18][18]=inversejacob[21][21]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
inversejacob[0][1]=inversejacob[3][4]=inversejacob[6][7]=inversejacob[9][10]=inversejacob[12][13]=inversejacob[15][16]=inversejacob[18][19]=inversejacob[21][22]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
inversejacob[0][2]=inversejacob[3][5]=inversejacob[6][8]=inversejacob[9][11]=inversejacob[12][14]=inversejacob[15][17]=inversejacob[18][20]=inversejacob[21][23]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj; 

inversejacob[1][0]=inversejacob[4][3]=inversejacob[7][6]=inversejacob[10][9]=inversejacob[13][12]=inversejacob[16][15]=inversejacob[19][18]=inversejacob[22][21]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
inversejacob[1][1]=inversejacob[4][4]=inversejacob[7][7]=inversejacob[10][10]=inversejacob[13][13]=inversejacob[16][16]=inversejacob[19][19]=inversejacob[22][22]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
inversejacob[1][2]=inversejacob[4][5]=inversejacob[7][8]=inversejacob[10][11]=inversejacob[13][14]=inversejacob[16][17]=inversejacob[19][20]=inversejacob[22][23]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;

inversejacob[2][0]=inversejacob[5][3]=inversejacob[8][6]=inversejacob[11][9]=inversejacob[14][12]=inversejacob[17][15]=inversejacob[20][18]=inversejacob[23][21]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;
inversejacob[2][1]=inversejacob[5][4]=inversejacob[8][7]=inversejacob[11][10]=inversejacob[14][13]=inversejacob[17][16]=inversejacob[20][19]=inversejacob[23][22]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;
inversejacob[2][2]=inversejacob[5][5]=inversejacob[8][8]=inversejacob[11][11]=inversejacob[14][14]=inversejacob[17][17]=inversejacob[20][20]=inversejacob[23][23]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;


for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
dNdX[i]=dNdX[i]+inversejacob[i][j]*dNdksi[j];


double first_v[32]; //left node driven 
double second_v[32]; //right node driven 



//they are all supposed to be initialized as zero 
for(int n=0;n<32;n++)
first_v[n]=0.0;

for(int n=0;n<32;n++)
second_v[n]=0.0;

for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
first_v[i]=first_v[i]+CauchyInv[j]*B[j][i];



for(int n=0;n<32;n++)
second_v[n]=-0.5*N_lambda[n];


for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
{
_K_Wlambda[i][j]=_K_Wlambda[i][j]+first_v[i]*second_v[j]*_coeff;
//_K_Wlambda[i][j]=_K_Wlambda[j][i];
}
/*
for(int j=0;j<8;j++)
for(int i=0;i<32;i++)
_K_Wlambda[j*4+3][i]=_K_Wlambda[i][j*4+3];
*/
/*
double trans[32][32];
double self[32][32];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
{
trans[j][i]=_K_Wlambda[i][j];
self[i][j]=_K_Wlambda[i][j];
}

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_Wlambda[i][j]=0.5*(self[i][j]+trans[i][j]);
*/

/*
double Ktrans[32][32];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
Ktrans[j][i]=_K_Wlambda[i][j];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_Wlambda[i][j]=0.5*(_K_Wlambda[i][j]+Ktrans[i][j]);
*/

/*
for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
if(i%4==3)
_K_Wlambda[i][j]=_K_Wlambda[j][i];
*/


//integration over the left and right nodes 

/*
//loop over the right node 
for(int L=0;L<8;L++)
{
//for the sake of easiness, there are in total 3 vectors splitted 
double first_v[32]; //left node driven 
double second_v[32]; //right node driven 
double third_v[32]; //right node driven 


//they are all supposed to be initialized as zero 
for(int n=0;n<32;n++)
first_v[n]=0.0;



for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
if(i/4==L) //only the corresponding nodes are taken 
first_v[i]=first_v[i]+CauchyInv[j]*B[j][i];






//loop over the right node
for(int M=0;M<8;M++)
{
for(int n=0;n<32;n++)
second_v[n]=0.0;

for(int n=0;n<4;n++)
second_v[M*4+n]=-0.5*N_lambda[M*4+n];




for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_Wlambda[i][j]=_K_Wlambda[i][j]+first_v[i]*second_v[j];

}
}
*/
};
























void calc_K_M1uJ_at_gauss(Node _nodes[],
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
                          double _coeff)
{


double beta_2=0.25;
double beta_1=0.5;
double Ftens[3][3];
double FtensT[3][3];
double Cauchytens[3][3];
double CauchytensInv[3][3];
double Cauchy_rate[6];
double CauchyInv[6];
double dNdX[24]; 
double detF;
double nodal_u_only[32];
double N_lambda[32];

//can: startof:=these values are for the geometric stiffnes similar part and its convergence 
double K_M1uJ_temp[24][24];
double K_M1uJ_temp_large[32][32];
double dNdksi[24]; 
double inversejacob[24][24]; 
double jacob[3][3]; 
//can: endof:=these values are for the geometric stiffnes similar part and its convergence 


Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];

transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchytens); 
inverse(Cauchytens,CauchytensInv);
detF=det(Ftens);


CauchyInv[0]=CauchytensInv[0][0];
CauchyInv[1]=CauchytensInv[1][1];
CauchyInv[2]=CauchytensInv[2][2];
CauchyInv[3]=CauchytensInv[0][1];
CauchyInv[4]=CauchytensInv[1][2];
CauchyInv[5]=CauchytensInv[2][0];

/*
Cauchy_rate[0]=1.0*_coeff*(Cauchytens[0][0]-1.0);
Cauchy_rate[1]=1.0*_coeff*(Cauchytens[1][1]-1.0);
Cauchy_rate[2]=1.0*_coeff*(Cauchytens[2][2]-1.0);
Cauchy_rate[3]=1.0*2.0*_coeff*Cauchytens[0][1];
Cauchy_rate[4]=1.0*2.0*_coeff*Cauchytens[1][2];
Cauchy_rate[5]=1.0*2.0*_coeff*Cauchytens[2][0];

*/

Cauchy_rate[0]=1.0*_coeff*(Cauchytens[0][0]-1.0);
Cauchy_rate[1]=1.0*_coeff*(Cauchytens[1][1]-1.0);
Cauchy_rate[2]=1.0*_coeff*(Cauchytens[2][2]-1.0);
Cauchy_rate[3]=1.0*_coeff*Cauchytens[0][1];
Cauchy_rate[4]=1.0*_coeff*Cauchytens[1][2];
Cauchy_rate[5]=1.0*_coeff*Cauchytens[2][0];




for(int n=0;n<32;n++)
nodal_u_only[n]=0.0;

for(int n=0;n<24;n++)
nodal_u_only[n+n/3]=_nodes[_elements[_currentelement].connid[n/3]].u_lamb[n%4];


//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N_lambda[shape_counter+0]=0.0;
N_lambda[shape_counter+1]=0.0;
N_lambda[shape_counter+2]=0.0;
N_lambda[shape_counter+3]=0.125*(1+sign_ksi*_ksi)*(1+sign_eta*_eta)*(1+sign_zeta*_zeta);
shape_counter=shape_counter+4;
}

double temp=N_lambda[11];
N_lambda[11]=N_lambda[15];
N_lambda[15]=temp;
temp=N_lambda[27];
N_lambda[27]=N_lambda[31];
N_lambda[31]=temp;
/*
for(int n=0;n<4;n++)
{
double temp=N_lambda[2*4+n];
N_lambda[2*4+n]=N_lambda[3*4+n];
N_lambda[3*4+n]=temp;
temp=N_lambda[6*4+n];
N_lambda[6*4+n]=N_lambda[7*4+n];
N_lambda[7*4+n]=temp;
}
*/
//end definition shape value vector



double somecoeff=0.0;
for(int n=0;n<6;n++)
somecoeff=somecoeff+Cauchy_rate[n]*CauchyInv[n];




//computation of the third part 
dNdksi[0]=-0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[1]=-0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[2]=-0.125*(1.0-_ksi)*(1.0-_eta);


dNdksi[3]=0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[4]=-0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[5]=-0.125*(1.0+_ksi)*(1.0-_eta);


dNdksi[6]=0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[7]=0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[8]=-0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[9]=-0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[10]=0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[11]=-0.125*(1.0-_ksi)*(1.0+_eta);

dNdksi[12]=-0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[13]=-0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[14]=0.125*(1.0-_ksi)*(1.0-_eta);

dNdksi[15]=0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[16]=-0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[17]=0.125*(1.0+_ksi)*(1.0-_eta);

dNdksi[18]=0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[19]=0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[20]=0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[21]=-0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[22]=0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[23]=0.125*(1.0-_ksi)*(1.0+_eta);




for(int i=0;i<24;i++)
{
for(int j=0;j<24;j++)
inversejacob[i][j]=0.0; 
dNdX[i]=0.0; 
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0;


for(int i=0;i<3;i++) // row label 
for(int j=0;j<8;j++) 
{
jacob[i][0]=jacob[i][0]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[0];
jacob[i][1]=jacob[i][1]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[1];
jacob[i][2]=jacob[i][2]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[2];
}

detj=jacob[0][0]*(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])-
     jacob[0][1]*(jacob[1][0]*jacob[2][2]-jacob[2][0]*jacob[1][2])+
     jacob[0][2]*(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1]);


inversejacob[0][0]=inversejacob[3][3]=inversejacob[6][6]=inversejacob[9][9]=inversejacob[12][12]=inversejacob[15][15]=inversejacob[18][18]=inversejacob[21][21]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
inversejacob[0][1]=inversejacob[3][4]=inversejacob[6][7]=inversejacob[9][10]=inversejacob[12][13]=inversejacob[15][16]=inversejacob[18][19]=inversejacob[21][22]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
inversejacob[0][2]=inversejacob[3][5]=inversejacob[6][8]=inversejacob[9][11]=inversejacob[12][14]=inversejacob[15][17]=inversejacob[18][20]=inversejacob[21][23]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj; 

inversejacob[1][0]=inversejacob[4][3]=inversejacob[7][6]=inversejacob[10][9]=inversejacob[13][12]=inversejacob[16][15]=inversejacob[19][18]=inversejacob[22][21]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
inversejacob[1][1]=inversejacob[4][4]=inversejacob[7][7]=inversejacob[10][10]=inversejacob[13][13]=inversejacob[16][16]=inversejacob[19][19]=inversejacob[22][22]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
inversejacob[1][2]=inversejacob[4][5]=inversejacob[7][8]=inversejacob[10][11]=inversejacob[13][14]=inversejacob[16][17]=inversejacob[19][20]=inversejacob[22][23]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;

inversejacob[2][0]=inversejacob[5][3]=inversejacob[8][6]=inversejacob[11][9]=inversejacob[14][12]=inversejacob[17][15]=inversejacob[20][18]=inversejacob[23][21]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;
inversejacob[2][1]=inversejacob[5][4]=inversejacob[8][7]=inversejacob[11][10]=inversejacob[14][13]=inversejacob[17][16]=inversejacob[20][19]=inversejacob[23][22]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;
inversejacob[2][2]=inversejacob[5][5]=inversejacob[8][8]=inversejacob[11][11]=inversejacob[14][14]=inversejacob[17][17]=inversejacob[20][20]=inversejacob[23][23]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;


for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
dNdX[i]=dNdX[i]+inversejacob[i][j]*dNdksi[j];



//for the sake of easiness, there are in total 3 vectors splitted 
double first_v[32]; //left node driven 
double second_v[32]; //right node driven 




for(int n=0;n<32;n++)
{
first_v[n]=0.0;
second_v[n]=0.0;
}



for(int n=0;n<32;n++)
first_v[n]=0.5*detF*somecoeff*N_lambda[n]; //it was different before 

for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
second_v[i]=second_v[i]+CauchyInv[j]*B[j][i];


for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
{
_K_M1uJ[j][i]=_K_M1uJ[j][i]+first_v[i]*second_v[j];
//_K_M1uJ[i][j]=_K_M1uJ[j][i];
}
/*
double Ktrans[32][32];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
Ktrans[j][i]=_K_M1uJ[i][j];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_M1uJ[i][j]=0.5*(_K_M1uJ[i][j]+Ktrans[i][j]);
*/

/*
for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
if(i%4==3)
_K_M1uJ[j][i]=_K_M1uJ[i][j];
*/

//integration over the left and right nodes 

/*
//loop over the right node 
for(int L=0;L<8;L++)
{
//for the sake of easiness, there are in total 3 vectors splitted 
double first_v[32]; //left node driven 
double second_v[32]; //right node driven 
double third_v[32]; //right node driven 

//they are all supposed to be initialized as zero 
for(int n=0;n<32;n++)
first_v[n]=0.0;


for(int n=0;n<4;n++)
first_v[L*4+n]=-0.5/detF*somecoeff*N_lambda[L*4+n];
//loop over the right node
for(int M=0;M<8;M++)
{
for(int n=0;n<32;n++)
{
second_v[n]=0.0;
third_v[n]=0.0;
}

for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
if(i/4==M) //only the corresponding nodes are taken 
second_v[i]=second_v[i]+CauchyInv[j]*B[j][i];

/*
for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
K_M1uJ_temp[i][j]=0.0;

// THIS HAS TO BE CHANGED MOST PROBABLY
for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
if(i/3==M && j/3==M) //only the dyadic product of the right node is done 
{
if(i%3==0 && j%3==0)
K_M1uJ_temp[i][j]=K_M1uJ_temp[i][j]+
(dNdX[i]*CauchyInv[0]*dNdX[j]
+dNdX[i+1]*CauchyInv[1]*dNdX[j+1]
+dNdX[i+2]*CauchyInv[2]*dNdX[j+2]
+dNdX[i]*CauchyInv[3]*dNdX[j+1]
+dNdX[i+1]*CauchyInv[3]*dNdX[j]
+dNdX[i+1]*CauchyInv[4]*dNdX[j+2]
+dNdX[i+2]*CauchyInv[4]*dNdX[j+1]
+dNdX[i]*CauchyInv[5]*dNdX[j+2]
+dNdX[i+2]*CauchyInv[5]*dNdX[j]);//*detj; 
else if(i%3!=0 && j>0)
K_M1uJ_temp[i][j]=K_M1uJ_temp[i-1][j-1]; 
}



for(int i=0;i<24;i++)
for(int k=0;k<24;k++)
K_M1uJ_temp_large[i+i/3][k+k/3]=K_M1uJ_temp[i][k];



//according to my nodes the third addend is zero 

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
if(j/4==M)//take only the right node 
third_v[i]=third_v[i]+nodal_u_only[j]*K_M1uJ_temp_large[j][i];
*/
/*
for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_M1uJ[i][j]=_K_M1uJ[i][j]+first_v[i]*(second_v[j]+third_v[j]);

}
}
*/
};
































void calc_K_M1uGeo_at_gauss(Node _nodes[],
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
			    double _coeff)
{


double beta_2=0.25;
double beta_1=0.5;
double Ftens[3][3];
double FtensT[3][3];
double Cauchytens[3][3];
double CauchytensInv[3][3];
double Cauchy_rate[6];
double CauchyInv[6];
double dNdX[24]; 
double detF;
double nodal_u_only[32];
double N_lambda[32];

//can: startof:=these values are for the geometric stiffnes similar part and its convergence 
double K_M1uGeo_temp[24][24];
double K_M1uGeo_temp_large[32][32];
double dNdksi[24]; 
double inversejacob[24][24]; 
double jacob[3][3]; 
//can: endof:=these values are for the geometric stiffnes similar part and its convergence 


Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];

transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchytens); 
inverse(Cauchytens,CauchytensInv);
detF=det(Ftens);


CauchyInv[0]=CauchytensInv[0][0];
CauchyInv[1]=CauchytensInv[1][1];
CauchyInv[2]=CauchytensInv[2][2];
CauchyInv[3]=CauchytensInv[0][1];
CauchyInv[4]=CauchytensInv[1][2];
CauchyInv[5]=CauchytensInv[2][0];




for(int n=0;n<32;n++)
nodal_u_only[n]=0.0;

for(int n=0;n<24;n++)
nodal_u_only[n+n/3]=_nodes[_elements[_currentelement].connid[n/3]].u_lamb[n%4];


//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N_lambda[shape_counter+0]=0.0;
N_lambda[shape_counter+1]=0.0;
N_lambda[shape_counter+2]=0.0;
N_lambda[shape_counter+3]=0.125*(1+sign_ksi*_ksi)*(1+sign_eta*_eta)*(1+sign_zeta*_zeta);
shape_counter=shape_counter+4;
}
/*
for(int n=0;n<4;n++)
{
double temp=N_lambda[2*4+n];
N_lambda[2*4+n]=N_lambda[3*4+n];
N_lambda[3*4+n]=temp;
temp=N_lambda[6*4+n];
N_lambda[6*4+n]=N_lambda[7*4+n];
N_lambda[7*4+n]=temp;
}
*/
double temp=N_lambda[11];
N_lambda[11]=N_lambda[15];
N_lambda[15]=temp;
temp=N_lambda[27];
N_lambda[27]=N_lambda[31];
N_lambda[31]=temp;
//end definition shape value vector







//computation of the third part 
dNdksi[0]=-0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[1]=-0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[2]=-0.125*(1.0-_ksi)*(1.0-_eta);


dNdksi[3]=0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[4]=-0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[5]=-0.125*(1.0+_ksi)*(1.0-_eta);


dNdksi[6]=0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[7]=0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[8]=-0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[9]=-0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[10]=0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[11]=-0.125*(1.0-_ksi)*(1.0+_eta);

dNdksi[12]=-0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[13]=-0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[14]=0.125*(1.0-_ksi)*(1.0-_eta);

dNdksi[15]=0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[16]=-0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[17]=0.125*(1.0+_ksi)*(1.0-_eta);

dNdksi[18]=0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[19]=0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[20]=0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[21]=-0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[22]=0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[23]=0.125*(1.0-_ksi)*(1.0+_eta);




for(int i=0;i<24;i++)
{
for(int j=0;j<24;j++)
inversejacob[i][j]=0.0; 
dNdX[i]=0.0; 
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0;


for(int i=0;i<3;i++) // row label 
for(int j=0;j<8;j++) 
{
jacob[i][0]=jacob[i][0]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[0];
jacob[i][1]=jacob[i][1]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[1];
jacob[i][2]=jacob[i][2]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[2];
}

detj=jacob[0][0]*(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])-
     jacob[0][1]*(jacob[1][0]*jacob[2][2]-jacob[2][0]*jacob[1][2])+
     jacob[0][2]*(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1]);


inversejacob[0][0]=inversejacob[3][3]=inversejacob[6][6]=inversejacob[9][9]=inversejacob[12][12]=inversejacob[15][15]=inversejacob[18][18]=inversejacob[21][21]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
inversejacob[0][1]=inversejacob[3][4]=inversejacob[6][7]=inversejacob[9][10]=inversejacob[12][13]=inversejacob[15][16]=inversejacob[18][19]=inversejacob[21][22]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
inversejacob[0][2]=inversejacob[3][5]=inversejacob[6][8]=inversejacob[9][11]=inversejacob[12][14]=inversejacob[15][17]=inversejacob[18][20]=inversejacob[21][23]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj; 

inversejacob[1][0]=inversejacob[4][3]=inversejacob[7][6]=inversejacob[10][9]=inversejacob[13][12]=inversejacob[16][15]=inversejacob[19][18]=inversejacob[22][21]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
inversejacob[1][1]=inversejacob[4][4]=inversejacob[7][7]=inversejacob[10][10]=inversejacob[13][13]=inversejacob[16][16]=inversejacob[19][19]=inversejacob[22][22]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
inversejacob[1][2]=inversejacob[4][5]=inversejacob[7][8]=inversejacob[10][11]=inversejacob[13][14]=inversejacob[16][17]=inversejacob[19][20]=inversejacob[22][23]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;

inversejacob[2][0]=inversejacob[5][3]=inversejacob[8][6]=inversejacob[11][9]=inversejacob[14][12]=inversejacob[17][15]=inversejacob[20][18]=inversejacob[23][21]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;
inversejacob[2][1]=inversejacob[5][4]=inversejacob[8][7]=inversejacob[11][10]=inversejacob[14][13]=inversejacob[17][16]=inversejacob[20][19]=inversejacob[23][22]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;
inversejacob[2][2]=inversejacob[5][5]=inversejacob[8][8]=inversejacob[11][11]=inversejacob[14][14]=inversejacob[17][17]=inversejacob[20][20]=inversejacob[23][23]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;


for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
dNdX[i]=dNdX[i]+inversejacob[i][j]*dNdksi[j];





double first_v[32]; //left node driven 
double second_v[32]; //right node driven 
double third_v[32]; //right node driven 
for(int n=0;n<32;n++)
first_v[n]=0.0;
for(int n=0;n<32;n++)
first_v[n]=_coeff*N_lambda[n];
for(int n=0;n<32;n++)
{
second_v[n]=0.0;
third_v[n]=0.0;
}
for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
second_v[i]=second_v[i]+CauchyInv[j]*B[j][i];


for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
{
_K_M1uGeo[j][i]=_K_M1uGeo[j][i]+first_v[i]*second_v[j];
//_K_M1uGeo[i][j]=_K_M1uGeo[j][i];
}
/*
double Ktrans[32][32];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
Ktrans[j][i]=_K_M1uGeo[i][j];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_M1uGeo[i][j]=0.5*(_K_M1uGeo[i][j]+Ktrans[i][j]);
*/

/*
for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
if(i%4==3)
_K_M1uGeo[j][i]=_K_M1uGeo[i][j];
*/

//integration over the left and right nodes 
/*

//loop over the right node 
for(int L=0;L<8;L++)
{
//for the sake of easiness, there are in total 3 vectors splitted 
double first_v[32]; //left node driven 
double second_v[32]; //right node driven 
double third_v[32]; //right node driven 

//they are all supposed to be initialized as zero 
for(int n=0;n<32;n++)
first_v[n]=0.0;



for(int n=0;n<4;n++)
first_v[L*4+n]=0.5/detF*beta_2/beta_1/_t*N_lambda[L*4+n];
//loop over the right node
for(int M=0;M<8;M++)
{
for(int n=0;n<32;n++)
{
second_v[n]=0.0;
third_v[n]=0.0;
}

for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
if(i/4==M) //only the corresponding nodes are taken 
second_v[i]=second_v[i]+CauchyInv[j]*B[j][i];

/*
for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
K_M1uGeo_temp[i][j]=0.0;

// THIS HAS TO BE CHANGED MOST PROBABLY
for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
if(i/3==M && j/3==M) //only the dyadic product of the right node is done 
{
if(i%3==0 && j%3==0)
K_M1uGeo_temp[i][j]=K_M1uGeo_temp[i][j]+
(dNdX[i]*CauchyInv[0]*dNdX[j]
+dNdX[i+1]*CauchyInv[1]*dNdX[j+1]
+dNdX[i+2]*CauchyInv[2]*dNdX[j+2]
+dNdX[i]*CauchyInv[3]*dNdX[j+1]
+dNdX[i+1]*CauchyInv[3]*dNdX[j]
+dNdX[i+1]*CauchyInv[4]*dNdX[j+2]
+dNdX[i+2]*CauchyInv[4]*dNdX[j+1]
+dNdX[i]*CauchyInv[5]*dNdX[j+2]
+dNdX[i+2]*CauchyInv[5]*dNdX[j]);//*detj; 
else if(i%3!=0 && j>0)
K_M1uGeo_temp[i][j]=K_M1uGeo_temp[i-1][j-1]; 
}



for(int i=0;i<24;i++)
for(int k=0;k<24;k++)
K_M1uGeo_temp_large[i+i/3][k+k/3]=K_M1uGeo_temp[i][k];
*/
//according to my notes the third addend is zero 
/*
for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
if(j/4==M)//take only the right node 
third_v[i]=third_v[i]+nodal_u_only[j]*K_M1uJ_temp_large[j][i];
*/
/*
for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_M1uGeo[i][j]=_K_M1uGeo[i][j]+first_v[i]*(second_v[j]+third_v[j]);

}
}
*/
};



void Pore_Elasticity_at_Gauss(double D_pore[][6],  
	    	              double F[], //deformation gradient 
			      ofstream &_debugfile,
			      double lambda)
{
double Ftens[3][3];
double FtensT[3][3];
double Cauchy[3][3]; //Cauchy-Green Strain Tensor 
double CauchyInv[3][3]; //Cauchy-Green Strain Tensor 
double C[3][3][3][3]; //bulk constitutive tensor 

Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];

transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchy); 
inverse(Cauchy,CauchyInv);




//effect of fluid pore pressure
for(int I=0;I<3;I++)
for(int J=0;J<3;J++)
for(int K=0;K<3;K++)
for(int L=0;L<3;L++)
C[I][J][K][L]=lambda*(CauchyInv[I][K]*CauchyInv[J][L]+CauchyInv[I][L]*CauchyInv[J][K]);





D_pore[0][0]=C[0][0][0][0];
D_pore[0][1]=C[0][0][1][1];
D_pore[0][2]=C[0][0][2][2];
D_pore[0][3]=C[0][0][0][1];
D_pore[0][4]=C[0][0][1][2];
D_pore[0][5]=C[0][0][2][0];
 
D_pore[1][0]=C[1][1][0][0];
D_pore[1][1]=C[1][1][1][1];
D_pore[1][2]=C[1][1][2][2];
D_pore[1][3]=C[1][1][0][1];
D_pore[1][4]=C[1][1][1][2];
D_pore[1][5]=C[1][1][2][0];

D_pore[2][0]=C[2][2][0][0];
D_pore[2][1]=C[2][2][1][1];
D_pore[2][2]=C[2][2][2][2];
D_pore[2][3]=C[2][2][0][1];
D_pore[2][4]=C[2][2][1][2];
D_pore[2][5]=C[2][2][2][0];

D_pore[3][0]=C[0][1][0][0];
D_pore[3][1]=C[0][1][1][1];
D_pore[3][2]=C[0][1][2][2];
D_pore[3][3]=C[0][1][0][1];
D_pore[3][4]=C[0][1][1][2];
D_pore[3][5]=C[0][1][2][0];

D_pore[4][0]=C[1][2][0][0];
D_pore[4][1]=C[1][2][1][1];
D_pore[4][2]=C[1][2][2][2];
D_pore[4][3]=C[1][2][0][1];
D_pore[4][4]=C[1][2][1][2];
D_pore[4][5]=C[1][2][2][0];

D_pore[5][0]=C[2][0][0][0];
D_pore[5][1]=C[2][0][1][1];
D_pore[5][2]=C[2][0][2][2];
D_pore[5][3]=C[2][0][0][1];
D_pore[5][4]=C[2][0][1][2];
D_pore[5][5]=C[2][0][2][0];



};






void calc_K_M1uMat_at_gauss(Node _nodes[],
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
			    double _coeff)
{


double beta_2=0.25;
double beta_1=0.5;
double Ftens[3][3];
double FtensT[3][3];
double Cauchytens[3][3];
double Cauchy_rate[6];
double dNdX[24]; 
double detF;
double N_lambda[32];
double D_pore[6][6];

//can: startof:=these values are for the geometric stiffnes similar part and its convergence 
double dNdksi[24]; 
double inversejacob[24][24]; 
double jacob[3][3]; 
//can: endof:=these values are for the geometric stiffnes similar part and its convergence 


Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];

transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchytens); 
detF=det(Ftens);


/*
Cauchy_rate[0]=1.0*_coeff*(Cauchytens[0][0]-1.0);
Cauchy_rate[1]=1.0*_coeff*(Cauchytens[1][1]-1.0);
Cauchy_rate[2]=1.0*_coeff*(Cauchytens[2][2]-1.0);
Cauchy_rate[3]=1.0*2.0*_coeff*Cauchytens[0][1];
Cauchy_rate[4]=1.0*2.0*_coeff*Cauchytens[1][2];
Cauchy_rate[5]=1.0*2.0*_coeff*Cauchytens[2][0];
*/

Cauchy_rate[0]=1.0*_coeff*(Cauchytens[0][0]-1.0);
Cauchy_rate[1]=1.0*_coeff*(Cauchytens[1][1]-1.0);
Cauchy_rate[2]=1.0*_coeff*(Cauchytens[2][2]-1.0);
Cauchy_rate[3]=1.0*_coeff*Cauchytens[0][1];
Cauchy_rate[4]=1.0*_coeff*Cauchytens[1][2];
Cauchy_rate[5]=1.0*_coeff*Cauchytens[2][0];



//compüuting the pore elasticity
Pore_Elasticity_at_Gauss(D_pore, F, _debugfile,1.0);

for(int i=0;i<6;i++)
for(int j=0;j<6;j++)
D_pore[i][j]=-0.5*D_pore[i][j];



//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N_lambda[shape_counter+0]=0.0;
N_lambda[shape_counter+1]=0.0;
N_lambda[shape_counter+2]=0.0;
N_lambda[shape_counter+3]=0.125*(1+sign_ksi*_ksi)*(1+sign_eta*_eta)*(1+sign_zeta*_zeta);
shape_counter=shape_counter+4;
}
/*
for(int n=0;n<4;n++)
{
double temp=N_lambda[2*4+n];
N_lambda[2*4+n]=N_lambda[3*4+n];
N_lambda[3*4+n]=temp;
temp=N_lambda[6*4+n];
N_lambda[6*4+n]=N_lambda[7*4+n];
N_lambda[7*4+n]=temp;
}
*/
double temp=N_lambda[11];
N_lambda[11]=N_lambda[15];
N_lambda[15]=temp;
temp=N_lambda[27];
N_lambda[27]=N_lambda[31];
N_lambda[31]=temp;
//end definition shape value vector





//computation of the third part 
dNdksi[0]=-0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[1]=-0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[2]=-0.125*(1.0-_ksi)*(1.0-_eta);


dNdksi[3]=0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[4]=-0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[5]=-0.125*(1.0+_ksi)*(1.0-_eta);


dNdksi[6]=0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[7]=0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[8]=-0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[9]=-0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[10]=0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[11]=-0.125*(1.0-_ksi)*(1.0+_eta);

dNdksi[12]=-0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[13]=-0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[14]=0.125*(1.0-_ksi)*(1.0-_eta);

dNdksi[15]=0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[16]=-0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[17]=0.125*(1.0+_ksi)*(1.0-_eta);

dNdksi[18]=0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[19]=0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[20]=0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[21]=-0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[22]=0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[23]=0.125*(1.0-_ksi)*(1.0+_eta);




for(int i=0;i<24;i++)
{
for(int j=0;j<24;j++)
inversejacob[i][j]=0.0; 
dNdX[i]=0.0; 
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0;


for(int i=0;i<3;i++) // row label 
for(int j=0;j<8;j++) 
{
jacob[i][0]=jacob[i][0]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[0];
jacob[i][1]=jacob[i][1]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[1];
jacob[i][2]=jacob[i][2]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[2];
}

detj=jacob[0][0]*(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])-
     jacob[0][1]*(jacob[1][0]*jacob[2][2]-jacob[2][0]*jacob[1][2])+
     jacob[0][2]*(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1]);


inversejacob[0][0]=inversejacob[3][3]=inversejacob[6][6]=inversejacob[9][9]=inversejacob[12][12]=inversejacob[15][15]=inversejacob[18][18]=inversejacob[21][21]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
inversejacob[0][1]=inversejacob[3][4]=inversejacob[6][7]=inversejacob[9][10]=inversejacob[12][13]=inversejacob[15][16]=inversejacob[18][19]=inversejacob[21][22]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
inversejacob[0][2]=inversejacob[3][5]=inversejacob[6][8]=inversejacob[9][11]=inversejacob[12][14]=inversejacob[15][17]=inversejacob[18][20]=inversejacob[21][23]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj; 

inversejacob[1][0]=inversejacob[4][3]=inversejacob[7][6]=inversejacob[10][9]=inversejacob[13][12]=inversejacob[16][15]=inversejacob[19][18]=inversejacob[22][21]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
inversejacob[1][1]=inversejacob[4][4]=inversejacob[7][7]=inversejacob[10][10]=inversejacob[13][13]=inversejacob[16][16]=inversejacob[19][19]=inversejacob[22][22]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
inversejacob[1][2]=inversejacob[4][5]=inversejacob[7][8]=inversejacob[10][11]=inversejacob[13][14]=inversejacob[16][17]=inversejacob[19][20]=inversejacob[22][23]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;

inversejacob[2][0]=inversejacob[5][3]=inversejacob[8][6]=inversejacob[11][9]=inversejacob[14][12]=inversejacob[17][15]=inversejacob[20][18]=inversejacob[23][21]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;
inversejacob[2][1]=inversejacob[5][4]=inversejacob[8][7]=inversejacob[11][10]=inversejacob[14][13]=inversejacob[17][16]=inversejacob[20][19]=inversejacob[23][22]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;
inversejacob[2][2]=inversejacob[5][5]=inversejacob[8][8]=inversejacob[11][11]=inversejacob[14][14]=inversejacob[17][17]=inversejacob[20][20]=inversejacob[23][23]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;


for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
dNdX[i]=dNdX[i]+inversejacob[i][j]*dNdksi[j];


double first_v[32]; //left node driven 
double second_v[6]; //right node driven 
double third_v[32]; //right node driven 

for(int n=0;n<32;n++)
first_v[n]=0.0;


for(int n=0;n<32;n++)
first_v[n]=1.0*detF*N_lambda[n];

for(int n=0;n<6;n++)
second_v[n]=0.0;
for(int n=0;n<32;n++)
third_v[n]=0.0;

for(int i=0;i<6;i++)
for(int j=0;j<6;j++)
second_v[i]=second_v[i]+Cauchy_rate[j]*D_pore[i][j];

for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
third_v[i]=third_v[i]+second_v[j]*B[j][i];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
{
_K_M1uMat[j][i]=_K_M1uMat[j][i]+first_v[i]*third_v[j];
//_K_M1uMat[i][j]=_K_M1uMat[j][i];
}

/*
double Ktrans[32][32];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
Ktrans[j][i]=_K_M1uMat[i][j];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_M1uMat[i][j]=0.5*(_K_M1uMat[i][j]+Ktrans[i][j]);
*/

/*
for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
if(i%4==3)
_K_M1uMat[j][i]=_K_M1uMat[i][j];
*/


//integration over the left and right nodes 
/*

//loop over the right node 
for(int L=0;L<8;L++)
{
//for the sake of easiness, there are in total 3 vectors splitted 
double first_v[32]; //left node driven 
double second_v[6]; //right node driven 
double third_v[32]; //right node driven 

//they are all supposed to be initialized as zero 
for(int n=0;n<32;n++)
first_v[n]=0.0;


for(int n=0;n<4;n++)
first_v[L*4+n]=1.0/detF*N_lambda[L*4+n];
//loop over the right node
for(int M=0;M<8;M++)
{

for(int n=0;n<6;n++)
second_v[n]=0.0;
for(int n=0;n<32;n++)
third_v[n]=0.0;


for(int i=0;i<6;i++)
for(int j=0;j<6;j++)
second_v[i]=second_v[i]+Cauchy_rate[j]*D_pore[i][j];

for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
if(i/4==M) //only the corresponding nodes are taken 
third_v[i]=third_v[i]+second_v[j]*B[j][i];



for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_M1uMat[i][j]=_K_M1uMat[i][j]+first_v[i]*third_v[j];

}
}


*/
};


void calc_K_M1uN_at_gauss(Node _nodes[],
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
			  double _coeff)
{


double beta_2=0.25;
double beta_1=0.5;
double Ftens[3][3];
double FtensT[3][3];
double FtensInv[3][3];
double Cauchytens[3][3];
double CauchytensInv[3][3];
double Cauchy_rate[6];
double CauchyInv[6];
double dNdX[24]; 
double detF;



//can: startof:=these values are for the geometric stiffnes similar part and its convergence 
double dNdksi[24]; 
double inversejacob[24][24]; 
double jacob[3][3]; 
//can: endof:=these values are for the geometric stiffnes similar part and its convergence 


Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];



transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchytens); 
inverse(Cauchytens,CauchytensInv);
detF=det(Ftens);


CauchyInv[0]=CauchytensInv[0][0];
CauchyInv[1]=CauchytensInv[1][1];
CauchyInv[2]=CauchytensInv[2][2];
CauchyInv[3]=CauchytensInv[0][1];
CauchyInv[4]=CauchytensInv[1][2];
CauchyInv[5]=CauchytensInv[2][0];



Cauchy_rate[0]=1.0*_coeff*(Cauchytens[0][0]-1.0);
Cauchy_rate[1]=1.0*_coeff*(Cauchytens[1][1]-1.0);
Cauchy_rate[2]=1.0*_coeff*(Cauchytens[2][2]-1.0);
Cauchy_rate[3]=1.0*_coeff*Cauchytens[0][1];
Cauchy_rate[4]=1.0*_coeff*Cauchytens[1][2];
Cauchy_rate[5]=1.0*_coeff*Cauchytens[2][0];

double somecoeff=0.0;
for(int n=0;n<6;n++)
somecoeff=somecoeff+Cauchy_rate[n]*CauchyInv[n]*detF;






//computation of the third part 
dNdksi[0]=-0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[1]=-0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[2]=-0.125*(1.0-_ksi)*(1.0-_eta);


dNdksi[3]=0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[4]=-0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[5]=-0.125*(1.0+_ksi)*(1.0-_eta);


dNdksi[6]=0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[7]=0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[8]=-0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[9]=-0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[10]=0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[11]=-0.125*(1.0-_ksi)*(1.0+_eta);

dNdksi[12]=-0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[13]=-0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[14]=0.125*(1.0-_ksi)*(1.0-_eta);

dNdksi[15]=0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[16]=-0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[17]=0.125*(1.0+_ksi)*(1.0-_eta);

dNdksi[18]=0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[19]=0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[20]=0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[21]=-0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[22]=0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[23]=0.125*(1.0-_ksi)*(1.0+_eta);




for(int i=0;i<24;i++)
{
for(int j=0;j<24;j++)
inversejacob[i][j]=0.0; 
dNdX[i]=0.0; 
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0;


for(int i=0;i<3;i++) // row label 
for(int j=0;j<8;j++) 
{
jacob[i][0]=jacob[i][0]+dNdksi[j*3+i]*(_nodes[_elements[_currentelement].connid[j]].X_lamb[0]+_nodes[_elements[_currentelement].connid[j]].u_lamb[0]);
jacob[i][1]=jacob[i][1]+dNdksi[j*3+i]*(_nodes[_elements[_currentelement].connid[j]].X_lamb[1]+_nodes[_elements[_currentelement].connid[j]].u_lamb[1]);
jacob[i][2]=jacob[i][2]+dNdksi[j*3+i]*(_nodes[_elements[_currentelement].connid[j]].X_lamb[2]+_nodes[_elements[_currentelement].connid[j]].u_lamb[2]);
}

detj=jacob[0][0]*(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])-
     jacob[0][1]*(jacob[1][0]*jacob[2][2]-jacob[2][0]*jacob[1][2])+
     jacob[0][2]*(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1]);


inversejacob[0][0]=inversejacob[3][3]=inversejacob[6][6]=inversejacob[9][9]=inversejacob[12][12]=inversejacob[15][15]=inversejacob[18][18]=inversejacob[21][21]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
inversejacob[0][1]=inversejacob[3][4]=inversejacob[6][7]=inversejacob[9][10]=inversejacob[12][13]=inversejacob[15][16]=inversejacob[18][19]=inversejacob[21][22]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
inversejacob[0][2]=inversejacob[3][5]=inversejacob[6][8]=inversejacob[9][11]=inversejacob[12][14]=inversejacob[15][17]=inversejacob[18][20]=inversejacob[21][23]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj; 

inversejacob[1][0]=inversejacob[4][3]=inversejacob[7][6]=inversejacob[10][9]=inversejacob[13][12]=inversejacob[16][15]=inversejacob[19][18]=inversejacob[22][21]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
inversejacob[1][1]=inversejacob[4][4]=inversejacob[7][7]=inversejacob[10][10]=inversejacob[13][13]=inversejacob[16][16]=inversejacob[19][19]=inversejacob[22][22]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
inversejacob[1][2]=inversejacob[4][5]=inversejacob[7][8]=inversejacob[10][11]=inversejacob[13][14]=inversejacob[16][17]=inversejacob[19][20]=inversejacob[22][23]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;

inversejacob[2][0]=inversejacob[5][3]=inversejacob[8][6]=inversejacob[11][9]=inversejacob[14][12]=inversejacob[17][15]=inversejacob[20][18]=inversejacob[23][21]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;
inversejacob[2][1]=inversejacob[5][4]=inversejacob[8][7]=inversejacob[11][10]=inversejacob[14][13]=inversejacob[17][16]=inversejacob[20][19]=inversejacob[23][22]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;
inversejacob[2][2]=inversejacob[5][5]=inversejacob[8][8]=inversejacob[11][11]=inversejacob[14][14]=inversejacob[17][17]=inversejacob[20][20]=inversejacob[23][23]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;


for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
dNdX[i]=dNdX[i]+inversejacob[i][j]*dNdksi[j];

/*
inverse(Ftens,FtensInv);

for(int i=0;i<8;i++)
{
double currentdNdX[3];
double currentdNdu[3];
for(int j=0;j<3;j++)
currentdNdu[j]=0.0;

currentdNdX[0]=dNdX[i*3+0];
currentdNdX[1]=dNdX[i*3+1];
currentdNdX[2]=dNdX[i*3+2];

for(int j=0;j<3;j++)
for(int k=0;k<3;k++)
currentdNdu[j]=currentdNdu[j]+FtensInv[j][k]*currentdNdX[k];

dNdX[i*3+0]=currentdNdu[0];
dNdX[i*3+1]=currentdNdu[1];
dNdX[i*3+2]=currentdNdu[2];



}
*/

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
{
if(i%4==3 && j%4==0)
{
_K_M1uN[i][j+0]=_K_M1uN[i][j+0]+somecoeff*dNdX[i-i/3];
_K_M1uN[i][j+1]=_K_M1uN[i][j+1]+somecoeff*dNdX[i-i/3+1];
_K_M1uN[i][j+2]=_K_M1uN[i][j+2]+somecoeff*dNdX[i-i/3+2];
}
}


};




void calc_K_M2uJ_at_gauss(Node _nodes[],
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
			  double _current_layup)
{


double beta_2=0.25;
double beta_1=0.5;
double Ftens[3][3];
double FtensT[3][3];
double Cauchytens[3][3];
double CauchytensInv[3][3];
double CauchyInv[6];
double dNdX[24]; 
double detF;
double N[8];
double beta=0.0;
double alfa;

//can: startof:=these values are for the geometric stiffnes similar part and its convergence 
double dNdksi[24]; 
double inversejacob[24][24]; 
double jacob[3][3]; 
//can: endof:=these values are for the geometric stiffnes similar part and its convergence 
double dlambda_dX[3];
double Lambda_var[6][32];
double current_nS=0.0;
double current_nf;
double current_ns;

double dlambda_dX_f[3];
double dlambda_dX_pf[3];
double M_vec[3];
double M_met[3][3];

M_vec[0]=1.0;
M_vec[1]=0.0;
M_vec[2]=0.0;


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
M_met[i][j]=M_vec[i]*M_vec[j];


Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];

transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchytens); 
inverse(Cauchytens,CauchytensInv);
detF=det(Ftens);



CauchyInv[0]=CauchytensInv[0][0];
CauchyInv[1]=CauchytensInv[1][1];
CauchyInv[2]=CauchytensInv[2][2];
CauchyInv[3]=CauchytensInv[0][1];
CauchyInv[4]=CauchytensInv[1][2];
CauchyInv[5]=CauchytensInv[2][0];





//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N[shape_counter]=0.125*(1+sign_ksi*_ksi)*(1+sign_eta*_eta)*(1+sign_zeta*_zeta);
shape_counter++;
}
double temp=N[2];
N[2]=N[3];
N[3]=temp;
temp=N[6];
N[6]=N[7];
N[7]=temp;
//end definition shape value vector

for(int n=0;n<8;n++)
beta=beta+N[n]*_nodes[_elements[_currentelement].connid[n]].inpermeability;

//impermeability in fiber direction
alfa=beta/100.0;

if(_current_layup==1.0)
{
alfa=1000.0;
beta=1000.0;
}


for(int s=0;s<8;s++)
current_nS=current_nS+N[s]*_nodes[_elements[_currentelement].connid[s]].nS;
current_ns=current_nS/detF;
current_nf=1.0-current_ns;


if(_zeta<0.0 && _eta<0.0 && _ksi<0.0)
_elements[_currentelement].ns[0]=current_ns;
if(_zeta<0.0 && _eta<0.0 && _ksi>0.0)
_elements[_currentelement].ns[1]=current_ns;
if(_zeta<0.0 && _eta>0.0 && _ksi>0.0)
_elements[_currentelement].ns[2]=current_ns;
if(_zeta<0.0 && _eta>0.0 && _ksi<0.0)
_elements[_currentelement].ns[3]=current_ns;
if(_zeta>0.0 && _eta<0.0 && _ksi<0.0)
_elements[_currentelement].ns[4]=current_ns;
if(_zeta>0.0 && _eta<0.0 && _ksi>0.0)
_elements[_currentelement].ns[5]=current_ns;
if(_zeta>0.0 && _eta>0.0 && _ksi>0.0)
_elements[_currentelement].ns[6]=current_ns;
if(_zeta>0.0 && _eta>0.0 && _ksi<0.0)
_elements[_currentelement].ns[7]=current_ns;






double some_deriv=0.0;
some_deriv=current_nf*current_nf+2.0*current_nf*current_nS/detF;
//some_deriv=-2.0*current_nS*(1.0-current_nS/detF)/detF+current_nf*current_nf;

double coeff=0.0;
coeff=(0.5*detF)*some_deriv;
//coeff=(1.0/beta)*(-0.5/detF)*some_deriv;




//computation of the third part 
dNdksi[0]=-0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[1]=-0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[2]=-0.125*(1.0-_ksi)*(1.0-_eta);


dNdksi[3]=0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[4]=-0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[5]=-0.125*(1.0+_ksi)*(1.0-_eta);


dNdksi[6]=0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[7]=0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[8]=-0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[9]=-0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[10]=0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[11]=-0.125*(1.0-_ksi)*(1.0+_eta);

dNdksi[12]=-0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[13]=-0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[14]=0.125*(1.0-_ksi)*(1.0-_eta);

dNdksi[15]=0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[16]=-0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[17]=0.125*(1.0+_ksi)*(1.0-_eta);

dNdksi[18]=0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[19]=0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[20]=0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[21]=-0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[22]=0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[23]=0.125*(1.0-_ksi)*(1.0+_eta);




for(int i=0;i<24;i++)
{
for(int j=0;j<24;j++)
inversejacob[i][j]=0.0; 
dNdX[i]=0.0; 
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0;


for(int i=0;i<3;i++) // row label 
for(int j=0;j<8;j++) 
{
jacob[i][0]=jacob[i][0]+dNdksi[j*3+i]*(_nodes[_elements[_currentelement].connid[j]].X_lamb[0]+_nodes[_elements[_currentelement].connid[j]].u_lamb[0]);
jacob[i][1]=jacob[i][1]+dNdksi[j*3+i]*(_nodes[_elements[_currentelement].connid[j]].X_lamb[1]+_nodes[_elements[_currentelement].connid[j]].u_lamb[1]);
jacob[i][2]=jacob[i][2]+dNdksi[j*3+i]*(_nodes[_elements[_currentelement].connid[j]].X_lamb[2]+_nodes[_elements[_currentelement].connid[j]].u_lamb[2]);
}

detj=jacob[0][0]*(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])-
     jacob[0][1]*(jacob[1][0]*jacob[2][2]-jacob[2][0]*jacob[1][2])+
     jacob[0][2]*(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1]);


inversejacob[0][0]=inversejacob[3][3]=inversejacob[6][6]=inversejacob[9][9]=inversejacob[12][12]=inversejacob[15][15]=inversejacob[18][18]=inversejacob[21][21]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
inversejacob[0][1]=inversejacob[3][4]=inversejacob[6][7]=inversejacob[9][10]=inversejacob[12][13]=inversejacob[15][16]=inversejacob[18][19]=inversejacob[21][22]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
inversejacob[0][2]=inversejacob[3][5]=inversejacob[6][8]=inversejacob[9][11]=inversejacob[12][14]=inversejacob[15][17]=inversejacob[18][20]=inversejacob[21][23]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj; 

inversejacob[1][0]=inversejacob[4][3]=inversejacob[7][6]=inversejacob[10][9]=inversejacob[13][12]=inversejacob[16][15]=inversejacob[19][18]=inversejacob[22][21]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
inversejacob[1][1]=inversejacob[4][4]=inversejacob[7][7]=inversejacob[10][10]=inversejacob[13][13]=inversejacob[16][16]=inversejacob[19][19]=inversejacob[22][22]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
inversejacob[1][2]=inversejacob[4][5]=inversejacob[7][8]=inversejacob[10][11]=inversejacob[13][14]=inversejacob[16][17]=inversejacob[19][20]=inversejacob[22][23]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;

inversejacob[2][0]=inversejacob[5][3]=inversejacob[8][6]=inversejacob[11][9]=inversejacob[14][12]=inversejacob[17][15]=inversejacob[20][18]=inversejacob[23][21]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;
inversejacob[2][1]=inversejacob[5][4]=inversejacob[8][7]=inversejacob[11][10]=inversejacob[14][13]=inversejacob[17][16]=inversejacob[20][19]=inversejacob[23][22]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;
inversejacob[2][2]=inversejacob[5][5]=inversejacob[8][8]=inversejacob[11][11]=inversejacob[14][14]=inversejacob[17][17]=inversejacob[20][20]=inversejacob[23][23]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;


for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
dNdX[i]=dNdX[i]+inversejacob[i][j]*dNdksi[j];




for(int i=0;i<3;i++)
{
dlambda_dX[i]=0.0;
dlambda_dX_f[i]=0.0;
dlambda_dX_pf[i]=0.0;
}





for(int i=0;i<24;i++) 
dlambda_dX[i%3]=dlambda_dX[i%3]+dNdX[i]*_nodes[_elements[_currentelement].connid[i/3]].u_lamb[3];

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
dlambda_dX_f[i]=dlambda_dX_f[i]+M_met[i][j]*dlambda_dX[j];

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
dlambda_dX_pf[i]=dlambda_dX_pf[i]+(cronecker(i,j)-M_met[i][j])*dlambda_dX[j];


double seepage[3];
for(int i=0;i<3;i++)
seepage[i]=-1.0/beta*current_nf*dlambda_dX_pf[i]-1.0/alfa*current_nf*dlambda_dX_f[i];

double v_solid[3];
for(int i=0;i<3;i++)
v_solid[i]=0.0;

for(int s=0;s<8;s++)
for(int i=0;i<3;i++)
v_solid[i]=v_solid[i]+_coeff*N[s]*_nodes[_elements[_currentelement].connid[s]].u_lamb[i];


if(_zeta<0.0 && _eta<0.0 && _ksi<0.0)
for(int i=0;i<3;i++)
{
_elements[_currentelement].v_fluid[0*3+i]=seepage[i]+v_solid[i];
_elements[_currentelement].v_fs[0*3+i]=seepage[i];
}
if(_zeta<0.0 && _eta<0.0 && _ksi>0.0)
for(int i=0;i<3;i++)
{
_elements[_currentelement].v_fluid[1*3+i]=seepage[i]+v_solid[i];
_elements[_currentelement].v_fs[1*3+i]=seepage[i];
}
if(_zeta<0.0 && _eta>0.0 && _ksi>0.0)
for(int i=0;i<3;i++)
{
_elements[_currentelement].v_fluid[2*3+i]=seepage[i]+v_solid[i];
_elements[_currentelement].v_fs[2*3+i]=seepage[i];
}
if(_zeta<0.0 && _eta>0.0 && _ksi<0.0)
for(int i=0;i<3;i++)
{
_elements[_currentelement].v_fluid[3*3+i]=seepage[i]+v_solid[i];
_elements[_currentelement].v_fs[3*3+i]=seepage[i];
}
if(_zeta>0.0 && _eta<0.0 && _ksi<0.0)
for(int i=0;i<3;i++)
{
_elements[_currentelement].v_fluid[4*3+i]=seepage[i]+v_solid[i];
_elements[_currentelement].v_fs[4*3+i]=seepage[i];
}
if(_zeta>0.0 && _eta<0.0 && _ksi>0.0)
for(int i=0;i<3;i++)
{
_elements[_currentelement].v_fluid[5*3+i]=seepage[i]+v_solid[i];
_elements[_currentelement].v_fs[5*3+i]=seepage[i];
}
if(_zeta>0.0 && _eta>0.0 && _ksi>0.0)
for(int i=0;i<3;i++)
{
_elements[_currentelement].v_fluid[6*3+i]=seepage[i]+v_solid[i];
_elements[_currentelement].v_fs[6*3+i]=seepage[i];
}
if(_zeta>0.0 && _eta>0.0 && _ksi<0.0)
for(int i=0;i<3;i++)
{
_elements[_currentelement].v_fluid[7*3+i]=seepage[i]+v_solid[i];
_elements[_currentelement].v_fs[7*3+i]=seepage[i];
}














for(int i=0;i<6;i++)
for(int j=0;j<32;j++)
Lambda_var[i][j]=0.0;


for(int i=0;i<8;i++)
{
//the first three columns are set to be zero 
Lambda_var[0][i*4+3]=(1.0*dNdX[i*3+0]*dlambda_dX_f[0])/alfa;
Lambda_var[1][i*4+3]=(1.0*dNdX[i*3+1]*dlambda_dX_f[1])/alfa;
Lambda_var[2][i*4+3]=(1.0*dNdX[i*3+2]*dlambda_dX_f[2])/alfa;
Lambda_var[3][i*4+3]=(0.5*dNdX[i*3+0]*dlambda_dX_f[1]+0.5*dNdX[i*3+1]*dlambda_dX_f[0])/alfa;
Lambda_var[4][i*4+3]=(0.5*dNdX[i*3+1]*dlambda_dX_f[2]+0.5*dNdX[i*3+2]*dlambda_dX_f[1])/alfa;
Lambda_var[5][i*4+3]=(0.5*dNdX[i*3+2]*dlambda_dX_f[0]+0.5*dNdX[i*3+0]*dlambda_dX_f[2])/alfa;

Lambda_var[0][i*4+3]=Lambda_var[0][i*4+3]+(1.0*dNdX[i*3+0]*dlambda_dX_pf[0])/beta;
Lambda_var[1][i*4+3]=Lambda_var[1][i*4+3]+(1.0*dNdX[i*3+1]*dlambda_dX_pf[1])/beta;
Lambda_var[2][i*4+3]=Lambda_var[2][i*4+3]+(1.0*dNdX[i*3+2]*dlambda_dX_pf[2])/beta;
Lambda_var[3][i*4+3]=Lambda_var[3][i*4+3]+(0.5*dNdX[i*3+0]*dlambda_dX_pf[1]+0.5*dNdX[i*3+1]*dlambda_dX_pf[0])/beta;
Lambda_var[4][i*4+3]=Lambda_var[4][i*4+3]+(0.5*dNdX[i*3+1]*dlambda_dX_pf[2]+0.5*dNdX[i*3+2]*dlambda_dX_pf[1])/beta;
Lambda_var[5][i*4+3]=Lambda_var[5][i*4+3]+(0.5*dNdX[i*3+2]*dlambda_dX_pf[0]+0.5*dNdX[i*3+0]*dlambda_dX_pf[2])/beta;
} 






double first_v[32]; //left node driven 
double second_v[32]; //right node driven 

for(int n=0;n<32;n++)
first_v[n]=0.0;

for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
first_v[i]=first_v[i]+coeff*CauchyInv[j]*Lambda_var[j][i];

for(int n=0;n<32;n++)
second_v[n]=0.0;


for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
second_v[i]=second_v[i]+CauchyInv[j]*B[j][i];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
{
_K_M2uJ[j][i]=_K_M2uJ[j][i]+first_v[i]*second_v[j];
//_K_M2uJ[i][j]=_K_M2uJ[j][i];
}
/*
double Ktrans[32][32];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
Ktrans[j][i]=_K_M2uJ[i][j];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_M2uJ[i][j]=0.5*(_K_M2uJ[i][j]+Ktrans[i][j]);
*/
//HERE-> VOLUME FRACTIONS 

/*
for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
if(i%4==3)
_K_M2uJ[j][i]=_K_M2uJ[i][j];
*/


//integration over the left and right nodes 

/*
//loop over the right node 
for(int L=0;L<8;L++)
{
//for the sake of easiness, there are in total 3 vectors splitted 
double first_v[32]; //left node driven 
double second_v[32]; //right node driven 

//they are all supposed to be initialized as zero 
for(int n=0;n<32;n++)
first_v[n]=0.0;



for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
if(i/4==L) //only the corresponding nodes are taken 
first_v[i]=first_v[i]+coeff*CauchyInv[j]*Lambda_var[j][i];

//loop over the right node
for(int M=0;M<8;M++)
{

for(int n=0;n<32;n++)
second_v[n]=0.0;


for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
if(i/4==M) //only the corresponding nodes are taken 
second_v[i]=second_v[i]+CauchyInv[j]*B[j][i];



for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_M2uJ[i][j]=_K_M2uJ[i][j]+first_v[i]*second_v[j];

}
}

*/

};



//to be corrected 



void calc_K_M2uMat_at_gauss(Node _nodes[],
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
			    double _current_layup)
{


double beta_2=0.25;
double beta_1=0.5;
double Ftens[3][3];
double FtensT[3][3];
double dNdX[24]; 
double detF;
double N[8];
double beta=0.0;
double alfa;
double D_pore[6][6];

//can: startof:=these values are for the geometric stiffnes similar part and its convergence 
double dNdksi[24]; 
double inversejacob[24][24]; 
double jacob[3][3]; 
//can: endof:=these values are for the geometric stiffnes similar part and its convergence 
double dlambda_dX[3];
double Lambda_var[6][32];
double current_nS=0.0;
double current_nf;

double dlambda_dX_f[3];
double dlambda_dX_pf[3];
double M_vec[3];
double M_met[3][3];


M_vec[0]=1.0;
M_vec[1]=0.0;
M_vec[2]=0.0;


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
M_met[i][j]=M_vec[i]*M_vec[j];

Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];

detF=det(Ftens);


//compüuting the pore elasticity
Pore_Elasticity_at_Gauss(D_pore, F, _debugfile,1.0);

for(int i=0;i<6;i++)
for(int j=0;j<6;j++)
D_pore[i][j]=-0.5*D_pore[i][j];




//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N[shape_counter]=0.125*(1+sign_ksi*_ksi)*(1+sign_eta*_eta)*(1+sign_zeta*_zeta);
shape_counter++;
}
double temp=N[2];
N[2]=N[3];
N[3]=temp;
temp=N[6];
N[6]=N[7];
N[7]=temp;
//end definition shape value vector

for(int n=0;n<8;n++)
beta=beta+N[n]*_nodes[_elements[_currentelement].connid[n]].inpermeability;
//impermeability in fiber direction
alfa=beta/100.0;
/*
if(_current_layup==1.0)
{
alfa=1000.0;
beta=1000.0;
}
*/

for(int s=0;s<8;s++)
current_nS=current_nS+N[s]*_nodes[_elements[_currentelement].connid[s]].nS;
current_nf=1.0-current_nS/detF;





//computation of the third part 
dNdksi[0]=-0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[1]=-0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[2]=-0.125*(1.0-_ksi)*(1.0-_eta);


dNdksi[3]=0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[4]=-0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[5]=-0.125*(1.0+_ksi)*(1.0-_eta);


dNdksi[6]=0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[7]=0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[8]=-0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[9]=-0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[10]=0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[11]=-0.125*(1.0-_ksi)*(1.0+_eta);

dNdksi[12]=-0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[13]=-0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[14]=0.125*(1.0-_ksi)*(1.0-_eta);

dNdksi[15]=0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[16]=-0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[17]=0.125*(1.0+_ksi)*(1.0-_eta);

dNdksi[18]=0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[19]=0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[20]=0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[21]=-0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[22]=0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[23]=0.125*(1.0-_ksi)*(1.0+_eta);




for(int i=0;i<24;i++)
{
for(int j=0;j<24;j++)
inversejacob[i][j]=0.0; 
dNdX[i]=0.0; 
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0;


for(int i=0;i<3;i++) // row label 
for(int j=0;j<8;j++) 
{
jacob[i][0]=jacob[i][0]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[0];
jacob[i][1]=jacob[i][1]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[1];
jacob[i][2]=jacob[i][2]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[2];
}

detj=jacob[0][0]*(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])-
     jacob[0][1]*(jacob[1][0]*jacob[2][2]-jacob[2][0]*jacob[1][2])+
     jacob[0][2]*(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1]);


inversejacob[0][0]=inversejacob[3][3]=inversejacob[6][6]=inversejacob[9][9]=inversejacob[12][12]=inversejacob[15][15]=inversejacob[18][18]=inversejacob[21][21]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
inversejacob[0][1]=inversejacob[3][4]=inversejacob[6][7]=inversejacob[9][10]=inversejacob[12][13]=inversejacob[15][16]=inversejacob[18][19]=inversejacob[21][22]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
inversejacob[0][2]=inversejacob[3][5]=inversejacob[6][8]=inversejacob[9][11]=inversejacob[12][14]=inversejacob[15][17]=inversejacob[18][20]=inversejacob[21][23]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj; 

inversejacob[1][0]=inversejacob[4][3]=inversejacob[7][6]=inversejacob[10][9]=inversejacob[13][12]=inversejacob[16][15]=inversejacob[19][18]=inversejacob[22][21]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
inversejacob[1][1]=inversejacob[4][4]=inversejacob[7][7]=inversejacob[10][10]=inversejacob[13][13]=inversejacob[16][16]=inversejacob[19][19]=inversejacob[22][22]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
inversejacob[1][2]=inversejacob[4][5]=inversejacob[7][8]=inversejacob[10][11]=inversejacob[13][14]=inversejacob[16][17]=inversejacob[19][20]=inversejacob[22][23]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;

inversejacob[2][0]=inversejacob[5][3]=inversejacob[8][6]=inversejacob[11][9]=inversejacob[14][12]=inversejacob[17][15]=inversejacob[20][18]=inversejacob[23][21]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;
inversejacob[2][1]=inversejacob[5][4]=inversejacob[8][7]=inversejacob[11][10]=inversejacob[14][13]=inversejacob[17][16]=inversejacob[20][19]=inversejacob[23][22]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;
inversejacob[2][2]=inversejacob[5][5]=inversejacob[8][8]=inversejacob[11][11]=inversejacob[14][14]=inversejacob[17][17]=inversejacob[20][20]=inversejacob[23][23]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;


for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
dNdX[i]=dNdX[i]+inversejacob[i][j]*dNdksi[j];




for(int i=0;i<3;i++)
{
dlambda_dX_f[i]=0.0;
dlambda_dX[i]=0.0;
dlambda_dX_pf[i]=0.0;
}



for(int i=0;i<24;i++) 
dlambda_dX[i%3]=dlambda_dX[i%3]+dNdX[i]*_nodes[_elements[_currentelement].connid[i/3]].u_lamb[3];

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
dlambda_dX_f[i]=dlambda_dX_f[i]+M_met[i][j]*dlambda_dX[j];

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
dlambda_dX_pf[i]=dlambda_dX_pf[i]+(cronecker(i,j)-M_met[i][j])*dlambda_dX[j];


for(int i=0;i<6;i++)
for(int j=0;j<32;j++)
Lambda_var[i][j]=0.0;


for(int i=0;i<8;i++)
{
//the first three columns are set to be zero 
Lambda_var[0][i*4+3]=(1.0*dNdX[i*3+0]*dlambda_dX_f[0])/alfa;
Lambda_var[1][i*4+3]=(1.0*dNdX[i*3+1]*dlambda_dX_f[1])/alfa;
Lambda_var[2][i*4+3]=(1.0*dNdX[i*3+2]*dlambda_dX_f[2])/alfa;
Lambda_var[3][i*4+3]=(0.5*dNdX[i*3+0]*dlambda_dX_f[1]+0.5*dNdX[i*3+1]*dlambda_dX_f[0])/alfa;
Lambda_var[4][i*4+3]=(0.5*dNdX[i*3+1]*dlambda_dX_f[2]+0.5*dNdX[i*3+2]*dlambda_dX_f[1])/alfa;
Lambda_var[5][i*4+3]=(0.5*dNdX[i*3+2]*dlambda_dX_f[0]+0.5*dNdX[i*3+0]*dlambda_dX_f[2])/alfa;

Lambda_var[0][i*4+3]=Lambda_var[0][i*4+3]+(1.0*dNdX[i*3+0]*dlambda_dX_pf[0])/beta;
Lambda_var[1][i*4+3]=Lambda_var[1][i*4+3]+(1.0*dNdX[i*3+1]*dlambda_dX_pf[1])/beta;
Lambda_var[2][i*4+3]=Lambda_var[2][i*4+3]+(1.0*dNdX[i*3+2]*dlambda_dX_pf[2])/beta;
Lambda_var[3][i*4+3]=Lambda_var[3][i*4+3]+(0.5*dNdX[i*3+0]*dlambda_dX_pf[1]+0.5*dNdX[i*3+1]*dlambda_dX_pf[0])/beta;
Lambda_var[4][i*4+3]=Lambda_var[4][i*4+3]+(0.5*dNdX[i*3+1]*dlambda_dX_pf[2]+0.5*dNdX[i*3+2]*dlambda_dX_pf[1])/beta;
Lambda_var[5][i*4+3]=Lambda_var[5][i*4+3]+(0.5*dNdX[i*3+2]*dlambda_dX_pf[0]+0.5*dNdX[i*3+0]*dlambda_dX_pf[2])/beta;
} 


double first_v[32][6]; //left node driven 
double second_v[6][32]; //right node driven 

//they are all supposed to be initialized as zero 
for(int n=0;n<32;n++)
for(int j=0;j<6;j++)
first_v[n][j]=0.0;



for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
first_v[i][j]=1.0*current_nf*current_nf*detF*Lambda_var[j][i];



for(int n=0;n<32;n++)
for(int j=0;j<6;j++)
second_v[j][n]=0.0;


for(int j=0;j<6;j++)
for(int k=0;k<6;k++)
for(int i=0;i<32;i++)
second_v[j][i]=second_v[j][i]+D_pore[j][k]*B[k][i];


for(int i=0;i<32;i++)
for(int k=0;k<6;k++)
for(int j=0;j<32;j++)
{
_K_M2uMat[j][i]=_K_M2uMat[j][i]+first_v[i][k]*second_v[k][j];
//_K_M2uMat[i][j]=_K_M2uMat[j][i];
}
/*
double Ktrans[32][32];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
Ktrans[j][i]=_K_M2uMat[i][j];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_M2uMat[i][j]=0.5*(_K_M2uMat[i][j]+Ktrans[i][j]);
*/

/*
for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
if(i%4==3)
_K_M2uMat[j][i]=_K_M2uMat[i][j];
*/

//HERE-> VOLUME FRACTIONS 


//integration over the left and right nodes 

/*
//loop over the right node 
for(int L=0;L<8;L++)
{
//for the sake of easiness, there are in total 3 vectors splitted 
double first_v[32][6]; //left node driven 
double second_v[6][32]; //right node driven 

//they are all supposed to be initialized as zero 
for(int n=0;n<32;n++)
for(int j=0;j<6;j++)
first_v[n][j]=0.0;





for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
if(i/4==L) //only the corresponding nodes are taken 
first_v[i][j]=1.0/beta*current_nf*current_nf/detF*Lambda_var[j][i];

//loop over the right node
for(int M=0;M<8;M++)
{

for(int n=0;n<32;n++)
for(int j=0;j<6;j++)
second_v[j][n]=0.0;


for(int j=0;j<6;j++)
for(int k=0;k<6;k++)
for(int i=0;i<32;i++)
if(i/4==M) //only the corresponding nodes are taken 
second_v[j][i]=second_v[j][i]+D_pore[j][k]*B[k][i];



for(int i=0;i<32;i++)
for(int k=0;k<6;k++)
for(int j=0;j<32;j++)
_K_M2uMat[i][j]=_K_M2uMat[i][j]+first_v[i][k]*second_v[k][j];

}
}

*/

};


void calc_K_M2lambda_at_gauss(Node _nodes[],
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
			      double _current_layup)
{
double beta_2=0.25;
double beta_1=0.005;
double dNdX[24]; 
double detF;
double N[8];
double beta=0.0;
double alfa;
double Ftens[3][3];
double FtensT[3][3];
double Cauchytens[3][3];
double CauchytensInv[3][3];
double CauchyInv[6];

//can: startof:=these values are for the geometric stiffnes similar part and its convergence 
double dNdksi[24]; 
double inversejacob[24][24]; 
double jacob[3][3]; 
//can: endof:=these values are for the geometric stiffnes similar part and its convergence 
double dlambda_dX[3];
double Lambda_var[6][32][32];
double current_nS=0.0;
double current_nf;

double dlambda_dX_f[3];
double dlambda_dX_pf[3];
double M_vec[3];
double M_met[3][3];
M_vec[0]=1.0;
M_vec[1]=0.0;
M_vec[2]=0.0;


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
M_met[i][j]=M_vec[i]*M_vec[j];






Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];
transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchytens); 
inverse(Cauchytens,CauchytensInv);
detF=det(Ftens);






CauchyInv[0]=CauchytensInv[0][0];
CauchyInv[1]=CauchytensInv[1][1];
CauchyInv[2]=CauchytensInv[2][2];
CauchyInv[3]=CauchytensInv[0][1];
CauchyInv[4]=CauchytensInv[1][2];
CauchyInv[5]=CauchytensInv[2][0];











//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N[shape_counter]=0.125*(1+sign_ksi*_ksi)*(1+sign_eta*_eta)*(1+sign_zeta*_zeta);
shape_counter++;
}
double temp=N[2];
N[2]=N[3];
N[3]=temp;
temp=N[6];
N[6]=N[7];
N[7]=temp;
//end definition shape value vector

for(int n=0;n<8;n++)
beta=beta+N[n]*_nodes[_elements[_currentelement].connid[n]].inpermeability;
//impermeability in fiber direction
alfa=beta/100.0;
/*
if(_current_layup==1.0)
{
alfa=1000.0;
beta=1000.0;
}
*/


for(int s=0;s<8;s++)
current_nS=current_nS+N[s]*_nodes[_elements[_currentelement].connid[s]].nS;
current_nf=1.0-current_nS/detF;





//computation of the third part 
dNdksi[0]=-0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[1]=-0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[2]=-0.125*(1.0-_ksi)*(1.0-_eta);


dNdksi[3]=0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[4]=-0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[5]=-0.125*(1.0+_ksi)*(1.0-_eta);


dNdksi[6]=0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[7]=0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[8]=-0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[9]=-0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[10]=0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[11]=-0.125*(1.0-_ksi)*(1.0+_eta);

dNdksi[12]=-0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[13]=-0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[14]=0.125*(1.0-_ksi)*(1.0-_eta);

dNdksi[15]=0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[16]=-0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[17]=0.125*(1.0+_ksi)*(1.0-_eta);

dNdksi[18]=0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[19]=0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[20]=0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[21]=-0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[22]=0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[23]=0.125*(1.0-_ksi)*(1.0+_eta);




for(int i=0;i<24;i++)
{
for(int j=0;j<24;j++)
inversejacob[i][j]=0.0; 
dNdX[i]=0.0; 
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0;


for(int i=0;i<3;i++) // row label 
for(int j=0;j<8;j++) 
{
jacob[i][0]=jacob[i][0]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[0];
jacob[i][1]=jacob[i][1]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[1];
jacob[i][2]=jacob[i][2]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[2];
}

detj=jacob[0][0]*(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])-
     jacob[0][1]*(jacob[1][0]*jacob[2][2]-jacob[2][0]*jacob[1][2])+
     jacob[0][2]*(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1]);


inversejacob[0][0]=inversejacob[3][3]=inversejacob[6][6]=inversejacob[9][9]=inversejacob[12][12]=inversejacob[15][15]=inversejacob[18][18]=inversejacob[21][21]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
inversejacob[0][1]=inversejacob[3][4]=inversejacob[6][7]=inversejacob[9][10]=inversejacob[12][13]=inversejacob[15][16]=inversejacob[18][19]=inversejacob[21][22]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
inversejacob[0][2]=inversejacob[3][5]=inversejacob[6][8]=inversejacob[9][11]=inversejacob[12][14]=inversejacob[15][17]=inversejacob[18][20]=inversejacob[21][23]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj; 

inversejacob[1][0]=inversejacob[4][3]=inversejacob[7][6]=inversejacob[10][9]=inversejacob[13][12]=inversejacob[16][15]=inversejacob[19][18]=inversejacob[22][21]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
inversejacob[1][1]=inversejacob[4][4]=inversejacob[7][7]=inversejacob[10][10]=inversejacob[13][13]=inversejacob[16][16]=inversejacob[19][19]=inversejacob[22][22]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
inversejacob[1][2]=inversejacob[4][5]=inversejacob[7][8]=inversejacob[10][11]=inversejacob[13][14]=inversejacob[16][17]=inversejacob[19][20]=inversejacob[22][23]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;

inversejacob[2][0]=inversejacob[5][3]=inversejacob[8][6]=inversejacob[11][9]=inversejacob[14][12]=inversejacob[17][15]=inversejacob[20][18]=inversejacob[23][21]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;
inversejacob[2][1]=inversejacob[5][4]=inversejacob[8][7]=inversejacob[11][10]=inversejacob[14][13]=inversejacob[17][16]=inversejacob[20][19]=inversejacob[23][22]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;
inversejacob[2][2]=inversejacob[5][5]=inversejacob[8][8]=inversejacob[11][11]=inversejacob[14][14]=inversejacob[17][17]=inversejacob[20][20]=inversejacob[23][23]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;





for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
dNdX[i]=dNdX[i]+inversejacob[i][j]*dNdksi[j];


double dlambda_dX_fi[3];
double dlambda_dX_pfi[3];
for(int i=0;i<3;i++)
{
dlambda_dX_f[i]=0.0;
dlambda_dX[i]=0.0;
dlambda_dX_pf[i]=0.0;
dlambda_dX_fi[i]=0.0;
dlambda_dX_pfi[i]=0.0;
}



for(int i=0;i<24;i++) 
dlambda_dX[i%3]=dlambda_dX[i%3]+dNdX[i]*_nodes[_elements[_currentelement].connid[i/3]].u_lamb[3];



for(int i=0;i<6;i++)
for(int j=0;j<32;j++)
for(int k=0;k<32;k++)
Lambda_var[i][j][k]=0.0;



for(int i=0;i<8;i++)
for(int j=0;j<8;j++)
{
//the first three columns are set to be zero 
//Lambda_var[0][i*4+3][j*4+3]=1.0*dNdX[i*3+0]*dNdX[j*3+0];
//Lambda_var[1][i*4+3][j*4+3]=1.0*dNdX[i*3+1]*dNdX[j*3+1];
//Lambda_var[2][i*4+3][j*4+3]=1.0*dNdX[i*3+2]*dNdX[j*3+2];
//Lambda_var[3][i*4+3][j*4+3]=0.5*dNdX[i*3+0]*dNdX[j*3+1]+0.5*dNdX[i*3+1]*dNdX[j*3+0];
//Lambda_var[4][i*4+3][j*4+3]=0.5*dNdX[i*3+1]*dNdX[j*3+2]+0.5*dNdX[i*3+2]*dNdX[j*3+1];
//Lambda_var[5][i*4+3][j*4+3]=0.5*dNdX[i*3+2]*dNdX[j*3+0]+0.5*dNdX[i*3+0]*dNdX[j*3+2];


double current_convert[3];
current_convert[0]=dNdX[j*3+0];
current_convert[1]=dNdX[j*3+1];
current_convert[2]=dNdX[j*3+2];

double current_convert_i[3];
current_convert_i[0]=dNdX[i*3+0];
current_convert_i[1]=dNdX[i*3+1];
current_convert_i[2]=dNdX[i*3+2];



for(int k=0;k<3;k++)
for(int l=0;l<3;l++)
dlambda_dX_f[k]=dlambda_dX_f[k]+M_met[k][l]*current_convert[l];

for(int k=0;k<3;k++)
for(int l=0;l<3;l++)
dlambda_dX_pf[k]=dlambda_dX_pf[k]+(cronecker(k,l)-M_met[k][l])*current_convert[l];


for(int k=0;k<3;k++)
for(int l=0;l<3;l++)
dlambda_dX_fi[k]=dlambda_dX_fi[k]+M_met[k][l]*current_convert_i[l];

for(int k=0;k<3;k++)
for(int l=0;l<3;l++)
dlambda_dX_pfi[k]=dlambda_dX_pfi[k]+(cronecker(k,l)-M_met[k][l])*current_convert_i[l];


/*
Lambda_var[0][i*4+3][j*4+3]=(1.0*dNdX[i*3+0]*dlambda_dX_f[0])/alfa;
Lambda_var[1][i*4+3][j*4+3]=(1.0*dNdX[i*3+1]*dlambda_dX_f[1])/alfa;
Lambda_var[2][i*4+3][j*4+3]=(1.0*dNdX[i*3+2]*dlambda_dX_f[2])/alfa;
Lambda_var[3][i*4+3][j*4+3]=(0.5*dNdX[i*3+0]*dlambda_dX_f[1]+0.5*dNdX[i*3+1]*dlambda_dX_f[0])/alfa;
Lambda_var[4][i*4+3][j*4+3]=(0.5*dNdX[i*3+1]*dlambda_dX_f[2]+0.5*dNdX[i*3+2]*dlambda_dX_f[1])/alfa;
Lambda_var[5][i*4+3][j*4+3]=(0.5*dNdX[i*3+2]*dlambda_dX_f[0]+0.5*dNdX[i*3+0]*dlambda_dX_f[2])/alfa;

Lambda_var[0][i*4+3][j*4+3]=Lambda_var[0][i*4+3][j*4+3]+(1.0*dNdX[i*3+0]*dlambda_dX_pf[0])/beta;
Lambda_var[1][i*4+3][j*4+3]=Lambda_var[1][i*4+3][j*4+3]+(1.0*dNdX[i*3+1]*dlambda_dX_pf[1])/beta;
Lambda_var[2][i*4+3][j*4+3]=Lambda_var[2][i*4+3][j*4+3]+(1.0*dNdX[i*3+2]*dlambda_dX_pf[2])/beta;
Lambda_var[3][i*4+3][j*4+3]=Lambda_var[3][i*4+3][j*4+3]+(0.5*dNdX[i*3+0]*dlambda_dX_pf[1]+0.5*dNdX[i*3+1]*dlambda_dX_pf[0])/beta;
Lambda_var[4][i*4+3][j*4+3]=Lambda_var[4][i*4+3][j*4+3]+(0.5*dNdX[i*3+1]*dlambda_dX_pf[2]+0.5*dNdX[i*3+2]*dlambda_dX_pf[1])/beta;
Lambda_var[5][i*4+3][j*4+3]=Lambda_var[5][i*4+3][j*4+3]+(0.5*dNdX[i*3+2]*dlambda_dX_pf[0]+0.5*dNdX[i*3+0]*dlambda_dX_pf[2])/beta;
*/



Lambda_var[0][i*4+3][j*4+3]=Lambda_var[0][i*4+3][j*4+3]+0.5*(1.0*dNdX[i*3+0]*dlambda_dX_f[0])/alfa;
Lambda_var[1][i*4+3][j*4+3]=Lambda_var[0][i*4+3][j*4+3]+0.5*(1.0*dNdX[i*3+1]*dlambda_dX_f[1])/alfa;
Lambda_var[2][i*4+3][j*4+3]=Lambda_var[0][i*4+3][j*4+3]+0.5*(1.0*dNdX[i*3+2]*dlambda_dX_f[2])/alfa;
Lambda_var[3][i*4+3][j*4+3]=Lambda_var[0][i*4+3][j*4+3]+0.5*(0.5*dNdX[i*3+0]*dlambda_dX_f[1]+0.5*dNdX[i*3+1]*dlambda_dX_f[0])/alfa;
Lambda_var[4][i*4+3][j*4+3]=Lambda_var[0][i*4+3][j*4+3]+0.5*(0.5*dNdX[i*3+1]*dlambda_dX_f[2]+0.5*dNdX[i*3+2]*dlambda_dX_f[1])/alfa;
Lambda_var[5][i*4+3][j*4+3]=Lambda_var[0][i*4+3][j*4+3]+0.5*(0.5*dNdX[i*3+2]*dlambda_dX_f[0]+0.5*dNdX[i*3+0]*dlambda_dX_f[2])/alfa;

Lambda_var[0][i*4+3][j*4+3]=Lambda_var[0][i*4+3][j*4+3]+0.5*(1.0*dNdX[i*3+0]*dlambda_dX_pf[0])/beta;
Lambda_var[1][i*4+3][j*4+3]=Lambda_var[1][i*4+3][j*4+3]+0.5*(1.0*dNdX[i*3+1]*dlambda_dX_pf[1])/beta;
Lambda_var[2][i*4+3][j*4+3]=Lambda_var[2][i*4+3][j*4+3]+0.5*(1.0*dNdX[i*3+2]*dlambda_dX_pf[2])/beta;
Lambda_var[3][i*4+3][j*4+3]=Lambda_var[3][i*4+3][j*4+3]+0.5*(0.5*dNdX[i*3+0]*dlambda_dX_pf[1]+0.5*dNdX[i*3+1]*dlambda_dX_pf[0])/beta;
Lambda_var[4][i*4+3][j*4+3]=Lambda_var[4][i*4+3][j*4+3]+0.5*(0.5*dNdX[i*3+1]*dlambda_dX_pf[2]+0.5*dNdX[i*3+2]*dlambda_dX_pf[1])/beta;
Lambda_var[5][i*4+3][j*4+3]=Lambda_var[5][i*4+3][j*4+3]+0.5*(0.5*dNdX[i*3+2]*dlambda_dX_pf[0]+0.5*dNdX[i*3+0]*dlambda_dX_pf[2])/beta;


Lambda_var[0][i*4+3][j*4+3]=Lambda_var[0][i*4+3][j*4+3]+0.5*(1.0*dNdX[j*3+0]*dlambda_dX_fi[0])/alfa;
Lambda_var[1][i*4+3][j*4+3]=Lambda_var[1][i*4+3][j*4+3]+0.5*(1.0*dNdX[j*3+1]*dlambda_dX_fi[1])/alfa;
Lambda_var[2][i*4+3][j*4+3]=Lambda_var[2][i*4+3][j*4+3]+0.5*(1.0*dNdX[j*3+2]*dlambda_dX_fi[2])/alfa;
Lambda_var[3][i*4+3][j*4+3]=Lambda_var[3][i*4+3][j*4+3]+0.5*(0.5*dNdX[j*3+0]*dlambda_dX_fi[1]+0.5*dNdX[j*3+1]*dlambda_dX_fi[0])/alfa;
Lambda_var[4][i*4+3][j*4+3]=Lambda_var[4][i*4+3][j*4+3]+0.5*(0.5*dNdX[j*3+1]*dlambda_dX_fi[2]+0.5*dNdX[j*3+2]*dlambda_dX_fi[1])/alfa;
Lambda_var[5][i*4+3][j*4+3]=Lambda_var[5][i*4+3][j*4+3]+0.5*(0.5*dNdX[j*3+2]*dlambda_dX_fi[0]+0.5*dNdX[j*3+0]*dlambda_dX_fi[2])/alfa;

Lambda_var[0][i*4+3][j*4+3]=Lambda_var[0][i*4+3][j*4+3]+0.5*(1.0*dNdX[j*3+0]*dlambda_dX_pfi[0])/beta;
Lambda_var[1][i*4+3][j*4+3]=Lambda_var[1][i*4+3][j*4+3]+0.5*(1.0*dNdX[j*3+1]*dlambda_dX_pfi[1])/beta;
Lambda_var[2][i*4+3][j*4+3]=Lambda_var[2][i*4+3][j*4+3]+0.5*(1.0*dNdX[j*3+2]*dlambda_dX_pfi[2])/beta;
Lambda_var[3][i*4+3][j*4+3]=Lambda_var[3][i*4+3][j*4+3]+0.5*(0.5*dNdX[j*3+0]*dlambda_dX_pfi[1]+0.5*dNdX[j*3+1]*dlambda_dX_pfi[0])/beta;
Lambda_var[4][i*4+3][j*4+3]=Lambda_var[4][i*4+3][j*4+3]+0.5*(0.5*dNdX[j*3+1]*dlambda_dX_pfi[2]+0.5*dNdX[j*3+2]*dlambda_dX_pfi[1])/beta;
Lambda_var[5][i*4+3][j*4+3]=Lambda_var[5][i*4+3][j*4+3]+0.5*(0.5*dNdX[j*3+2]*dlambda_dX_pfi[0]+0.5*dNdX[j*3+0]*dlambda_dX_pfi[2])/beta;




/*
Lambda_var[0][i*4+3][j*4+3]=(1.0*dlambda_dX_fi[0]*dlambda_dX_f[0])/alfa;
Lambda_var[1][i*4+3][j*4+3]=(1.0*dlambda_dX_fi[1]*dlambda_dX_f[1])/alfa;
Lambda_var[2][i*4+3][j*4+3]=(1.0*dlambda_dX_fi[2]*dlambda_dX_f[2])/alfa;
Lambda_var[3][i*4+3][j*4+3]=(0.5*dlambda_dX_fi[0]*dlambda_dX_f[1]+0.5*dlambda_dX_fi[1]*dlambda_dX_f[0])/alfa;
Lambda_var[4][i*4+3][j*4+3]=(0.5*dlambda_dX_fi[1]*dlambda_dX_f[2]+0.5*dlambda_dX_fi[2]*dlambda_dX_f[1])/alfa;
Lambda_var[5][i*4+3][j*4+3]=(0.5*dlambda_dX_fi[2]*dlambda_dX_f[0]+0.5*dlambda_dX_fi[0]*dlambda_dX_f[2])/alfa;

Lambda_var[0][i*4+3][j*4+3]=Lambda_var[0][i*4+3][j*4+3]+(1.0*dlambda_dX_pfi[0]*dlambda_dX_pf[0])/beta;
Lambda_var[1][i*4+3][j*4+3]=Lambda_var[1][i*4+3][j*4+3]+(1.0*dlambda_dX_pfi[1]*dlambda_dX_pf[1])/beta;
Lambda_var[2][i*4+3][j*4+3]=Lambda_var[2][i*4+3][j*4+3]+(1.0*dlambda_dX_pfi[2]*dlambda_dX_pf[2])/beta;
Lambda_var[3][i*4+3][j*4+3]=Lambda_var[3][i*4+3][j*4+3]+(0.5*dlambda_dX_pfi[0]*dlambda_dX_pf[1]+0.5*dlambda_dX_pfi[1]*dlambda_dX_pf[0])/beta;
Lambda_var[4][i*4+3][j*4+3]=Lambda_var[4][i*4+3][j*4+3]+(0.5*dlambda_dX_pfi[1]*dlambda_dX_pf[2]+0.5*dlambda_dX_pfi[2]*dlambda_dX_pf[1])/beta;
Lambda_var[5][i*4+3][j*4+3]=Lambda_var[5][i*4+3][j*4+3]+(0.5*dlambda_dX_pfi[2]*dlambda_dX_pf[0]+0.5*dlambda_dX_pfi[0]*dlambda_dX_pf[2])/beta;
*/
for(int k=0;k<3;k++)
{
dlambda_dX_f[k]=0.0;
dlambda_dX_fi[k]=0.0;
dlambda_dX_pf[k]=0.0;
dlambda_dX_pfi[k]=0.0;
}


} 

/*
double Kgeo_temp[24][24];

for(int i=0;i<24;i++)
for(int k=0;k<24;k++)
Kgeo_temp[i][k]=0.0;

for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
{
if(i%3==0 && j%3==0)
Kgeo_temp[i][j]=Kgeo_temp[i][j]+
(dNdX[i]*CauchyInv[0]*dNdX[j]
+dNdX[i+1]*CauchyInv[1]*dNdX[j+1]
+dNdX[i+2]*CauchyInv[2]*dNdX[j+2]
+dNdX[i]*CauchyInv[3]*dNdX[j+1]
+dNdX[i+1]*CauchyInv[3]*dNdX[j]
+dNdX[i+1]*CauchyInv[4]*dNdX[j+2]
+dNdX[i+2]*CauchyInv[4]*dNdX[j+1]
+dNdX[i]*CauchyInv[5]*dNdX[j+2]
+dNdX[i+2]*CauchyInv[5]*dNdX[j]);//*detj; 
else if(i%3!=0 && j>0)
Kgeo_temp[i][j]=Kgeo_temp[i-1][j-1]; 
}

*/






for(int i=0;i<32;i++)
for(int k=0;k<6;k++)
for(int j=0;j<32;j++)
{
_K_M2lambda[j][i]=_K_M2lambda[j][i]+1.0*current_nf*current_nf*detF*Lambda_var[k][i][j]*CauchyInv[k];
//_K_M2lambda[i][j]=_K_M2lambda[j][i];
}
/*
for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
{
if(i%4==3 && j%4==3)
_K_M2lambda[j][i]=_K_M2lambda[j][i]+
1.0/beta*current_nf*current_nf*detF*
(dNdX[i]*CauchyInv[0]*dNdX[j]
+dNdX[i+1]*CauchyInv[1]*dNdX[j+1]
+dNdX[i+2]*CauchyInv[2]*dNdX[j+2]
+dNdX[i]*CauchyInv[3]*dNdX[j+1]
+dNdX[i+1]*CauchyInv[3]*dNdX[j]
+dNdX[i+1]*CauchyInv[4]*dNdX[j+2]
+dNdX[i+2]*CauchyInv[4]*dNdX[j+1]
+dNdX[i]*CauchyInv[5]*dNdX[j+2]
+dNdX[i+2]*CauchyInv[5]*dNdX[j]);//*detj; 
//else if(i%4!=0 && j>0)
//_K_M2lambda[i][j]=_K_M2lambda[i-1][j-1]; 
}
*/

/*
double Ktrans[32][32];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
Ktrans[j][i]=_K_M2lambda[i][j];

for(int i=0;i<32;i++)
for(int j=0;j<32;j++)
_K_M2lambda[i][j]=0.5*(_K_M2lambda[i][j]+Ktrans[i][j]);
*/



/*
for(int i=0;i<24;i++)
for(int k=0;k<24;k++)
_K_M2lambda[i+i/3][k+k/3]=_K_M2lambda[i+i/3][k+k/3]+1.0/beta*current_nf*current_nf/detF*Kgeo_temp[i][k];
*/


};void calc_internal_M1(Node _nodes[],
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
		      double _coeff)
{

double beta_2=0.25;
double beta_1=0.5;
double Ftens[3][3];
double FtensT[3][3];
double Cauchytens[3][3];
double CauchytensInv[3][3];
double Cauchy_rate[6];
double CauchyInv[6];
double dNdX[24]; 
double detF;
double nodal_u_only[32];
double N_lambda[32];

//can: startof:=these values are for the geometric stiffnes similar part and its convergence 
double dNdksi[24]; 
double inversejacob[24][24]; 
double jacob[3][3]; 
//can: endof:=these values are for the geometric stiffnes similar part and its convergence 






Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];

transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchytens); 
inverse(Cauchytens,CauchytensInv);
detF=det(Ftens);




CauchyInv[0]=CauchytensInv[0][0];
CauchyInv[1]=CauchytensInv[1][1];
CauchyInv[2]=CauchytensInv[2][2];
CauchyInv[3]=CauchytensInv[0][1];
CauchyInv[4]=CauchytensInv[1][2];
CauchyInv[5]=CauchytensInv[2][0];

/*
Cauchy_rate[0]=1.0*_coeff*(Cauchytens[0][0]-1.0);
Cauchy_rate[1]=1.0*_coeff*(Cauchytens[1][1]-1.0);
Cauchy_rate[2]=1.0*_coeff*(Cauchytens[2][2]-1.0);
Cauchy_rate[3]=1.0*2.0*_coeff*Cauchytens[0][1];
Cauchy_rate[4]=1.0*2.0*_coeff*Cauchytens[1][2];
Cauchy_rate[5]=1.0*2.0*_coeff*Cauchytens[2][0];
*/

Cauchy_rate[0]=1.0*_coeff*(Cauchytens[0][0]-1.0);
Cauchy_rate[1]=1.0*_coeff*(Cauchytens[1][1]-1.0);
Cauchy_rate[2]=1.0*_coeff*(Cauchytens[2][2]-1.0);
Cauchy_rate[3]=1.0*_coeff*Cauchytens[0][1];
Cauchy_rate[4]=1.0*_coeff*Cauchytens[1][2];
Cauchy_rate[5]=1.0*_coeff*Cauchytens[2][0];






//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N_lambda[shape_counter+0]=0.0;
N_lambda[shape_counter+1]=0.0;
N_lambda[shape_counter+2]=0.0;
N_lambda[shape_counter+3]=0.125*(1+sign_ksi*_ksi)*(1+sign_eta*_eta)*(1+sign_zeta*_zeta);
shape_counter=shape_counter+4;
}
double temp=N_lambda[11];
N_lambda[11]=N_lambda[15];
N_lambda[15]=temp;
temp=N_lambda[27];
N_lambda[27]=N_lambda[31];
N_lambda[31]=temp;
/*
for(int n=0;n<4;n++)
{
double temp=N_lambda[2*4+n];
N_lambda[2*4+n]=N_lambda[3*4+n];
N_lambda[3*4+n]=temp;
temp=N_lambda[6*4+n];
N_lambda[6*4+n]=N_lambda[7*4+n];
N_lambda[7*4+n]=temp;
}
*/
//end definition shape value vector



double somecoeff=0.0;
for(int n=0;n<6;n++)
somecoeff=somecoeff+Cauchy_rate[n]*CauchyInv[n]*detF; //it was a division bevore 

for(int i=0;i<32;i++)
_internal_M1[i]=_internal_M1[i]+somecoeff*N_lambda[i];




};








void calc_internal_M2(Node _nodes[],
	              Element _elements[],
	              double _ksi, double _eta, double _zeta,
         	      double B[][32], //this is necessary for the computation of the green-lagrange strain tensor and the second piola kirchoff stress tensor 
	              double D[][6], //this is necessary for the computation of the second piola kirchoff stress tensor 
                      double F[],
	              double S[],
	              double _internal_M2[32],//internalforce of the conservation of mass 1 
	              double& detj,
                      int _currentelement,
                      ofstream &_debugfile,
	              double _t,
	              double _current_layup)
{



double beta_2=0.25;
double beta_1=0.5;
double Ftens[3][3];
double Cauchytens[3][3];
double CauchytensInv[3][3];
double CauchyInv[6];

double FtensT[3][3];
double dNdX[24]; 
double detF;
double N[8];
double beta=0.0;
double alfa=0.0;

//can: startof:=these values are for the geometric stiffnes similar part and its convergence 
double dNdksi[24]; 
double inversejacob[24][24]; 
double jacob[3][3]; 
//can: endof:=these values are for the geometric stiffnes similar part and its convergence 
double dlambda_dX[3];
double dlambda_dX_f[3];
double dlambda_dX_pf[3];
double M_vec[3];
double M_met[3][3];
double Lambda_var[6][32];
double current_nS=0.0;
double current_nf;

M_vec[0]=1.0;
M_vec[1]=0.0;
M_vec[2]=0.0;


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
M_met[i][j]=M_vec[i]*M_vec[j];

Ftens[0][0]=F[0];
Ftens[0][1]=F[3];
Ftens[0][2]=F[6];
Ftens[1][0]=F[1];
Ftens[1][1]=F[4];
Ftens[1][2]=F[7];
Ftens[2][0]=F[2];
Ftens[2][1]=F[5];
Ftens[2][2]=F[8];

transpose(Ftens,FtensT); 
matrix_product(FtensT,Ftens,Cauchytens); 
inverse(Cauchytens,CauchytensInv);
detF=det(Ftens);

CauchyInv[0]=CauchytensInv[0][0];
CauchyInv[1]=CauchytensInv[1][1];
CauchyInv[2]=CauchytensInv[2][2];
CauchyInv[3]=CauchytensInv[0][1];
CauchyInv[4]=CauchytensInv[1][2];
CauchyInv[5]=CauchytensInv[2][0];






//can: starting to fill in the shape value vector 
int shape_counter=0;
for(double sign_zeta=-1.0;sign_zeta<2.0;sign_zeta=sign_zeta+2.0)
for(double sign_eta=-1.0;sign_eta<2.0;sign_eta=sign_eta+2.0)
for(double sign_ksi=-1.0;sign_ksi<2.0;sign_ksi=sign_ksi+2.0)
{
N[shape_counter]=0.125*(1+sign_ksi*_ksi)*(1+sign_eta*_eta)*(1+sign_zeta*_zeta);
shape_counter++;
}
double temp=N[2];
N[2]=N[3];
N[3]=temp;
temp=N[6];
N[6]=N[7];
N[7]=temp;
//end definition shape value vector


//impermeability in perpendicular fiber direction 
for(int n=0;n<8;n++)
beta=beta+N[n]*_nodes[_elements[_currentelement].connid[n]].inpermeability;


//impermeability in fiber direction
alfa=beta/100.0;

/*
if(_current_layup==1.0)
{
alfa=1000.0;
beta=1000.0;
}
*/

for(int s=0;s<8;s++)
current_nS=current_nS+N[s]*_nodes[_elements[_currentelement].connid[s]].nS;
current_nf=1.0-current_nS/detF;





//computation of the third part 
dNdksi[0]=-0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[1]=-0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[2]=-0.125*(1.0-_ksi)*(1.0-_eta);


dNdksi[3]=0.125*(1.0-_eta)*(1.0-_zeta);
dNdksi[4]=-0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[5]=-0.125*(1.0+_ksi)*(1.0-_eta);


dNdksi[6]=0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[7]=0.125*(1.0+_ksi)*(1.0-_zeta);
dNdksi[8]=-0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[9]=-0.125*(1.0+_eta)*(1.0-_zeta);
dNdksi[10]=0.125*(1.0-_ksi)*(1.0-_zeta);
dNdksi[11]=-0.125*(1.0-_ksi)*(1.0+_eta);

dNdksi[12]=-0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[13]=-0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[14]=0.125*(1.0-_ksi)*(1.0-_eta);

dNdksi[15]=0.125*(1.0-_eta)*(1.0+_zeta);
dNdksi[16]=-0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[17]=0.125*(1.0+_ksi)*(1.0-_eta);

dNdksi[18]=0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[19]=0.125*(1.0+_ksi)*(1.0+_zeta);
dNdksi[20]=0.125*(1.0+_ksi)*(1.0+_eta);

dNdksi[21]=-0.125*(1.0+_eta)*(1.0+_zeta);
dNdksi[22]=0.125*(1.0-_ksi)*(1.0+_zeta);
dNdksi[23]=0.125*(1.0-_ksi)*(1.0+_eta);




for(int i=0;i<24;i++)
{
for(int j=0;j<24;j++)
inversejacob[i][j]=0.0; 
dNdX[i]=0.0; 
}


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
jacob[i][j]=0.0;


for(int i=0;i<3;i++) // row label 
for(int j=0;j<8;j++) 
{
jacob[i][0]=jacob[i][0]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[0];
jacob[i][1]=jacob[i][1]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[1];
jacob[i][2]=jacob[i][2]+dNdksi[j*3+i]*_nodes[_elements[_currentelement].connid[j]].X_lamb[2];
}

detj=jacob[0][0]*(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])-
     jacob[0][1]*(jacob[1][0]*jacob[2][2]-jacob[2][0]*jacob[1][2])+
     jacob[0][2]*(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1]);


inversejacob[0][0]=inversejacob[3][3]=inversejacob[6][6]=inversejacob[9][9]=inversejacob[12][12]=inversejacob[15][15]=inversejacob[18][18]=inversejacob[21][21]=(jacob[1][1]*jacob[2][2]-jacob[2][1]*jacob[1][2])/detj;
inversejacob[0][1]=inversejacob[3][4]=inversejacob[6][7]=inversejacob[9][10]=inversejacob[12][13]=inversejacob[15][16]=inversejacob[18][19]=inversejacob[21][22]=(jacob[0][2]*jacob[2][1]-jacob[2][2]*jacob[0][1])/detj;
inversejacob[0][2]=inversejacob[3][5]=inversejacob[6][8]=inversejacob[9][11]=inversejacob[12][14]=inversejacob[15][17]=inversejacob[18][20]=inversejacob[21][23]=(jacob[0][1]*jacob[1][2]-jacob[1][1]*jacob[0][2])/detj; 

inversejacob[1][0]=inversejacob[4][3]=inversejacob[7][6]=inversejacob[10][9]=inversejacob[13][12]=inversejacob[16][15]=inversejacob[19][18]=inversejacob[22][21]=(jacob[1][2]*jacob[2][0]-jacob[2][2]*jacob[1][0])/detj;
inversejacob[1][1]=inversejacob[4][4]=inversejacob[7][7]=inversejacob[10][10]=inversejacob[13][13]=inversejacob[16][16]=inversejacob[19][19]=inversejacob[22][22]=(jacob[0][0]*jacob[2][2]-jacob[2][0]*jacob[0][2])/detj;
inversejacob[1][2]=inversejacob[4][5]=inversejacob[7][8]=inversejacob[10][11]=inversejacob[13][14]=inversejacob[16][17]=inversejacob[19][20]=inversejacob[22][23]=(jacob[0][2]*jacob[1][0]-jacob[1][2]*jacob[0][0])/detj;

inversejacob[2][0]=inversejacob[5][3]=inversejacob[8][6]=inversejacob[11][9]=inversejacob[14][12]=inversejacob[17][15]=inversejacob[20][18]=inversejacob[23][21]=(jacob[1][0]*jacob[2][1]-jacob[2][0]*jacob[1][1])/detj;
inversejacob[2][1]=inversejacob[5][4]=inversejacob[8][7]=inversejacob[11][10]=inversejacob[14][13]=inversejacob[17][16]=inversejacob[20][19]=inversejacob[23][22]=(jacob[0][1]*jacob[2][0]-jacob[2][1]*jacob[0][0])/detj;
inversejacob[2][2]=inversejacob[5][5]=inversejacob[8][8]=inversejacob[11][11]=inversejacob[14][14]=inversejacob[17][17]=inversejacob[20][20]=inversejacob[23][23]=(jacob[0][0]*jacob[1][1]-jacob[1][0]*jacob[0][1])/detj;


for(int i=0;i<24;i++)
for(int j=0;j<24;j++)
dNdX[i]=dNdX[i]+inversejacob[i][j]*dNdksi[j];




for(int i=0;i<3;i++)
{
dlambda_dX[i]=0.0;
dlambda_dX_f[i]=0.0;
dlambda_dX_pf[i]=0.0;
}



for(int i=0;i<24;i++) 
dlambda_dX[i%3]=dlambda_dX[i%3]+dNdX[i]*_nodes[_elements[_currentelement].connid[i/3]].u_lamb[3];


for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
dlambda_dX_f[i]=dlambda_dX_f[i]+M_met[i][j]*dlambda_dX[j];

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
dlambda_dX_pf[i]=dlambda_dX_pf[i]+(cronecker(i,j)-M_met[i][j])*dlambda_dX[j];



for(int i=0;i<6;i++)
for(int j=0;j<32;j++)
Lambda_var[i][j]=0.0;




for(int i=0;i<8;i++)
{
//the first three columns are set to be zero 
Lambda_var[0][i*4+3]=(1.0*dNdX[i*3+0]*dlambda_dX_f[0])/alfa;
Lambda_var[1][i*4+3]=(1.0*dNdX[i*3+1]*dlambda_dX_f[1])/alfa;
Lambda_var[2][i*4+3]=(1.0*dNdX[i*3+2]*dlambda_dX_f[2])/alfa;
Lambda_var[3][i*4+3]=(0.5*dNdX[i*3+0]*dlambda_dX_f[1]+0.5*dNdX[i*3+1]*dlambda_dX_f[0])/alfa;
Lambda_var[4][i*4+3]=(0.5*dNdX[i*3+1]*dlambda_dX_f[2]+0.5*dNdX[i*3+2]*dlambda_dX_f[1])/alfa;
Lambda_var[5][i*4+3]=(0.5*dNdX[i*3+2]*dlambda_dX_f[0]+0.5*dNdX[i*3+0]*dlambda_dX_f[2])/alfa;

Lambda_var[0][i*4+3]=Lambda_var[0][i*4+3]+(1.0*dNdX[i*3+0]*dlambda_dX_pf[0])/beta;
Lambda_var[1][i*4+3]=Lambda_var[1][i*4+3]+(1.0*dNdX[i*3+1]*dlambda_dX_pf[1])/beta;
Lambda_var[2][i*4+3]=Lambda_var[2][i*4+3]+(1.0*dNdX[i*3+2]*dlambda_dX_pf[2])/beta;
Lambda_var[3][i*4+3]=Lambda_var[3][i*4+3]+(0.5*dNdX[i*3+0]*dlambda_dX_pf[1]+0.5*dNdX[i*3+1]*dlambda_dX_pf[0])/beta;
Lambda_var[4][i*4+3]=Lambda_var[4][i*4+3]+(0.5*dNdX[i*3+1]*dlambda_dX_pf[2]+0.5*dNdX[i*3+2]*dlambda_dX_pf[1])/beta;
Lambda_var[5][i*4+3]=Lambda_var[5][i*4+3]+(0.5*dNdX[i*3+2]*dlambda_dX_pf[0]+0.5*dNdX[i*3+0]*dlambda_dX_pf[2])/beta;

} 







for(int i=0;i<32;i++)
for(int j=0;j<6;j++)
_internal_M2[i]=_internal_M2[i]+current_nf*current_nf*detF*CauchyInv[j]*Lambda_var[j][i]; //detF was division before 

};











void set_layup_integration(Node _nodes[], Element _elements[], int _numberofnodes, int _numberofelements)
{
for(int e=0;e<_numberofelements;e++)
{
int current_counter=0;
for(double g3=-1.0;g3<2.0;g3=g3+2.0)
for(double g2=-1.0;g2<2.0;g2=g2+2.0)
for(double g1=-1.0;g1<2.0;g1=g1+2.0)
{
int current_node=current_counter;
if(current_counter==3)
current_node=4;
if(current_counter==4)
current_node=3;
if(current_counter==7)
current_node=8;
if(current_counter==8)
current_node=7;
_elements[e].layup[current_counter]=_nodes[_elements[e].connid[current_node]].layup;

current_counter++;
}
}

};



void apply_inpermeability_bc(Node _nodes[], ifstream& _inpermeability_bc, int _numberofnodes)
{


char star='+'; 
int currentnodelabel;
int currentid; 
do
{

	_inpermeability_bc>>star; 
	_inpermeability_bc.putback(star); 
	_inpermeability_bc>>currentnodelabel; 
for(int n=0;n<_numberofnodes;n++)
if(currentnodelabel==_nodes[n].nodelabel)
currentid=n;

	if(star!='*')
	{
	//_nodes[currentid].inpermeability=20.0;
	_nodes[currentid].layup=1.0;
	}

}while(star!='*'); 
 }; 






void apply_natural_bc(Node _nodes[],ifstream& _natural_bc,int _numberofnodes, double _u[], ofstream & _debugfile)
{
char star='+'; 
int currentnodelabel; 
int currentid; 
do
{
	_natural_bc>>star; 
	_natural_bc.putback(star); 
	_natural_bc>>currentnodelabel; 
for(int n=0;n<_numberofnodes;n++)
if(currentnodelabel==_nodes[n].nodelabel)
currentid=n;

	if(star!='*')
	{
		_natural_bc>>_nodes[currentid].natural_bc[0]; 
		_natural_bc>>_nodes[currentid].natural_bc[1];
		_natural_bc>>_nodes[currentid].natural_bc[2];
	        _natural_bc>>_nodes[currentid].natural_bc[3];

		_natural_bc>>_u[currentid*4+0];
                _natural_bc>>_u[currentid*4+1];
 		_natural_bc>>_u[currentid*4+2];
 		_natural_bc>>_u[currentid*4+3];

		
                _nodes[currentid].natural_val[0]=_u[currentid*4+0];
		_nodes[currentid].natural_val[1]=_u[currentid*4+1];
		_nodes[currentid].natural_val[2]=_u[currentid*4+2];
		_nodes[currentid].natural_val[3]=_u[currentid*4+3];


//pre assignment of the neumann boundary conditions 
for(int i=0;i<4;i++)
if(_nodes[currentid].natural_val[i]!=0.0)
{
_nodes[currentid].nonzero_natural_bc[i]=1;
//_nodes[currentid].neumann_bc[i]=1;
}


	}

}while(star!='*'); 


//_debugfile<<"u_after_reading"<<endl;
//for(int i=0;i<4*_numberofnodes;i++)
//_debugfile<<_u[i]<<'\t'; 
//_debugfile<<endl; 


};  


void apply_neumann_bc(double _force[], Node _nodes[], ifstream& _force_bc, int _numberofnodes)
{

for(int i=0;i<_numberofnodes*3;i++)
_force[i]=0.0; 

char star='+'; 
int currentnodelabel;
int currentid; 
do
{
	_force_bc>>star; 
	_force_bc.putback(star); 
	_force_bc>>currentnodelabel; 
for(int n=0;n<_numberofnodes;n++)
if(currentnodelabel==_nodes[n].nodelabel)
currentid=n;

	if(star!='*')
	{
		_force_bc>>_nodes[currentid].neumann_bc[0]; 
		_force_bc>>_nodes[currentid].neumann_bc[1];
		_force_bc>>_nodes[currentid].neumann_bc[2];
		_force_bc>>_nodes[currentid].neumann_bc[3];


		_force_bc>>_force[currentid*4+0];
                _force_bc>>_force[currentid*4+1];
 		_force_bc>>_force[currentid*4+2];
 		_force_bc>>_force[currentid*4+3];

		
                _nodes[currentid].neumann_val[0]=_force[currentid*4+0];
		_nodes[currentid].neumann_val[1]=_force[currentid*4+1];
		_nodes[currentid].neumann_val[2]=_force[currentid*4+2];
		_nodes[currentid].neumann_val[3]=_force[currentid*4+3];


//		_force[currentid*4+3]=0.0;
//		_nodes[currentid].neumann_bc[3]=1; // in fact does nothing 
	}

}while(star!='*'); 
 }; 



void initialize(Node _nodes[], double _numberofloadsteps, ofstream& _debugfile, int _numberofnodes,double _du[], double _u[], double _incremental_force[], double _force[], double _l)
{
double beta_2=0.25;
double beta_1=0.5;

for(int i=0;i<4*_numberofnodes;i++)
{
if(_nodes[i/4].nonzero_natural_bc[i%4]==1)
{
_du[i]=_u[i]/_numberofloadsteps;
_nodes[i/4].du_lamb[i%4]=_du[i];
}
else 
_nodes[i/4].du_lamb[i%4]=0.0;

//if(_nodes[i/4].neumann_bc[i%4]==1)
//_incremental_force[i]=_incremental_force[i]+beta_2/beta_1/_l*(_du[i-1]+_du[i-2]+_du[i-3]);
//if(_nodes[i/4].neumann_bc[i%4]==1)
//_incremental_force[i]=-beta_2/beta_1/_l*(_incremental_force[i]*_l*beta_1/beta_2+_du[i-1]+_du[i-2]+_du[i-3]);
}

};


void K_global(Node _nodes[],
	      Element _elements[],
	      int _numberofelements,
              ofstream &_debugfile,
              double t)
{
for(int e=0;e<_numberofelements;e++)
K_element(_nodes,_elements,e,_debugfile,t);

}; 


void eftable(Element _elements[], int** _eft,ofstream& _debugfile,int _numberofelements)
{
for(int e=0;e<_numberofelements;e++)
for(int c=0;c<8;c++)
{
_eft[4*c][e]=_elements[e].connid[c]*4;
_eft[4*c+1][e]=_elements[e].connid[c]*4+1;
_eft[4*c+2][e]=_elements[e].connid[c]*4+2;
_eft[4*c+3][e]=_elements[e].connid[c]*4+3;
}
/*
debugfile<<"element freedom table:"<<endl;
for(int i=0;i<8;i++)
{
for(int j=0;j<2;j++)
debugfile<<eft[i][j]<<'\t';
debugfile<<endl; 
}
*/ 
};


void K_assembly(Element _elements[], double** _Kglobal, int** _eft,ofstream& _debugfile,int _numberofnodes, int _numberofelements)
{

for(int i=0;i<4*_numberofnodes;i++)
for(int j=0;j<4*_numberofnodes;j++)	
_Kglobal[i][j]=0.0; 	

for(int e=0;e<_numberofelements;e++)
for(int m=0;m<32;m++)
for(int n=0;n<32;n++)
{
_Kglobal[_eft[m][e]][_eft[n][e]]=_Kglobal[_eft[m][e]][_eft[n][e]]+_elements[e].Kelem[m][n]; 
//_elements[e].Kelem[m][n]=0.0;
}



/*
_debugfile<<"Unreduced Global stiffness matrix:"<<endl; 
_debugfile<<setw(10); 
for(int i=0;i<12;i++)
{
for(int j=0;j<12;j++)
_debugfile<<_Kglobal[i][j]<<setw(18); 
_debugfile<<endl; 
}
*/ 
}; 


void Kglobal_reducer(Node _nodes[], double** _Kglobal, double _fi[], ofstream& _debugfile, int _numberofnodes)
{
for(int n=0;n<_numberofnodes;n++)
{
if(_nodes[n].natural_bc[0]==1 && _nodes[n].natural_val[0]==0)
for(int i=0;i<4*_numberofnodes;i++)
{
_Kglobal[n*4][i]=0.0;
_Kglobal[i][n*4]=0.0;
_Kglobal[n*4][n*4]=1.0;
_fi[n*4]=0.0; 
}


if(_nodes[n].natural_bc[1]==1 && _nodes[n].natural_val[1]==0)
for(int i=0;i<4*_numberofnodes;i++)
{
_Kglobal[n*4+1][i]=0.0;
_Kglobal[i][n*4+1]=0.0;
_Kglobal[n*4+1][n*4+1]=1.0;
_fi[n*4+1]=0.0; 
}


if(_nodes[n].natural_bc[2]==1 && _nodes[n].natural_val[2]==0)
for(int i=0;i<4*_numberofnodes;i++)
{
_Kglobal[n*4+2][i]=0.0;
_Kglobal[i][n*4+2]=0.0;
_Kglobal[n*4+2][n*4+2]=1.0;
_fi[n*4+2]=0.0; 
}


if(_nodes[n].natural_bc[3]==1 && _nodes[n].natural_val[3]==0)
for(int i=0;i<4*_numberofnodes;i++)
{
_Kglobal[n*4+3][i]=0.0;
_Kglobal[i][n*4+3]=0.0;
_Kglobal[n*4+3][n*4+3]=1.0;
_fi[n*4+3]=0.0; 
}




}



};



void Kglobal_ff_reducer(Node _nodes[], double** _Kglobal, double _fi[], ofstream& _debugfile, int _numberofnodes)
{
for(int n=0;n<_numberofnodes;n++)
{

if(_nodes[n].neumann_bc[0]==1)
for(int i=0;i<4*_numberofnodes;i++)
{
_Kglobal[n*4][i]=0.0;
_Kglobal[i][n*4]=0.0;
_Kglobal[n*4][n*4]=1.0;
_fi[n*4]=0.0; 
}



if(_nodes[n].neumann_bc[1]==1)
for(int i=0;i<4*_numberofnodes;i++)
{
_Kglobal[n*4+1][i]=0.0;
_Kglobal[i][n*4+1]=0.0;
_Kglobal[n*4+1][n*4+1]=1.0;
_fi[n*4+1]=0.0; 
}



if(_nodes[n].neumann_bc[2]==1)
for(int i=0;i<4*_numberofnodes;i++)
{
_Kglobal[n*4+2][i]=0.0;
_Kglobal[i][n*4+2]=0.0;
_Kglobal[n*4+2][n*4+2]=1.0;
_fi[n*4+2]=0.0; 
}

//for this stage, namely, where the element stiffness matrices not embedded, we have to impose this everywhere 
//this below is only outcommented for the reducer_ff:freefree, because this matrix is the one which is solved 
//if(_nodes[n].neumann_bc[3]==1)
for(int i=0;i<4*_numberofnodes;i++)
{
_Kglobal[n*4+3][i]=0.0;
_Kglobal[i][n*4+3]=0.0;
_Kglobal[n*4+3][n*4+3]=1.0;
_fi[n*4+3]=0.0; 
}

}


};


void gausselim(double rhs[],double** K,int numberofdatapoints,ofstream& outfile,Node nodes[], double du[])
{
  double** M=new double*[numberofdatapoints]; 

	for (int i = 0; i < numberofdatapoints; i++) 
	{
		du[i]=rhs[i]; 
		M[i]=new double[numberofdatapoints];
		for (int j = 0; j < numberofdatapoints; j++) 
		{
		M[i][j]=K[i][j]; 
		}
	}
	
	// LU decomposition without pivoting 
  for (int k = 0; k < numberofdatapoints - 1; k++) {
    if (M[k][k] == 0) 
      cout<<"pivot is zero in Mtx::GaussElim()"<<endl;
    for (int i = k + 1; i < numberofdatapoints; i++) {
      if (M[i][k] != 0) {   // M[i][k] can be complex
        double mult = M[i][k]/M[k][k];
        M[i][k] = mult;
        for (int j = k + 1; j < numberofdatapoints; j++) 
          M[i][j] = M[i][j]-mult*M[k][j];
      }
    }
  }

  // forwad substitution for L y = b. y still stored in rhs
  for (int i = 1; i < numberofdatapoints; i++) 
  for (int j = 0; j < i; j++) du[i] =du[i]-M[i][j]*du[j];

  // back substitution for U x = y. x still stored in rhs
  for (int i = numberofdatapoints - 1; i >= 0; i--) {
  for (int j=i+1; j<numberofdatapoints; j++) du[i] = du[i]-M[i][j]*du[j];
  du[i] =du[i]/ M[i][i];
  }

	
  //outfile<<endl; 
  //outfile<<"nodal displacement vector:"<<endl; 
  for(int i=0;i<numberofdatapoints;i++)
  {


          if(nodes[i/4].neumann_bc[i%4]==0)
	  nodes[i/4].du_lamb[i%4]=du[i];
	  else 
	  nodes[i/4].du_lamb[i%4]=0.0;

// this one below is certainly to be investigated 
//          if(i%4==3 && nodes[i/4].neumann_bc[i%4]==0)
//           nodes[i/4].X_lamb[3]=nodes[i/4].X_lamb[3]+du[i];
  }


for(int i=0;i<numberofdatapoints;i++)
delete[] M[i]; 
delete[] M; 


};



void du_copy(double du[],double du_prev[],ofstream& debugfile,int numberofnodes)
{
for(int i=0;i<4*numberofnodes;i++)
du_prev[i]=du[i]; 
}; 



void fullincrementalforcevector(Node nodes[],double incrementalforcevector[], double du[],
				double** Kglobal, ofstream& debugfile, int numberofnodes, double fi[])
{

for(int i=0;i<4*numberofnodes;i++)
{
for(int j=0;j<4*numberofnodes;j++)
if(nodes[i/4].neumann_bc[i%4]==0 && nodes[j/4].neumann_bc[j%4]==1)
incrementalforcevector[i]=incrementalforcevector[i]+Kglobal[i][j]*du[j];  


if(nodes[i/4].neumann_bc[i%4]==0)
fi[i]=-incrementalforcevector[i]+fi[i];

}

};


void csrfullincrementalforcevector(Node nodes[],double incrementalforcevector[], double du[],
				  csrM K, ofstream& debugfile, int numberofnodes, double fi[])
{

for(int i=0;i<4*numberofnodes;i++)
{
for(int j=K.iM[i];j<K.iM[i+1];j++)
if(nodes[i/4].nonzero_natural_bc[i%4]==0 && nodes[K.jM[j]/4].nonzero_natural_bc[K.jM[j]%4]==1)
incrementalforcevector[i]=incrementalforcevector[i]-K.vM[j]*du[K.jM[j]];



//if(nodes[i/4].neumann_bc[i%4]==1)
//fi[i]=10.0+fi[i];


//if(nodes[i/4].neumann_bc[i%4]==0)
if(nodes[i/4].nonzero_natural_bc[i%4]==0)
fi[i]=incrementalforcevector[i]+fi[i];



}

debugfile<<"fi_fullinc"<<endl;
for(int i=0;i<4*numberofnodes;i++)
if(nodes[i/4].neumann_bc[i%4]==1)
debugfile<<fi[i]<<'\t'; 
debugfile<<endl; 


/*
debugfile<<"incremental_beginning"<<endl;
for(int k=0;k<4*numberofnodes;k++)
debugfile<<incrementalforcevector[k]<<'\t'; 
debugfile<<endl; 

debugfile<<"fi_beginning"<<endl;
for(int k=0;k<4*numberofnodes;k++)
debugfile<<fi[k]<<'\t'; 
debugfile<<endl; 
*/


};












double eulernorm_residual(double fi[],int numberofnodes)
{
double norm=0.0;
for(int i=0;i<4*numberofnodes;i++)
norm=norm+fi[i]*fi[i];

return sqrt(norm); 
};  


double error_indu(double du[],double du_prev[],ofstream& debugfile,int numberofnodes)
{
double max;
double maxprev;  

double norm_du=eulernorm_residual(du,numberofnodes);
double norm_du_prev=eulernorm_residual(du_prev,numberofnodes);



/*
for(int i=0;i<12;i++)
{
maxprev=fabs(du[i]/du_prev[i]); 
if(i==0)
max=maxprev; 
if(maxprev>max)
max=maxprev; 
}
return max; 
*/ 
//cout<<"norm_du: "<<norm_du<<endl;
//cout<<"norm_du_prev: "<<norm_du_prev<<endl;
//return fabs(norm_du-norm_du_prev)/norm_du_prev; 
return norm_du/norm_du_prev; 
}; 


void Nodes_disp_update(Node _nodes[], int _numberofnodes)
{
for(int i=0;i<4*_numberofnodes;i++)
_nodes[i/4].u_lamb[i%4]=_nodes[i/4].u_lamb[i%4]+_nodes[i/4].du_lamb[i%4];
};

void Nodes_disp_downdate(Node _nodes[], int _numberofnodes)
{
for(int i=0;i<4*_numberofnodes;i++)
_nodes[i/4].u_lamb[i%4]=_nodes[i/4].u_lamb[i%4]-_nodes[i/4].du_lamb[i%4];
};




void csrK_assembly(Element _elements[], csrM & _K,  int** _eft, int _numberofnodes, int _numberofelements)
{
double current_Kelem_node[4][4];
for(int e=0;e<_numberofelements;e++)
for(int m=0;m<32;m++)
for(int n=0;n<32;n++)
if(m%4==0 && n%4==0)
{
for(int i=0;i<4;i++)
for(int j=0;j<4;j++)
current_Kelem_node[i][j]= _elements[e].Kelem[m+i][n+j];
_K.Insert_Matrix_incsrM(current_Kelem_node, _eft[m][e], _eft[n][e]);
}


for(int e=0;e<_numberofelements;e++)
for(int m=0;m<32;m++)
for(int n=0;n<32;n++)
_elements[e].Kelem[m][n]=0.0;
};





void csrKglobal_reducer(Node _nodes[], csrM &_K, double _fi[], ofstream& _debugfile, int _numberofnodes)
{


vector<int> rowscolumns;
for(int n=0;n<_numberofnodes;n++)
{
for(int i=0;i<4;i++)
if(_nodes[n].natural_bc[i]==1 && _nodes[n].natural_val[i]==0)
{
rowscolumns.push_back(n*4+i);
_fi[n*4+i]=0.0; 
}

//rowscolumns.push_back(n*4+3);
//_fi[n*4+3]=0.0;
}
_K.Shrink_Columns(rowscolumns);
_K.Shrink_Rows(rowscolumns);

//_K.writeCsrM_arrayform(_debugfile);

/*
for(int n=0;n<_numberofnodes;n++)
{
if(_nodes[n].natural_bc[0]==1 && _nodes[n].natural_val[0]==0)
{
_K.Shrink_Column(n*4);
_K.Shrink_Row(n*4);
_fi[n*4]=0.0; 
}

if(_nodes[n].natural_bc[1]==1 && _nodes[n].natural_val[1]==0)
{
_K.Shrink_Column(n*4+1);
_K.Shrink_Row(n*4+1);
_fi[n*4+1]=0.0; 
}


if(_nodes[n].natural_bc[2]==1 && _nodes[n].natural_val[2]==0)
{
_K.Shrink_Column(n*4+2);
_K.Shrink_Row(n*4+2);
_fi[n*4+2]=0.0; 
}


if(_nodes[n].natural_bc[3]==1 && _nodes[n].natural_val[3]==0)
{
_K.Shrink_Column(n*4+3);
_K.Shrink_Row(n*4+3);
_fi[n*4+3]=0.0; 
}




}
*/


};




void csrKglobal_ff_reducer(Node _nodes[], csrM &_K, double _fi[], ofstream& _debugfile, int _numberofnodes)
{


vector<int> rowscolumns;
for(int n=0;n<_numberofnodes;n++)
{
for(int i=0;i<4;i++)
if(_nodes[n].nonzero_natural_bc[i]==1)
{
rowscolumns.push_back(n*4+i);
_fi[n*4+i]=0.0; 
}
//rowscolumns.push_back(n*4+3);
}

_K.Shrink_Columns(rowscolumns);
_K.Shrink_Rows(rowscolumns);



/*
for(int n=0;n<_numberofnodes;n++)
for(int i=0;i<4;i++)
if(_nodes[n].neumann_bc[i]==1)
rowscolumns.push_back(n*4+i);

_K.Shrink_Rows(rowscolumns);
*/
//rowscolumns.push_back(n*4+3);
//_fi[n*4+3]=0.0; 


//_K.writeCsrM_arrayform(_debugfile);


/*
vector<int> additional;
for(int n=0;n<_numberofnodes;n++)
additional.push_back(n*4+3);
_K.Shrink_Rows(additional);
*/
/*
for(int n=0;n<_numberofnodes;n++)
{

if(_nodes[n].neumann_bc[0]==1)
{
_K.Shrink_Column(n*4);
_K.Shrink_Row(n*4);
_fi[n*4]=0.0; 
}



if(_nodes[n].neumann_bc[1]==1)
{
_K.Shrink_Column(n*4+1);
_K.Shrink_Row(n*4+1);
_fi[n*4+1]=0.0; 
}



if(_nodes[n].neumann_bc[2]==1)
{
_K.Shrink_Column(n*4+2);
_K.Shrink_Row(n*4+2);
_fi[n*4+2]=0.0; 
}

//for this stage, namely, where the element stiffness matrices not embedded, we have to impose this everywhere 
//this below is only outcommented for the reducer_ff:freefree, because this matrix is the one which is solved 
//if(_nodes[n].neumann_bc[3]==1)
{
_K.Shrink_Column(n*4+3);
_K.Shrink_Row(n*4+3);
_fi[n*4+3]=0.0; 
}

}
*/

};


void Nodes_neumann_update(Node nodes[], double du[], int numberofdatapoints) 
{
  for(int i=0;i<numberofdatapoints;i++)
  {

          if(nodes[i/4].nonzero_natural_bc[i%4]==0)
	  nodes[i/4].du_lamb[i%4]=du[i];
	  else 

	  nodes[i/4].du_lamb[i%4]=0.0;
	}
};


double det(double a[][3])
{
double result=0.0;
result=a[0][0]*(a[1][1]*a[2][2]-a[2][1]*a[1][2])-
       a[0][1]*(a[1][0]*a[2][2]-a[2][0]*a[1][2])+
       a[0][2]*(a[1][0]*a[2][1]-a[2][0]*a[1][1]);

return result;


};

double trace(double tens[3][3])
{
double result=0.0; 

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
if(i==j)
result=result+tens[i][j]; 

return result; 



}; 

void matrix_product(double tens1[3][3],double tens2[3][3],double result[3][3])
{
for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
result[i][j]=0.0; 

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
for(int k=0;k<3;k++)
result[i][j]=result[i][j]+tens1[i][k]*tens2[k][j]; 
}; 


void transpose(double tens[3][3], double tens_T[3][3])
{
for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
tens_T[i][j]=tens[j][i]; 
}; 

double cronecker(int i,int j)
{
double result; 
if(i==j)
result=1.0; 
else 
result=0.0; 
return result; 

}; 

double symmetry_tensor(int i,int j, int k, int l)
{
double result;
result=0.5*(cronecker(i,k)*cronecker(j,l)+cronecker(i,l)*cronecker(j,k));
return result;
};

void inverse(double tens[3][3], double inv_tens[3][3])
{
inv_tens[0][0]=tens[1][1]*tens[2][2]-tens[2][1]*tens[1][2];
inv_tens[0][1]=tens[0][2]*tens[2][1]-tens[2][2]*tens[0][1];
inv_tens[0][2]=tens[0][1]*tens[1][2]-tens[1][1]*tens[0][2];


inv_tens[1][0]=tens[1][2]*tens[2][0]-tens[2][2]*tens[1][0];
inv_tens[1][1]=tens[0][0]*tens[2][2]-tens[2][0]*tens[0][2];
inv_tens[1][2]=tens[0][2]*tens[1][0]-tens[1][2]*tens[0][0];


inv_tens[2][0]=tens[1][0]*tens[2][1]-tens[2][0]*tens[1][1];
inv_tens[2][1]=tens[0][1]*tens[2][0]-tens[2][1]*tens[0][0];
inv_tens[2][2]=tens[0][0]*tens[1][1]-tens[1][0]*tens[0][1];




for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
inv_tens[i][j]=inv_tens[i][j]/det(tens); 

};


void SetArrayZero(double _v[], int _m)
{
for(int i=0;i<_m;i++)
_v[i]=0.0; 
};



int node_counter(ifstream &_nodes_infile, int &_numberofnodes)
{
_numberofnodes=0;
char current_line[256];
_nodes_infile.getline(current_line,256);
int condition=0;
while(condition==0)
{
_nodes_infile.getline(current_line,256);
//cout<<current_line<<endl;
if(current_line[1]=='N' && current_line[2]=='O' && current_line[3]=='D' && current_line[4]=='E')
condition=1; 
}
_nodes_infile.getline(current_line,256);

while(current_line[0]!='*')
{
_nodes_infile.getline(current_line,256);
_numberofnodes++;
}


cout<<"The number of nodes in your model: "<<_numberofnodes<<endl;



return 0;
}

int element_counter(ifstream &_elements_infile, int &_numberofelements)
{
_numberofelements=0;
char current_line[256];
_elements_infile.getline(current_line,256);
int condition=0;
while(condition==0)
{
_elements_infile.getline(current_line,256);
//cout<<current_line<<endl;
if(current_line[1]=='E' && current_line[2]=='L' && current_line[3]=='E' && current_line[4]=='M' 
&& current_line[5]=='E' && current_line[6]=='N' && current_line[7]=='T')
condition=1; 
}
_elements_infile.getline(current_line,256);

while(current_line[0]!='*')
{
_elements_infile.getline(current_line,256);
_numberofelements++;
}


cout<<"The number of elements in your model: "<<_numberofelements<<endl;



return 0;
}


void result_writer(Node _nodes[], Element _elements[],
                  int _numberofnodes, int _numberofelements)
{


ofstream resultfile("example_result.inp");
resultfile<<"*NODE"<<endl;

double scale=1.0;

for(int n=0;n<_numberofnodes;n++)
{
/*
resultfile<<_nodes[n].nodelabel<<','<<setw(10);
resultfile<<_nodes[n].X[0]<<','<<setw(10);
resultfile<<_nodes[n].X[1]<<','<<setw(10);
resultfile<<_nodes[n].X[2]<<endl;
*/

resultfile<<_nodes[n].nodelabel<<','<<setw(10);
resultfile<<(_nodes[n].X_lamb[0])+scale*_nodes[n].u_lamb[0]<<','<<setw(10);
resultfile<<(_nodes[n].X_lamb[1])+scale*_nodes[n].u_lamb[1]<<','<<setw(10);
resultfile<<(_nodes[n].X_lamb[2])+scale*_nodes[n].u_lamb[2]<<endl;

}
resultfile<<"**"<<endl;



resultfile<<"*ELEMENT, TYPE=C3D8, ELSET=solids_example"<<endl;

for(int e=0;e<_numberofelements;e++)
{
resultfile<<_elements[e].elementlabel<<','<<setw(10);
for(int c=0;c<7;c++)
resultfile<<_elements[e].conname[c]<<','<<setw(10);
resultfile<<_elements[e].conname[7]<<endl;
}

resultfile<<"**"<<endl;

};
void odb_model_writer(Node _nodes[], Element _elements[],
                       int _numberofnodes, int _numberofelements,char *argv[],odb_Odb& odb)
{







  odb_Part& part1 = odb.Part("part-1",
       odb_Enum::THREE_D, odb_Enum::DEFORMABLE_BODY);

//writing the model data: nodes
  odb_SequenceInt nodeLabels;

for(int n=0;n<_numberofnodes;n++)
{
nodeLabels.append(_nodes[n].nodelabel);
}

  odb_SequenceSequenceFloat nodeCoor;
  for (int n=0; n<nodeLabels.size(); n++) {    
    odb_SequenceFloat loc;
    for (int i=0; i<3; i++)
      loc.append(_nodes[n].X_lamb[i]);
    nodeCoor.append(loc);
  }
  part1.addNodes(nodeLabels, nodeCoor, "Nodes");

//writing the model data: elements 
  odb_SequenceInt elLabels;

for(int e=0;e<_numberofelements;e++)
elLabels.append(_elements[e].elementlabel);

  odb_SequenceSequenceInt connect;
  const int numNodePerEl = 8;
  for (int e=0; e<elLabels.size(); e++) {
    odb_SequenceInt l;
    for (int i=0; i<numNodePerEl; i++)
      l.append(_elements[e].conname[i]);
    connect.append(l);
  }
	
    odb_String sectionCategoryName("C3D8");
    odb_String sectionCategoryDescription("my_hexahedras");
    odb_SectionCategory& sCat =
        odb.SectionCategory(sectionCategoryName,
                            sectionCategoryDescription);


  part1.addElements(elLabels, connect, "C3D8",
		    "Hexahedras",sCat);

    odb_Assembly& rootAssy = odb.rootAssembly();
    odb_Instance& instance1 =
    odb.rootAssembly().Instance("part-1-1", part1);

  // An ElementSet on an instance  
  odb_SequenceInt eLabelsA;
for(int e=0;e<_numberofelements;e++)
eLabelsA.append(_elements[e].elementlabel);
  instance1.ElementSet("elSetA", eLabelsA);
  
  // A NodeSet on the rootAssembly

  odb_SequenceSequenceInt nodeLabelsRA;
  odb_SequenceString namesRA;
  namesRA.append("part-1-1");
  odb_SequenceInt nodeLabelsRA_A;
  for (int n=0; n<nodeLabels.size(); n++)
  nodeLabelsRA_A.append(_nodes[n].nodelabel);
  nodeLabelsRA.append(nodeLabelsRA_A);
  const odb_Set& nSetRA = rootAssy.NodeSet("nodeSetRA",
					   namesRA, nodeLabelsRA);  


    odb_String stepName("step-1");
    odb_String stepDescription("first analysis step");
    odb_Step& step1 = odb.Step(stepName,
                               stepDescription,
                               odb_Enum::TIME,
			       1.0);

 };






void odb_result_writer(Node _nodes[], Element _elements[],
                       int _numberofnodes, int _numberofelements,char *argv[],odb_Odb& odb, double l, double numberofloadsteps)

{

odb_Assembly& rootAssy = odb.rootAssembly();

odb_Instance& instance1 = 
        rootAssy.instances()["part-1-1"];

    odb_Step& step1 = odb.steps()["step-1"];

    int incrementNumber = l+1;
    float analysisTime = (l+1)/numberofloadsteps;

    odb_String frameDescription("results frame for time");
    frameDescription.append(analysisTime);
    odb_Frame frame1 = step1.Frame(incrementNumber,
                                   analysisTime,
                                   frameDescription);


  odb_SequenceInt nodeLabels;
	for(int n=0;n<_numberofnodes;n++)
  	nodeLabels.append(_nodes[n].nodelabel);
  

  odb_SequenceInt eLabelsA;
	for(int e=0;e<_numberofelements;e++)
	eLabelsA.append(_elements[e].elementlabel);



// Field data:
    // Create a step and a frame.


    // Write nodal displacements.
    odb_String fieldName("U");
    odb_String fieldDescription("Displacements");
    odb_FieldOutput& uField =
        frame1.FieldOutput(fieldName,
                           fieldDescription,
			   odb_Enum::VECTOR);

    odb_SequenceSequenceFloat dispData;

    // create some displacement values
  for (int n=0; n<_numberofnodes; n++)
{
       odb_SequenceFloat current_dispData;
	current_dispData.append(_nodes[n].u_lamb[0]);
	current_dispData.append(_nodes[n].u_lamb[1]);
	current_dispData.append(_nodes[n].u_lamb[2]);
	dispData.append(current_dispData);
}   
    uField.addData(odb_Enum::NODAL,
                   instance1,
                   nodeLabels,
                   dispData);

    // Make this the default deformed field for visualization.

    step1.setDefaultDeformedField(uField);




    // Write nodal pore pressure 
    odb_String porepres("Pore-Pressure");
    odb_String porepresDescription("Porepressure");
    odb_FieldOutput& porepresField =
        frame1.FieldOutput(porepres,
                           porepresDescription,
			   odb_Enum::SCALAR);

    odb_SequenceSequenceFloat porepresData;

    // create some displacement values
  for (int n=0; n<_numberofnodes; n++)
{
        odb_SequenceFloat current_porepresData;
	current_porepresData.append(_nodes[n].u_lamb[3]);
	porepresData.append(current_porepresData);
}   
    porepresField.addData(odb_Enum::NODAL,
                         instance1,
                         nodeLabels,
                         porepresData);

//write integration solid fraction
   odb_String solidfraction("Solid-Fraction");
    odb_String sfDescription("Solidfraction");
    odb_FieldOutput& sfField =
        frame1.FieldOutput(solidfraction,
                           sfDescription,
			   odb_Enum::SCALAR);

    odb_SequenceSequenceFloat sfData;

    // create some displacement values
  for (int e=0; e<_numberofelements; e++)
	for(int i=0;i<8;i++)
{
        odb_SequenceFloat current_sfData;
	current_sfData.append(_elements[e].ns[i]);
	sfData.append(current_sfData);
}   
    sfField.addData(odb_Enum::INTEGRATION_POINT,
                         instance1,
                         eLabelsA,
                         sfData);

//write integration fluid fraction
   odb_String fluidfraction("Fluid-Fraction");
   odb_String ffDescription("Fluidfraction");
    odb_FieldOutput& ffField =
        frame1.FieldOutput(fluidfraction,
                           ffDescription,
			   odb_Enum::SCALAR);

    odb_SequenceSequenceFloat ffData;

    // create some displacement values
  for (int e=0; e<_numberofelements; e++)
	for(int i=0;i<8;i++)
{
        odb_SequenceFloat current_ffData;
	current_ffData.append(1.0-_elements[e].ns[i]);
	ffData.append(current_ffData);
}   
    ffField.addData(odb_Enum::INTEGRATION_POINT,
                         instance1,
                         eLabelsA,
                         ffData);


//write integration fluid velocity
   odb_String fluidvelocity("Fluid-Velocity");
   odb_String fvDescription("Fluidvelocity");
    odb_FieldOutput& fvField =
        frame1.FieldOutput(fluidvelocity,
                           fvDescription,
			   odb_Enum::VECTOR);

    odb_SequenceSequenceFloat fvData;

    // create some displacement values
  for (int e=0; e<_numberofelements; e++)
	for(int i=0;i<8;i++)
{
        odb_SequenceFloat current_fvData;
	current_fvData.append(_elements[e].v_fluid[i*3+0]);
	current_fvData.append(_elements[e].v_fluid[i*3+1]);
	current_fvData.append(_elements[e].v_fluid[i*3+2]);	
        fvData.append(current_fvData);
}   
    fvField.addData(odb_Enum::INTEGRATION_POINT,
                    instance1,
                    eLabelsA,
                    fvData);


//write integration seepage velocity
   odb_String seepagevelocity("Seepage-Velocity");
   odb_String spvDescription("Seepagevelocity");
    odb_FieldOutput& spvField =
        frame1.FieldOutput(seepagevelocity,
                           spvDescription,
			   odb_Enum::VECTOR);

    odb_SequenceSequenceFloat spvData;

    // create some displacement values
  for (int e=0; e<_numberofelements; e++)
	for(int i=0;i<8;i++)
{
        odb_SequenceFloat current_spvData;
	current_spvData.append(_elements[e].v_fs[i*3+0]);
	current_spvData.append(_elements[e].v_fs[i*3+1]);
	current_spvData.append(_elements[e].v_fs[i*3+2]);	
        spvData.append(current_spvData);
}   
    spvField.addData(odb_Enum::INTEGRATION_POINT,
                     instance1,
                     eLabelsA,
                     spvData);

//odb.save();
//odb.close();

};
