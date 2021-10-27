

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