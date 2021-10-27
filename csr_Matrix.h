#include <vector>
using namespace std;

struct csrM
{
vector<double> vM; //can: value of M 
vector<int> jM; //can: column index 
vector<int> iM; //can: row index 

void SetZero();
void CsrAddCsr(csrM L, csrM R, double aplha);
void CsrMultCsr(csrM L, csrM R, double _sparseprec);
void writeCsrM(); 
void writeCsrM(ofstream& outfile);
void writeCsrM_arrayform(ofstream& outfile);
void writeCsrM_arrayform();
void writeCsrM_row(int row, ofstream& outfile);
void ConvertToIdentity(int _m);
void CreateOfHollow(int _m);
void SortCsrM(double _saparseprec);
void KillCsrM();
void AsignCsrM(csrM right);
void MinusCsrM(); 
void MaxNormCsrM(double& _norm);
void MaxNormCsrM_local(double &_norm, int& vMloc, int iter);
void NeumannInverter(csrM _M, double _tol, char complement,
                     char killhim, double &_sparsity_loss, double _sparseprec);
void SparseSolverInverter(csrM _M, double _tol, char complement, char killhim,
                          double &_sparsity_loss, double _sparseprec);
void Pardiso_ExpInverter(csrM _K, double _tol, char complement, char killhim,
                         double &_sparsity_loss, double _sparseprec);
void Normalize_to_stochastic(ofstream& outfile);
void Push_Value_incsrM(double val, int row, int column);
void Insert_Value_incsrM(double val, int row, int column);
void Shrink_Column(int column);
void Shrink_Columns(vector<int> columns);
void Shrink_Row(int row);
void Shrink_Rows(vector<int> rows);
void ConvertFromCSC(csrM thisbutcsc);
void CsrMultVec(double _v[], double _result[],int _size);
void ScalarToTensorial(csrM Mscalar);
void Insert_Matrix_incsrM(double Kelem_node[4][4], int row, int column);
friend int Pardiso_UnsymSolver(csrM _K, double _lhs[], double _rhs[], int _size, char comment);

//void Insert_Matrix_incsrM(ublas::c_matrix<double,3,3> val, int row, int column);
};


