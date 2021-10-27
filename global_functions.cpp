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