






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




