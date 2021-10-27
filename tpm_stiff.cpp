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


};