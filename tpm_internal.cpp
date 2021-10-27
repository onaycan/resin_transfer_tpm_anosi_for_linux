void calc_internal_M1(Node _nodes[],
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











