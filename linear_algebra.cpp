


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
