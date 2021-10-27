



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