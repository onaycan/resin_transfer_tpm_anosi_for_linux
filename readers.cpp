


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