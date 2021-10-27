


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
