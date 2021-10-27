void odb_model_writer(Node _nodes[], Element _elements[],
                       int _numberofnodes, int _numberofelements,char *argv[],odb_Odb& odb)
{







  odb_Part& part1 = odb.Part("part-1",
       odb_Enum::THREE_D, odb_Enum::DEFORMABLE_BODY);

//writing the model data: nodes
  odb_SequenceInt nodeLabels;

for(int n=0;n<_numberofnodes;n++)
{
nodeLabels.append(_nodes[n].nodelabel);
}

  odb_SequenceSequenceFloat nodeCoor;
  for (int n=0; n<nodeLabels.size(); n++) {    
    odb_SequenceFloat loc;
    for (int i=0; i<3; i++)
      loc.append(_nodes[n].X_lamb[i]);
    nodeCoor.append(loc);
  }
  part1.addNodes(nodeLabels, nodeCoor, "Nodes");

//writing the model data: elements 
  odb_SequenceInt elLabels;

for(int e=0;e<_numberofelements;e++)
elLabels.append(_elements[e].elementlabel);

  odb_SequenceSequenceInt connect;
  const int numNodePerEl = 8;
  for (int e=0; e<elLabels.size(); e++) {
    odb_SequenceInt l;
    for (int i=0; i<numNodePerEl; i++)
      l.append(_elements[e].conname[i]);
    connect.append(l);
  }
	
    odb_String sectionCategoryName("C3D8");
    odb_String sectionCategoryDescription("my_hexahedras");
    odb_SectionCategory& sCat =
        odb.SectionCategory(sectionCategoryName,
                            sectionCategoryDescription);


  part1.addElements(elLabels, connect, "C3D8",
		    "Hexahedras",sCat);

    odb_Assembly& rootAssy = odb.rootAssembly();
    odb_Instance& instance1 =
    odb.rootAssembly().Instance("part-1-1", part1);

  // An ElementSet on an instance  
  odb_SequenceInt eLabelsA;
for(int e=0;e<_numberofelements;e++)
eLabelsA.append(_elements[e].elementlabel);
  instance1.ElementSet("elSetA", eLabelsA);
  
  // A NodeSet on the rootAssembly

  odb_SequenceSequenceInt nodeLabelsRA;
  odb_SequenceString namesRA;
  namesRA.append("part-1-1");
  odb_SequenceInt nodeLabelsRA_A;
  for (int n=0; n<nodeLabels.size(); n++)
  nodeLabelsRA_A.append(_nodes[n].nodelabel);
  nodeLabelsRA.append(nodeLabelsRA_A);
  const odb_Set& nSetRA = rootAssy.NodeSet("nodeSetRA",
					   namesRA, nodeLabelsRA);  


    odb_String stepName("step-1");
    odb_String stepDescription("first analysis step");
    odb_Step& step1 = odb.Step(stepName,
                               stepDescription,
                               odb_Enum::TIME,
			       1.0);

 };






void odb_result_writer(Node _nodes[], Element _elements[],
                       int _numberofnodes, int _numberofelements,char *argv[],odb_Odb& odb, double l, double numberofloadsteps)

{

odb_Assembly& rootAssy = odb.rootAssembly();

odb_Instance& instance1 = 
        rootAssy.instances()["part-1-1"];

    odb_Step& step1 = odb.steps()["step-1"];

    int incrementNumber = l+1;
    float analysisTime = (l+1)/numberofloadsteps;

    odb_String frameDescription("results frame for time");
    frameDescription.append(analysisTime);
    odb_Frame frame1 = step1.Frame(incrementNumber,
                                   analysisTime,
                                   frameDescription);


  odb_SequenceInt nodeLabels;
	for(int n=0;n<_numberofnodes;n++)
  	nodeLabels.append(_nodes[n].nodelabel);
  

  odb_SequenceInt eLabelsA;
	for(int e=0;e<_numberofelements;e++)
	eLabelsA.append(_elements[e].elementlabel);



// Field data:
    // Create a step and a frame.


    // Write nodal displacements.
    odb_String fieldName("U");
    odb_String fieldDescription("Displacements");
    odb_FieldOutput& uField =
        frame1.FieldOutput(fieldName,
                           fieldDescription,
			   odb_Enum::VECTOR);

    odb_SequenceSequenceFloat dispData;

    // create some displacement values
  for (int n=0; n<_numberofnodes; n++)
{
       odb_SequenceFloat current_dispData;
	current_dispData.append(_nodes[n].u_lamb[0]);
	current_dispData.append(_nodes[n].u_lamb[1]);
	current_dispData.append(_nodes[n].u_lamb[2]);
	dispData.append(current_dispData);
}   
    uField.addData(odb_Enum::NODAL,
                   instance1,
                   nodeLabels,
                   dispData);

    // Make this the default deformed field for visualization.

    step1.setDefaultDeformedField(uField);




    // Write nodal pore pressure 
    odb_String porepres("Pore-Pressure");
    odb_String porepresDescription("Porepressure");
    odb_FieldOutput& porepresField =
        frame1.FieldOutput(porepres,
                           porepresDescription,
			   odb_Enum::SCALAR);

    odb_SequenceSequenceFloat porepresData;

    // create some displacement values
  for (int n=0; n<_numberofnodes; n++)
{
        odb_SequenceFloat current_porepresData;
	current_porepresData.append(_nodes[n].u_lamb[3]);
	porepresData.append(current_porepresData);
}   
    porepresField.addData(odb_Enum::NODAL,
                         instance1,
                         nodeLabels,
                         porepresData);

//write integration solid fraction
   odb_String solidfraction("Solid-Fraction");
    odb_String sfDescription("Solidfraction");
    odb_FieldOutput& sfField =
        frame1.FieldOutput(solidfraction,
                           sfDescription,
			   odb_Enum::SCALAR);

    odb_SequenceSequenceFloat sfData;

    // create some displacement values
  for (int e=0; e<_numberofelements; e++)
	for(int i=0;i<8;i++)
{
        odb_SequenceFloat current_sfData;
	current_sfData.append(_elements[e].ns[i]);
	sfData.append(current_sfData);
}   
    sfField.addData(odb_Enum::INTEGRATION_POINT,
                         instance1,
                         eLabelsA,
                         sfData);

//write integration fluid fraction
   odb_String fluidfraction("Fluid-Fraction");
   odb_String ffDescription("Fluidfraction");
    odb_FieldOutput& ffField =
        frame1.FieldOutput(fluidfraction,
                           ffDescription,
			   odb_Enum::SCALAR);

    odb_SequenceSequenceFloat ffData;

    // create some displacement values
  for (int e=0; e<_numberofelements; e++)
	for(int i=0;i<8;i++)
{
        odb_SequenceFloat current_ffData;
	current_ffData.append(1.0-_elements[e].ns[i]);
	ffData.append(current_ffData);
}   
    ffField.addData(odb_Enum::INTEGRATION_POINT,
                         instance1,
                         eLabelsA,
                         ffData);


//write integration fluid velocity
   odb_String fluidvelocity("Fluid-Velocity");
   odb_String fvDescription("Fluidvelocity");
    odb_FieldOutput& fvField =
        frame1.FieldOutput(fluidvelocity,
                           fvDescription,
			   odb_Enum::VECTOR);

    odb_SequenceSequenceFloat fvData;

    // create some displacement values
  for (int e=0; e<_numberofelements; e++)
	for(int i=0;i<8;i++)
{
        odb_SequenceFloat current_fvData;
	current_fvData.append(_elements[e].v_fluid[i*3+0]);
	current_fvData.append(_elements[e].v_fluid[i*3+1]);
	current_fvData.append(_elements[e].v_fluid[i*3+2]);	
        fvData.append(current_fvData);
}   
    fvField.addData(odb_Enum::INTEGRATION_POINT,
                    instance1,
                    eLabelsA,
                    fvData);


//write integration seepage velocity
   odb_String seepagevelocity("Seepage-Velocity");
   odb_String spvDescription("Seepagevelocity");
    odb_FieldOutput& spvField =
        frame1.FieldOutput(seepagevelocity,
                           spvDescription,
			   odb_Enum::VECTOR);

    odb_SequenceSequenceFloat spvData;

    // create some displacement values
  for (int e=0; e<_numberofelements; e++)
	for(int i=0;i<8;i++)
{
        odb_SequenceFloat current_spvData;
	current_spvData.append(_elements[e].v_fs[i*3+0]);
	current_spvData.append(_elements[e].v_fs[i*3+1]);
	current_spvData.append(_elements[e].v_fs[i*3+2]);	
        spvData.append(current_spvData);
}   
    spvField.addData(odb_Enum::INTEGRATION_POINT,
                     instance1,
                     eLabelsA,
                     spvData);

//odb.save();
//odb.close();

};
