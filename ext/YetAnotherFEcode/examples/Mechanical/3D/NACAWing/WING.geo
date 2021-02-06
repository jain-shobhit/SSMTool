Include "AIRFOIL.geo";  // First Airfoil
SF=0.7; //Scale Factor
L=5;  //Length of wing
n = 25; // number of airfoil ribs
l = L/(n-1); // distance between each airfoil.
La =6; // number of layers in extrusion

DMN1=1;//Domain ID # 1
DMN2=2;//Domain ID # 2
DMN3=3;//Domain ID # 3
EBC=101;//EBC ID
NBC1=201;//NBC ID
NBC2=202;//NBC ID

Physical Surface(DMN1) = {}; 	// Curved Panel group
Physical Surface(DMN2) = {}; 	// longitunal stiffner Panel group
Physical Surface(DMN3) = {s1,s2,s3,s4,s5,s6}; 	// Ribs, Stiffeners group 
Physical Surface(EBC) = {s1,s2,s3,s4,s5,s6}; 	// Clamped surface 
Physical Surface(NBC1) = {}; 	// Force Area

For i In {1:(n-1)}
	If (i == 1)
		out[] = Extrude{0, 0, l} { Line{l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19}; Layers{La}; };
	Else
		out[] = Extrude{0, 0, l} { Line{out[0],out[4],out[8],out[12],out[16],out[20],out[24],out[28],out[32],out[36],out[40],out[44],out[48],out[52],out[56],out[60],out[64],out[68],out[72]}; Layers{La}; };
	EndIf
	
	//Surface loop 1	
	ll1 = newll;Line Loop(ll1) = {out[28], out[4], out[72], -out[0]};
	s1 = news;Plane Surface(s1) = {ll1};
	//Surface loop 2
	ll2 = newll;Line Loop(ll2) = {out[32], out[8], out[68], -out[4]};
	s2 = news;Plane Surface(s2) = {ll2};
	//Surface loop 3
	ll3 = newll;Line Loop(ll3) = {out[36], out[12], out[64], -out[8]};
	s3 = news;Plane Surface(s3) = {ll3};
	//Surface loop 4
	ll4 = newll;Line Loop(ll4) = {out[40], out[16], out[60], -out[12]};
	s4 = news;Plane Surface(s4) = {ll4};
	//Surface loop 5
	ll5 = newll;Line Loop(ll5) = {out[44], out[20], out[56], -out[16]};
	s5 = news;Plane Surface(s5) = {ll5};
	//Surface loop 6
	ll6 = newll;Line Loop(ll6) = {out[48], out[24], out[52], -out[20]};
	s6 = news;Plane Surface(s6) = {ll6};
	

	// Make surfaces transfinite for structured mesh
	Transfinite Line {out[0],out[4],out[8],out[12],out[16],out[20],out[24]} = a1 Using Progression 1;
	Transfinite Line {out[28],out[32],out[36],out[40],out[44],out[48],out[52],out[56],out[60],out[64],out[68],out[72]} = a2 Using Progression 1;
	Transfinite Surface {s1,s2,s3,s4,s5,s6};

	Physical Surface(DMN1) += {out[29],out[33],out[37],out[41],out[45],out[49],out[53],out[57],out[61],out[65],out[69],out[73]};
	Physical Surface(DMN2) += {out[1],out[5],out[9],out[13],out[17],out[21],out[25]};
	Physical Surface(DMN3) += {s1,s2,s3,s4,s5,s6};
	Physical Surface(NBC1) += {out[29]};
EndFor
Coherence;

