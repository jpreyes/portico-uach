function [] = meshgen()
clear all
close all
disp('************************************************');
disp('*               PROGRAM MESHGEN                *');
disp('*  MESH GENERATOR FOR TWO DIMENSIONAL REGIONS  *');
disp('*     (c) T.R.CHANDRUPATLA & A.D.BELEGUNDU     *');
disp('************************************************');

InputData
GlobalNode
CoordConnect
Output

%------------------------  function InputData  ---------------------------
function [] = InputData();
global NGCN NNS
global NEN NS NW NSJ NSW NGN NNW NSD NWD
global IDBLK MERG
global NODE 
global SR WR XB
global TITLE
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW

%FILE1 = input('Input Data File Name <DOS file name>','s');
FILE1 = input('Meshgen Input Data File Name (with extension) ','s');
LINP  = fopen(FILE1,'r');

% Read Dummy ---> Meshgen Inputfile
DUMMY = fgets(LINP);

% Read Title --->2D Quad (Title)
TITLE = fgets(LINP);

% Read Dummy ---> Element type(NEN): Triangular(3), Quatrilateral(4)
DUMMY = fgets(LINP);

% NEN = 3 for Triangle  4 for Quad
NEN = str2num(fgets(LINP));
if NEN < 3; NEN = 3; end;
if NEN > 4; NEN = 4; end;

%Hints:  A region is divided into 4-cornered blocks viewed as a
%        mapping from a Checkerboard pattern of S- and W- Sides
%        S- Side is one with lower number of final divisions
%        Blocks, Corners, S- and W- Sides are labeled as shown in Fig. 12.2
%        Make a sketch and identify void blocks and merging sides

% Read Dummy ---> BLOCK DATA
DUMMY = fgets(LINP);

% Read Dummy ---> #S-Spans(NS)  #W-Spans(NW)  #PairsOfEdgesMerged(NSJ)
DUMMY = fgets(LINP);

TMP = str2num(fgets(LINP));
[NS, NW, NSJ] = deal(TMP(1),TMP(2),TMP(3));

% Initial Data
NSW = NS * NW;
NGN = (NS + 1) * (NW + 1);
NM = 1;

IDBLK = zeros(NSW,1);
NSD = zeros(NS,1);
NWD = zeros(NW,1);
NGCN = zeros(NGN,1);
SH = zeros(8,1);

% Read Dummy ---> SPAN DIVISIONS
DUMMY = fgets(LINP);
NNS = 1; 
NNW = 1;

% Read Dummy ---> S-Direction
%DUMMY = fgets(LINP)

% Read Dummy ---> Number of divisions for each S-Span
DUMMY = fgets(LINP);

for KS = 1:NS
	TMP = str2num(fgets(LINP));
	[N, NSD(N)] = deal(TMP(1), TMP(2));
   NNS = NNS + NSD(N);
end

% Read Dummy ---> W-Direction
%DUMMY = fgets(LINP)

% Read Dummy ---> Number of divisions for each W-Span
DUMMY = fgets(LINP);

for KW = 1:NW
	TMP = str2num(fgets(LINP));
	[N, NWD(N)] = deal(TMP(1), TMP(2));
   NNW = NNW + NWD(N);
end

% Read Dummy ---> BLOCK MATERIAL DATA
DUMMY = fgets(LINP);
% Read Dummy ---> Block# 	Material Number(Void = 0)
DUMMY = fgets(LINP);


%-------- Block Identifier / Material# (Default# is 1) --------
for I = 1:NSW
   IDBLK(I) = 1;
end
while 1
   TMP = str2num(fgets(LINP));
   NTMP = TMP(1);
   if NTMP == 0; break; end;
   IDBLK(NTMP)= TMP(2);
   if NM < IDBLK(NTMP); NM = IDBLK(NTMP); end;
end

NSR = NS * (NW + 1);
NWR = NW * (NS + 1);

% Define Matrices
XB = zeros(NGN, 2);
SR = zeros(NSR, 2);
WR = zeros(NWR, 2);

% Read Dummy ---> BLOCK CORNER DATA
DUMMY = fgets(LINP);


% Read Dummy ---> Corner#	X-Coord.	Y-Coord.
DUMMY = fgets(LINP);


while 1
   TMP = str2num(fgets(LINP));
   NTMP = TMP(1);
   if NTMP == 0; break; end;
   XB(NTMP,:) = TMP(2:3);
end

%---------- Evaluate Mid-points of S-Sides -------------
for I = 1:NW + 1
   for J = 1:NS
      IJ = (I - 1) * NS + J;
      SR(IJ, 1) = 0.5 * (XB(IJ + I - 1, 1) + XB(IJ + I, 1));
      SR(IJ, 2) = 0.5 * (XB(IJ + I - 1, 2) + XB(IJ + I, 2));
   end
end

%---------- Evaluate Mid-points of W-Sides -------------
for I = 1:NW
   for J = 1:NS + 1
      IJ = (I - 1) * (NS + 1) + J;
      WR(IJ, 1) = 0.5 * (XB(IJ, 1) + XB(IJ + NS + 1, 1));
      WR(IJ, 2) = 0.5 * (XB(IJ, 2) + XB(IJ + NS + 1, 2));
   end
end

%------ Mid Points for Sides that are curved or graded ------
% Read Dummy --->MIDPOINT DATA
DUMMY = fgets(LINP);

% Read Dummy --->S-Side#	X-Coord. Y-Coord.
DUMMY = fgets(LINP);

%--- S-Sides
while 1
   TMP = str2num(fgets(LINP));
   NTMP = TMP(1);
   if NTMP == 0; break; end;
   SR(NTMP,:) = TMP(2:3);
end

% Read Dummy ---> MIDPOINT DATA
%DUMMY = fgets(LINP);
%DUMMY
% Read Dummy --->W-Side#	X-Coord. Y-Coord.
DUMMY = fgets(LINP);


%--- W-Sides
while 1
   TMP = str2num(fgets(LINP));
   NTMP = TMP(1);
   if NTMP == 0; break; end;
   WR(NTMP,:) = TMP(2:3);
end

%------- Global Node Locations of Corner Nodes ---------
NTMPI = 1;
for I = 1 : NW + 1
   if I == 1
      IINC = 0;
   else
      IINC = NNS * NWD(I - 1);
   end
   NTMPI = NTMPI + IINC;
   NTMPJ = 0;
   for J = 1 : NS + 1
      IJ = (NS + 1) * (I - 1) + J;
      if J == 1
         JINC = 0;
      else
         JINC = NSD(J - 1);
      end
      NTMPJ = NTMPJ + JINC;
      NGCN(IJ) = NTMPI + NTMPJ;
   end
end

%--------- Merging Sides ----------
if NSJ > 0
	% Read Dummy ---> MERGING SIDES
	DUMMY = fgets(LINP);

	% Read Dummy ---> 			SIDE 1			SIDE 2
	%DUMMY = fgets(LINP);
    %DUMMY
	% Read Dummy ---> Pair#	Node 1, Node 2	Node 1, Node 2
	DUMMY = fgets(LINP);

   MERG = zeros(NSJ, 4);
   for I = 1:NSJ
		TMP = str2num(fgets(LINP));
		[N, L1, L2, L3, L4] = deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5));
      IDIV1 = SideDiv(L1, L2);
      IDIV2 = SideDiv(L3, L4);
      if IDIV1 ~= IDIV2
         disp('#Div don''t match. Check merge data.');
         break;
      end
      MERG(I, 1) = L1;
      MERG(I, 2) = L2;
      MERG(I, 3) = L3; 
      MERG(I, 4) = L4;
   end     
end
        
fclose all;   

%------------------------  function SideDiv  ---------------------------
function  IDIV = SideDiv(I1, I2);
global NGCN NNS

%===========  Number of Divisions  for Side I1,I2  ===========
IMIN = I1;
IMAX = I2;
if IMIN > I2
   IMIN = I2;
   IMAX = I1;
end
if (IMAX - IMIN) == 1
   IDIV = NGCN(IMAX) - NGCN(IMIN);
else
   IDIV = (NGCN(IMAX) - NGCN(IMIN)) / NNS;
end

%------------------------  function GlobalNode  ---------------------------
function [] = GlobalNode();
global NGCN NNS
global NEN NS NW NSJ NSW NGN NNW NSD NWD
global IDBLK MERG
global NNAR
global X NOC MAT
global NODE

%---------------- Node Point Array --------------------
NNT = NNS * NNW;
for I = 1:NNT
   NNAR(I) = -1;
end
%--------- Zero Non-Existing Node Locations ---------
for KW = 1 : NW
   for KS = 1 : NS
      KSW = NS * (KW - 1) + KS;
      if IDBLK(KSW) <= 0
         %-------- Operation within an Empty Block --------
         K1 = (KW - 1) * (NS + 1) + KS;
         N1 = NGCN(K1);
         NS1 = 2; if KS == 1; NS1 = 1; end;
         NW1 = 2; if KW == 1; NW1 = 1; end;
         NS2 = NSD(KS) + 1;
         if KS < NS
            if IDBLK(KSW + 1) > 0; NS2 = NSD(KS); end;
         end
         NW2 = NWD(KW) + 1;
         if KW < NW
            if IDBLK(KSW + NS) > 0; NW2 = NWD(KW); end;
         end
         for I = NW1:NW2
            IN1 = N1 + (I - 1) * NNS;
            for J = NS1:NS2
               IJ = IN1 + J - 1;
               NNAR(IJ) = 0;
            end
         end
         ICT = 0;
         if (NS2 == NSD(KS))|(NW2 == NWD(KW)); ICT = 1; end;
         if (KS == NS) | (KW == NW); ICT = 1; end;
         if ICT == 0
            if IDBLK(KSW + NS + 1) > 0 ; NNAR(IJ) = -1; end;
         end 
      end
   end
end

%--------  Node Identification for Side Merging ------
if NSJ > 0
   for I = 1:NSJ
      I1 = MERG(I, 1);
      I2 = MERG(I, 2);
      IDIV = SideDiv(I1, I2);
      IA1 = NGCN(I1);
      IA2 = NGCN(I2);
      IASTP = (IA2 - IA1) / IDIV;
      I1 = MERG(I, 3);
      I2 = MERG(I, 4);
      IDIV = SideDiv(I1, I2);
      IB1 = NGCN(I1);
      IB2 = NGCN(I2);
      IBSTP = (IB2 - IB1) / IDIV;
      IAA = IA1 - IASTP;
      for IBB = IB1:IBSTP:IB2
         IAA = IAA + IASTP;
         if IBB == IAA
            NNAR(IAA) = -1;
         else
            NNAR(IBB) = IAA;
         end
      end
   end
end
 
%----------  Final Node Numbers in the Array  --------
NODE = 0;
for I = 1 : NNT
   if NNAR(I) > 0
      II = NNAR(I);
      NNAR(I) = NNAR(II);
   elseif NNAR(I) < 0
      NODE = NODE + 1;
      NNAR(I) = NODE;
   end
end

%------------------------  function CoordConnect  ---------------------------
function [] = CoordConnect();
global NGCN NNS
global NEN NS NW NSJ NSW NGN NNW NSD NWD
global IDBLK MERG
global NODE
global NNAR
global X NOC MAT
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
%------------  Nodal Coordinates  ---------------
NN = NODE; 
NELM = 0;

% X(NN, 2), XP(8, 2), NOC(2 * NNT, NEN), MAT(2 * NNT)
for KW = 1 : NW
   for KS = 1 : NS
      KSW = NS * (KW - 1) + KS;
      if IDBLK(KSW) ~= 0 
         %---------  Extraction of Block Data  ----------
         NODW = NGCN(KSW + KW - 1) - NNS - 1;
         for JW = 1 : NWD(KW) + 1
            ETA = -1 + 2 * (JW - 1) / NWD(KW);
            NODW = NODW + NNS;
            NODS = NODW;
            for JS = 1 : NSD(KS) + 1
               XI = -1 + 2 * (JS - 1) / NSD(KS);
               NODS = NODS + 1;
               NODE = NNAR(NODS);
               XP = BlockXY(KW, KSW);
               SH = Shape(XI, ETA);
               for J = 1 : 2
                  C1 = 0;
                  for I = 1 : 8
                     C1 = C1 + SH(I) * XP(I, J);
                  end
                  X(NODE, J) = C1;
               end
               %-----------------  Connectivity  ----------------
               if (JS ~= NSD(KS) + 1)&(JW ~= NWD(KW) + 1)
                  N1 = NODE;
                  N2 = NNAR(NODS + 1);
                  N3 = NNAR(NODS + NNS + 1);
                  N4 = NNAR(NODS + NNS);
                  NELM = NELM + 1;
                  if NEN == 3 
                     %------------- Triangular Elements ------------
                     NOC(NELM, 1) = N1;
                     NOC(NELM, 2) = N2;
                     NOC(NELM, 3) = N3;
                     MAT(NELM) = IDBLK(KSW);
                     NELM = NELM + 1; 
                     NOC(NELM, 1) = N3;
                     NOC(NELM, 2) = N4;
                     NOC(NELM, 3) = N1; 
                     MAT(NELM) = IDBLK(KSW);
                  else
                     %------------- Quadrilateral Elements ----------
                     NOC(NELM, 1) = N1; 
                     NOC(NELM, 2) = N2; 
                     MAT(NELM) = IDBLK(KSW);
                     NOC(NELM, 3) = N3; 
                     NOC(NELM, 4) = N4;
                  end
               end
            end
         end
      end
   end
end
NE = NELM;

if NEN == 3 
%--------- Readjustment for Triangle Connectivity ----------
	NE2 = NE / 2;
   for I = 1 : NE2
      I1 = 2 * I - 1;
      N1 = NOC(I1, 1); 
      N2 = NOC(I1, 2);
      N3 = NOC(I1, 3); 
      N4 = NOC(2 * I, 2);
      X13 = X(N1, 1) - X(N3, 1); 
      Y13 = X(N1, 2) - X(N3, 2);
      X24 = X(N2, 1) - X(N4, 1); 
      Y24 = X(N2, 2) - X(N4, 2);
      if (X13 * X13 + Y13 * Y13) > 1.1 * (X24 * X24 + Y24 * Y24) 
         NOC(I1, 3) = N4; 
         NOC(2 * I, 3) = N2;
      end
   end
end

%------------------------  function BlockXY  ---------------------------
function XP = BlockXY(KW, KSW);
global NEN NS NW NSJ NSW NGN NNW NSD NWD
global SR WR XB

%======  Coordinates of 8-Nodes of the Block  ======
N1 = KSW + KW - 1;
XP(1, 1) = XB(N1, 1);
XP(1, 2) = XB(N1, 2);
XP(3, 1) = XB(N1 + 1, 1);
XP(3, 2) = XB(N1 + 1, 2);
XP(5, 1) = XB(N1 + NS + 2, 1);
XP(5, 2) = XB(N1 + NS + 2, 2);
XP(7, 1) = XB(N1 + NS + 1, 1);
XP(7, 2) = XB(N1 + NS + 1, 2);
XP(2, 1) = SR(KSW, 1);
XP(2, 2) = SR(KSW, 2);
XP(6, 1) = SR(KSW + NS, 1);
XP(6, 2) = SR(KSW + NS, 2);
XP(8, 1) = WR(N1, 1);
XP(8, 2) = WR(N1, 2);
XP(4, 1) = WR(N1 + 1, 1);
XP(4, 2) = WR(N1 + 1, 2);
     
%------------------------  function Shape---------------------------
function SH = Shape(XI, ETA)

%==============  Shape Functions  ================
SH(1) = -(1 - XI) * (1 - ETA) * (1 + XI + ETA) / 4;
SH(2) = (1 - XI * XI) * (1 - ETA) / 2;
SH(3) = -(1 + XI) * (1 - ETA) * (1 - XI + ETA) / 4;
SH(4) = (1 - ETA * ETA) * (1 + XI) / 2;
SH(5) = -(1 + XI) * (1 + ETA) * (1 - XI - ETA) / 4;
SH(6) = (1 - XI * XI) * (1 + ETA) / 2;
SH(7) = -(1 - XI) * (1 + ETA) * (1 + XI - ETA) / 4;
SH(8) = (1 - ETA * ETA) * (1 - XI) / 2;

%------------------------  function Output  ---------------------------
function [] = Output();
global NEN NS NW NSJ NSW NGN NNW NSD NWD
global X NOC MAT
global TITLE
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW

%===== Print Displacements, Stresses, and Reactions
%FILE2 = input('Output Data File Name <DOS file name>','s');
FILE2 = input('Give name of MESHGEN output file (with extension) ','s');
LOUT  = fopen(FILE2,'w');

NDIM = 2;
NDN = 2;
ND = 0;
NL = 0;
NMPC = 0;

fprintf(LOUT,'Program MESHGEN - CHANDRUPATLA & BELEGUNDU\n');
fprintf(LOUT,'%s',TITLE);
fprintf(LOUT,'  NN   NE   NM  NDIM  NEN  NDN\n');

TMP =[NN NE NM NDIM NEN NDN];
fprintf(LOUT,'%5d',TMP');
fprintf(LOUT,'\n');

fprintf(LOUT,'   ND   NL   NMPC\n');
TMP = [ND NL NMPC];
fprintf(LOUT,'%5d',TMP');
fprintf(LOUT,'\n');

fprintf(LOUT,' Node#        X        Y\n');
for i=1:NN
   fprintf(LOUT,'%5d %15.5e %15.5e\n',i,X(i,:));
end
if NEN == 3
   fprintf(LOUT,' Elem#    Node1   Node2   Node3 	Mat#\n');
   for i=1:NE
      fprintf(LOUT,'%5d %7d %7d %7d %7d\n ',i,  NOC(i,:), MAT(i));
   end
elseif NEN == 4
	fprintf(LOUT,' Elem#    Node1   Node2	 Node3 	Node4   Mat#\n');
   for i=1:NE
      fprintf(LOUT,'%5d %7d %7d %7d %7d %7d\n ',i, NOC(i,:), MAT(i));
   end
end
fclose(LOUT);