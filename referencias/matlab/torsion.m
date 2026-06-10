function [] = torsion();
clear
close
%------------------------ TORSION2  ---------------------------
disp('==========================================');
disp('         PROGRAM TORSION2                 ');
disp('    TORSION WITH 3-NODED TRIANGLES        ');
disp('   T.R.Chandrupatla and A.D.Belegundu     ');
disp('==========================================');
InputData;
Bandwidth;
Stiffness;
ModifyForBC;
BandSolver;
AngleOfTwist;
Output;

%------------------------  function InputData  ---------------------------
function [] = InputData();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT F S PM NU U
global TORQUE SFAC
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2 IPL

disp(blanks(1));
FILE1 = input('Input Data File Name ','s');
LINP  = fopen(FILE1,'r');
FILE2 = input('Output Data File Name ','s');
LOUT  = fopen(FILE2,'w');

DUMMY = fgets(LINP);
TITLE = fgets(LINP);
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[NN, NE, NM, NDIM, NEN, NDN] = deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5),TMP(6));

DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));

% NPR = 1;
% NMPC = 0;
% NM = 1;
%--- ND = NO. OF SPECIFIED STRESS FUNCTION VALUES
%--- NL = NO. OF GENERALIZED NODAL FORCES "0" HERE
%--- NPR =1 (SHEAR MODULUS) AND NMPC = 0 FOR THIS PROGRAM
%--- ELEMENT CHARACTERISTIC NOT USED
%--- NO. OF MATERIALS = 1 FOR THIS PROGRAM

% Dimensioned for minimum 3 properties
disp('PLOT CHOICE');
disp('  1) No Plot Data');
disp('  2) Create Data File Containing Stress Function Values');
IPL = input('  Choose 1(defalut) or 2 :');
%     --- default is no data
if isempty(IPL) | IPL<1 | IPL>2
   IPL = 1;
end
if IPL > 1 
    disp(blanks(1));
    FILE3 = input('Give Data File Name for Nodal Values of Theta ','s');
    LOUT2  = fopen(FILE3,'w');
end

%----- Coordinates -----
DUMMY = fgets(LINP);
for I=1:NN
   TMP = str2num(fgets(LINP));
   [N, X(N,:)]=deal(TMP(1),TMP(2:1+NDIM));
end
%----- Connectivity -----
DUMMY = fgets(LINP);
for I=1:NE
   TMP = str2num(fgets(LINP));
   [N,NOC(N,:), MAT(N,:)] = ...
      deal(TMP(1),TMP(2:1+NEN),TMP(2+NEN));
end

%----- Specified Displacements -----
DUMMY = fgets(LINP);
for I=1:ND
   TMP = str2num(fgets(LINP));
   [NU(I,:),U(I,:)] = deal(TMP(1), TMP(2));
end
%----- Component Loads -----
DUMMY = fgets(LINP);
F = zeros(NN,1);
for I=1:NL
   TMP = str2num(fgets(LINP));
   [N,F(N)]=deal(TMP(1),TMP(2));
end

%----- Shear Modulus of Material
DUMMY = fgets(LINP);
for I=1:NM
   TMP = str2num(fgets(LINP));
   [N, PM(N,:)] = deal(TMP(1), TMP(2));
end
TORQUE = input('TORQUE : ');
SFAC = input('SYMMETRY FACTOR (eg. if 1/4 symmetry, then=4.0)= ');

fclose(LINP);

%------------------------  function Bandwidth  ---------------------------
function []=Bandwidth();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT F S PM NU U
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2 IPL
%----- Bandwidth NBW from Connectivity NOC() and MPC
NBW = 0;
for I = 1:NE
   NMIN = NOC(I, 1);
   NMAX = NOC(I, 1);
   for J = 2:3
      if NMIN > NOC(I, J); NMIN = NOC(I, J); end
      if NMAX < NOC(I, J); NMAX = NOC(I, J); end
   end
   NTMP = NDN * (NMAX - NMIN + 1);
   if NBW < NTMP; NBW = NTMP; end
end
for I = 1:NMPC
   NABS = abs(MPC(I, 1) - MPC(I, 2)) + 1;
   if (NBW < NABS); NBW = NABS; end
end
disp(sprintf('Bandwidth = %d', NBW));

%------------------------  function Stiffness  ---------------------------
function []=Stiffness();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT F S PM NU U
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2 IPL
%----- Global Stiffness Matrix
S = zeros(NN,NBW);

for I = 1:NE
   I1 = NOC(I, 1);
   I2 = NOC(I, 2);
   I3 = NOC(I, 3);
   X32 = X(I3, 1) - X(I2, 1);
   X13 = X(I1, 1) - X(I3, 1);
   X21 = X(I2, 1) - X(I1, 1);
   Y23 = X(I2, 2) - X(I3, 2);
   Y31 = X(I3, 2) - X(I1, 2);
   Y12 = X(I1, 2) - X(I2, 2);
   DETJ = X13 * Y23 - X32 * Y31;
   AREA = .5 * abs(DETJ);
%  --- LOAD CALCULATION
   C = 2 * AREA / 3;
   F(I1) = F(I1) + C;
   F(I2) = F(I2) + C;
   F(I3) = F(I3) + C;
%  --- STIFFNESS FORMATION
   BT(1, 1) = Y23;
   BT(1, 2) = Y31;
   BT(1, 3) = Y12;
   BT(2, 1) = X32;
   BT(2, 2) = X13;
   BT(2, 3) = X21;
	for II = 1:3
	   for JJ = 1:2
         BT(JJ, II) = BT(JJ, II) / DETJ;
      end
   end
   for II = 1:3
	   for JJ = 1:3
         II1 = NOC(I, II);
         II2 = NOC(I, JJ);
         if (II1 <= II2)
            SUM = 0.;
		      for J = 1:2
               SUM = SUM + BT(J, II) * BT(J, JJ);
            end  
            IC = II2 - II1 + 1;
            S(II1, IC) = S(II1, IC) + SUM * AREA;
         end
      end
   end
end

%------------------------  function ModifyForBC  ---------------------------
function []=ModifyForBC();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT F S PM NU U
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2 IPL
%----- Decide Penalty Parameter CNST -----
CNST = 0;
for I = 1:NN
   if CNST < S(I, 1); CNST = S(I, 1); end
end
CNST = CNST * 1000000;

%----- Modify for Boundary Conditions -----
%    --- Displacement BC ---
for I = 1:ND
   N = NU(I);
   S(N, 1) = S(N, 1) + CNST;
   F(N) = F(N) + CNST * U(I);
end

%------------------------  function BandSolver  ---------------------------
function []=BandSolver();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT F S PM NU U
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2 IPL
%----- Equation Solving using Band Solver -----
disp('Solving using Band Solver(bansol.m)');
[F] = bansol(NN,NBW,S,F);

%------------------------  function AngleOfTwist  ---------------------------
function [] = AngleOfTwist();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT F S PM NU U
global TORQUE SFAC TAU ALPHA
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2 IPL
%----- ANGLE OF TWIST PER UNIT LENGTH
SUM = 0. ;
for I = 1:NE
   I1 = NOC(I, 1);
   I2 = NOC(I, 2);
   I3 = NOC(I, 3);
   X32 = X(I3, 1) - X(I2, 1);
   X13 = X(I1, 1) - X(I3, 1);
   X21 = X(I2, 1) - X(I1, 1);
   Y23 = X(I2, 2) - X(I3, 2);
   Y31 = X(I3, 2) - X(I1, 2);
   Y12 = X(I1, 2) - X(I2, 2);
   DETJ = X13 * Y23 - X32 * Y31;
   SUM = SUM + abs(DETJ) / 3 * (F(I1) + F(I2) + F(I3));
end
SMOD = PM(1, 1);
ALPHA = TORQUE / SMOD / SUM / SFAC;
for I = 1:NE
   I1 = NOC(I, 1);
   I2 = NOC(I, 2);
   I3 = NOC(I, 3);
   X32 = X(I3, 1) - X(I2, 1);
   X13 = X(I1, 1) - X(I3, 1);
   X21 = X(I2, 1) - X(I1, 1);
   Y23 = X(I2, 2) - X(I3, 2);
   Y31 = X(I3, 2) - X(I1, 2);
   Y12 = X(I1, 2) - X(I2, 2);
   DETJ = X13 * Y23 - X32 * Y31;
   BT(1, 1) = Y23;
   BT(1, 2) = Y31;
   BT(1, 3) = Y12;
   BT(2, 1) = X32;
   BT(2, 2) = X13;
   BT(2, 3) = X21;
   
	for II = 1:3
	   for JJ = 1:2
         BT(JJ, II) = BT(JJ, II) / DETJ;
      end
   end
   
	TAUYZ = -(BT(1, 1)*F(I1) + BT(1, 2)*F(I2) + BT(1, 3)*F(I3));
	TAUXZ = BT(2, 1)*F(I1) + BT(2, 2)*F(I2) + BT(2, 3)*F(I3);
   TAUYZ = TAUYZ * SMOD * ALPHA;
   TAUXZ = TAUXZ * SMOD * ALPHA;
   TAU(I, 1) = TAUYZ;
   TAU(I, 2) = TAUXZ;
end


%------------------------  function Output  ---------------------------
function [] = Output();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT F S PM NU U
global TORQUE SFAC TAU ALPHA
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2 IPL

%----------------------------------------------
disp(sprintf('Output for Input Data from file %s\n',FILE1));
fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);

disp(' NODE#   Stress Function Value');
fprintf(LOUT,' NODE#   Stress Function Value\n');
I=[1:NN]';
% print a matrix
disp(sprintf(' %4d %15.4E\n',[I,F(I)]'));
fprintf(LOUT,' %4d %15.4E\n',[I,F(I)]');


if IPL > 1
   fprintf(LOUT2,'Stress Function Value\n');
   fprintf(LOUT2,'%15.4E\n',[F(I)]');
   disp(sprintf('Stress Function Values in Data file %s ',FILE3));
   disp('Run CONTOUR1 or CONTOUR2');
   fclose(LOUT2);
end

%----- ANGLE OF TWIST PER UNIT LENGTH
disp(sprintf('TWIST PER UNIT LENGTH = %15.5E', ALPHA));
fprintf(LOUT,'TWIST PER UNIT LENGTH = %15.5E\n', ALPHA);
fprintf(LOUT,'-- SHEARING STRESSES TAUYZ, TAUXZ IN EACH ELEMENT\n');
fprintf(LOUT,'ELEMENT#    TAUYZ       TAUXZ\n');
for I = 1:NE
   disp(sprintf('%d %16.5E %16.5E',I, TAU(I,1), TAU(I,2)));
   fprintf(LOUT,'%d %16.5E %16.5E\n',I, TAU(I,1), TAU(I,2));
end
fclose(LOUT);
 
disp(sprintf('The Results are available in the text file %s', FILE2));
disp('View using a text processor');


