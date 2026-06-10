function [] = truss2n();
close all
clear all

disp('=======================================');
disp('         PROGRAM TRUSS                 ');
disp('   Two Dimensional Truss Analysis      ');
disp('   T.R.Chandrupatla and A.D.Belegundu  ');
disp('=======================================');

InputData;
Bandwidth;
Stiffness;
ModifyForBC;
BandSolver;
StressCalc;
ReactionCalc;
Output;

%------------------------  function InputData  ---------------------------
function [] = InputData();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

FILE1 = input('Input Data File Name ','s');
LINP  = fopen(FILE1,'r');
FILE2 = input('Output Data File Name ','s');
LOUT  = fopen(FILE2,'w');

DUMMY = fgets(LINP);
TITLE = fgets(LINP);
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[NN, NE, NM, NDIM, NEN, NDN] = deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5),TMP(6));

NQ = NDN * NN;

DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));

NPR=2; % E, ALPHA

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
   [N,NOC(N,:), MAT(N,:), AREA(N,:), DT(N,:)] = ...
      deal(TMP(1),TMP(2:1+NEN), TMP(2+NEN), TMP(3+NEN), TMP(4+NEN));
end

%----- Specified Displacements -----
DUMMY = fgets(LINP);
for I=1:ND
   TMP = str2num(fgets(LINP));
   [NU(I,:),U(I,:)] = deal(TMP(1), TMP(2));
end
%----- Component Loads -----
DUMMY = fgets(LINP);
F = zeros(NQ,1);
for I=1:NL
   TMP = str2num(fgets(LINP));
   [N,F(N)]=deal(TMP(1),TMP(2));
end

%----- Material Properties -----
DUMMY = fgets(LINP);
for I=1:NM
   TMP = str2num(fgets(LINP));
   [N, PM(N,:)] = deal(TMP(1), TMP(2:NPR+1));
end
%----- Multi-point Constraints B1*Qi+B2*Qj=B0
if NMPC > 0
   DUMMY = fgets(LINP);
   for I=1:NMPC
   	TMP = str2num(fgets(LINP));
      [BT(I,1), MPC(I,1), BT(I,2), MPC(I,2), BT(I,3)] = ...
         			deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5));
   end
end
fclose(LINP);

%------------------------  function Bandwidth  ---------------------------
function []=Bandwidth();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- Bandwidth Evaluation -----
NBW = 0;
for N=1:NE
   NABS = NDN*(abs(NOC(N, 1) - NOC(N, 2)) + 1);
   if (NBW < NABS)
      NBW = NABS;
   end
end
for I=1:NMPC
   NABS = abs(MPC(I, 1) - MPC(I, 2)) + 1;
   if (NBW < NABS)
      NBW = NABS;
   end
end
disp(sprintf('Bandwidth = %d', NBW));

%------------------------  function Stiffness  ---------------------------
function []=Stiffness();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- First initialize stiffness matrix
S = zeros(NQ,NBW);      

%-----  Stiffness Matrix -----
for N = 1:NE
   I1 = NOC(N, 1);
   I2 = NOC(N, 2);
   I3 = MAT(N);
   X21 = X(I2, 1) - X(I1, 1);
   Y21 = X(I2, 2) - X(I1, 2);
   EL = sqrt(X21 * X21 + Y21 * Y21);
   EAL = PM(I3, 1) * AREA(N) / EL;
   CS = X21 / EL;
   SN = Y21 / EL;
%----------- Element Stiffness Matrix SE() -----------
   SE(1, 1) = CS * CS * EAL;
   SE(1, 2) = CS * SN * EAL;
   SE(2, 1) = SE(1, 2);
   SE(1, 3) = -CS * CS * EAL;
   SE(3, 1) = SE(1, 3);
   SE(1, 4) = -CS * SN * EAL;
   SE(4, 1) = SE(1, 4);
   SE(2, 2) = SN * SN * EAL;
   SE(2, 3) = -CS * SN * EAL;
   SE(3, 2) = SE(2, 3);
   SE(2, 4) = -SN * SN * EAL;
   SE(4, 2) = SE(2, 4);
   SE(3, 3) = CS * CS * EAL;
   SE(3, 4) = CS * SN * EAL;
   SE(4, 3) = SE(3, 4);
	SE(4, 4) = SN * SN * EAL;
%-------------- Temperature Load TL() ---------------
   EE0 = PM(I3, 2) * DT(N) * PM(I3, 1) * AREA(N);
   TL(1) = -EE0 * CS;
   TL(2) = -EE0 * SN;
   TL(3) = EE0 * CS;
   TL(4) = EE0 * SN;
   disp(sprintf('..... Adding %dth Element Stiffness to Global Locations',N));
	for II = 1:NEN
      NRT = NDN * (NOC(N, II) - 1);
	   for IT = 1:NDN
         NR = NRT + IT;
         I = NDN * (II - 1) + IT;
         for JJ = 1:NEN
            NCT = NDN * (NOC(N, JJ) - 1);
            for JT = 1:NDN
               J = NDN * (JJ - 1) + JT;
               NC = NCT + JT - NR + 1;
               if NC > 0
                  S(NR, NC) = S(NR, NC) + SE(I, J);
               end
            end
         end
         F(NR) = F(NR) + TL(I);
      end
   end
end

%------------------------  function ModifyForBC  ---------------------------
function []=ModifyForBC();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- Decide Penalty Parameter CNST -----
CNST = 0;
for I = 1:NQ
   if CNST < S(I, 1)
      CNST = S(I, 1);
   end
end
CNST = CNST * 10000;

%----- Modify for Boundary Conditions -----
%    --- Displacement BC ---
for I = 1:ND
   N = NU(I);
   S(N, 1) = S(N, 1) + CNST;
   F(N) = F(N) + CNST * U(I);
end
%--- Multi-point Constraints ---
for I = 1:NMPC
   I1 = MPC(I, 1);
   I2 = MPC(I, 2);
   S(I1, 1) = S(I1, 1) + CNST * BT(I, 1) * BT(I, 1);
   S(I2, 1) = S(I2, 1) + CNST * BT(I, 2) * BT(I, 2);
   IR = I1;
   if IR > I2
      IR = I2;
   end
   IC = ABS(I2 - I1) + 1;
   S(IR, IC) = S(IR, IC) + CNST * BT(I, 1) * BT(I, 2);
   F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3);
   F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3);
end

%------------------------  function BandSolver  ---------------------------
function []=BandSolver();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- Equation Solving using Band Solver -----
disp('Solving using Band Solver(bansol.m)');
[F] = bansol(NQ,NBW,S,F);

%------------------------  function StressCalc  ---------------------------
function []=StressCalc();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- Stress Calculation -----

for I = 1:NE
   I1 = NOC(I, 1);
   I2 = NOC(I, 2);
   I3 = MAT(I);
   X21 = X(I2, 1) - X(I1, 1);
   Y21 = X(I2, 2) - X(I1, 2);
   EL = sqrt(X21 * X21 + Y21 * Y21);
   CS = X21 / EL;
   SN = Y21 / EL;
   J2 = 2 * I1;
   J1 = J2 - 1;
   K2 = 2 * I2;
   K1 = K2 - 1;
   DLT = (F(K1) - F(J1)) * CS + (F(K2) - F(J2)) * SN;
   STRESS(I) = PM(I3, 1) * (DLT / EL - PM(I3, 2) * DT(I));
end

%------------------------  function ReactionCalc  ---------------------------
function []=ReactionCalc();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- Reaction Calculation -----
for I = 1:ND
   N = NU(I);
   REACT(I) = CNST * (U(I) - F(N));
end

%------------------------  function Output  ---------------------------
function []=Output();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

disp(sprintf('Output for Input Data from file %s\n',FILE1));
fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);

%----- Displacement  -----
disp(' Node#    X-Displ         Y-Displ');
fprintf(LOUT,' Node#    X-Displ         Y-Displ\n');
for I=1:NN
   % print a matrix
	disp(sprintf(' %4d %15.4E %15.4E',[I,F(2*I-1),F(2*I)]'));
	fprintf(LOUT,' %4d %15.4E %15.4E\n',[I,F(2*I-1),F(2*I)]');
end
%----- Stress  -----
disp(' Elem#    Stress');
fprintf(LOUT,' Elem#    Stress\n');
for I = 1:NE
   disp(sprintf(' %4d %15.4E',I,STRESS(I)));
   fprintf(LOUT,' %4d %15.4E\n',I,STRESS(I));
end

%----- Reaction  -----
disp(sprintf('  DOF#     Reaction'));
fprintf(LOUT,'  DOF#     Reaction\n');
for I = 1:ND
   N = NU(I);
   disp(sprintf(' %4d %15.4E',N,REACT(I)));
   fprintf(LOUT,' %4d %15.4E\n',N,REACT(I));
end

fclose(LOUT);
disp(sprintf('Results are now available in the text file %s', FILE2));

