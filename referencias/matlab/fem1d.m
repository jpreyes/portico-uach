function [] = fem1d2n();
clear all;close all;
% main 

%-----------------------  FEM1D2  ------------------------
disp('***************************************');
disp('*          PROGRAM FEM1D2             *');
disp('*    WITH MULTI-POINT CONSTRAINTS     *');
disp('* T.R.Chandrupatla and A.D.Belegundu  *');
disp('***************************************');      

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
global ND NL NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NPR

FILE1 = input('Input Data File Name: ','s');
LINP  = fopen(FILE1,'r');
FILE2 = input('Output Data File Name: ','s');
LOUT  = fopen(FILE2,'w');
disp(blanks(1));
DUMMY = fgets(LINP);
TITLE = fgets(LINP);
DUMMY = fgets(LINP);

TMP = str2num(fgets(LINP));
[NN, NE, NM, NDIM, NEN, NDN] = deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5),TMP(6));

DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));

NPR = 2;  %E, ALPHA

%----- Coordinates -----
DUMMY = fgets(LINP);
for I=1:NN
   TMP = str2num(fgets(LINP));
   [N, X(N,:)]=deal(TMP(1),TMP(2:1+NDIM));
end
%----- Connectivity -----
DUMMY = fgets(LINP);
%TMP = fscanf(LINP,'%g\n',[6 NE]);	% 6 columns, NE rows
%TMP = TMP';
%[NOC, MAT, AREA, DT] = deal(TMP(:,2:3), TMP(:,4), TMP(:,5), TMP(:,6));
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
F = zeros(NN,1);
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
global ND NL NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT

%----- Bandwidth Evaluation -----
NBW = 0;
for N=1:NE
   NABS = floor(abs(NOC(N, 1) - NOC(N, 2))) + 1;
   if (NBW < NABS)
      NBW = NABS;
   end
end
for I=1:NMPC
   NABS = floor(abs(MPC(I, 1) - MPC(I, 2))) + 1;
   if (NBW < NABS)
      NBW = NABS;
   end
end

%------------------------  function Stiffness  ---------------------------
function []=Stiffness();
global NN NE NM NDIM NEN NDN
global ND NL NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NPR

%----- First initialize stiffness matrix
S = zeros(NN,NBW);      
%----- Stiffness Matrix -----
for N=1:NE
   N1 = NOC(N, 1);
   N2 = NOC(N, 2);
   N3 = MAT(N);
   X21 = X(N2) - X(N1);
   EL = abs(X21);
   EAL = PM(N3, 1) * AREA(N) / EL;
   if (NPR > 1)
      C = PM(N3, 2);
   else
   	C = 1;   
   end
   TL = PM(N3, 1) * C * DT(N) * AREA(N) * EL / X21;
%----- Temperature Loads -----
   F(N1) = F(N1) - TL;
   F(N2) = F(N2) + TL;
%----- Element Stiffness in Global Locations -----
   S(N1, 1) = S(N1, 1) + EAL;
   S(N2, 1) = S(N2, 1) + EAL;
   IR = N1;
   if (IR > N2)
   	IR = N2;
   end
   IC = floor(abs(N2 - N1)) + 1;
   S(IR, IC) = S(IR, IC) - EAL;
end

%------------------------  function ModifyForBC  ---------------------------
function []=ModifyForBC();
global NN NE NM NDIM NEN NDN
global ND NL NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT

%----- Decide Penalty Parameter CNST -----
CNST = 0;
for I=1:NN
   if (CNST < S(I, 1))
      CNST = S(I, 1);
   end
end
CNST = CNST * 10000;
%----- Modify for Boundary Conditions -----
%--- Displacement BC ---
for I=1:ND
   N = NU(I);
   S(N, 1) = S(N, 1) + CNST;
   F(N) = F(N) + CNST * U(I);
end
%--- Multi-point Constraints ---
for I=1:NMPC
   I1 = MPC(I, 1);
   I2 = MPC(I, 2);
   S(I1, 1) = S(I1, 1) + CNST * BT(I, 1) * BT(I, 1);
   S(I2, 1) = S(I2, 1) + CNST * BT(I, 2) * BT(I, 2);
   IR = I1;
   if (IR > I2)
      IR = I2;;
   end
   IC = floor(abs(I2 - I1)) + 1;
   S(IR, IC) = S(IR, IC) + CNST * BT(I, 1) * BT(I, 2);
   F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3);
   F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3);
end

%------------------------  function BandSolver  ---------------------------
function []=BandSolver();
global NN NE NM NDIM NEN NDN
global ND NL NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
disp('Solving using Band Solver(bansol.m)');
[F] = bansol(NN,NBW,S,F);

%------------------------  function StressCalc  ---------------------------
function []=StressCalc();
global NN NE NM NDIM NEN NDN
global ND NL NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NPR
%----- Stress Calculation -----
for N=1:NE
   N1 = NOC(N, 1);
	N2 = NOC(N, 2);
   N3 = MAT(N);
   EPS = (F(N2) - F(N1)) / (X(N2) - X(N1));
   if (NPR > 1)
      C = PM(N3, 2);
   end
   STRESS(N) = PM(N3, 1) * (EPS - C * DT(N));
end

%------------------------  function ReactionCalc  ---------------------------
function []=ReactionCalc();
global NN NE NM NDIM NEN NDN
global ND NL NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
%----- Reaction Calculation -----
for I=1:ND
   N = NU(I);
   REACT(I) = CNST * (U(I) - F(N));
end

%------------------------  function Output  ---------------------------
function []=Output();
global NN NE NM NDIM NEN NDN
global ND NL NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);
disp('NODE#   DISPLACEMENT');
fprintf(LOUT,'NODE#   DISPLACEMENT\n');
for I=1:NN
	disp(sprintf('%d %15.5G', I, F(I)));
	fprintf(LOUT,'%d %15.5G\n', I, F(I));
end
%----- Stress  -----
disp('ELEM#    STRESS');
fprintf(LOUT,'ELEM#    STRESS\n');
for N=1:NE
	disp(sprintf('%d %15.5G',N, STRESS(N)));
	fprintf(LOUT,'%d %15.5G\n',N, STRESS(N));
end
%----- Reaction  -----
disp('NODE#    REACTION');
fprintf(LOUT,'NODE#    REACTION\n');
for I=1:ND
   N = NU(I);
	disp(sprintf('%d %15.5G',N, REACT(I)));
	fprintf(LOUT,'%d %15.5G\n', N, REACT(I));
end

fclose(LOUT);
disp(sprintf('RESULTS ARE IN FILE : %s', FILE2));


