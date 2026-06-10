function []=frame2dn()
clear all
close all
%------------------------ FRAME2D  ---------------------------
disp('========================================');
disp('            PROGRAM FRAME2D              ');
disp('        2-D Frame Analysis              ');
disp('   T.R.Chandrupatla and A.D.Belegundu   ');
disp('========================================');

InputData;
Bandwidth;
Stiffness;
AddLoads;
ModifyForBC;
BandSolver;
EndActions;
ReactionCalc;
Output;

%------------------------  function InputData  ---------------------------
function [] = InputData();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT ARIN UDL S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT LOUT2
global NQ


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

NQ = NDN * NN;

DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));

NPR=1;  % E

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
   [N,NOC(N,:), MAT(N), ARIN(N,:),UDL(N)] = ...
      deal(TMP(1),TMP(2:1+NEN), TMP(2+NEN), TMP(3+NEN:4+NEN), TMP(5+NEN));
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
global X NOC F AREA MAT ARIN UDL S
global PM NU U MPC BT STRESS REACT

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
global X NOC F AREA MAT ARIN UDL S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- Global Stiffness Matrix
S = zeros(NQ,NBW);

for N = 1:NE
   disp(sprintf('Forming Stiffness Matrix of Element %d', N));
   
%--------  Element Stiffness  -----
   ISTF = 2;
   [ALMBDA,SE,SEP,EL]=elstiff(N,NOC,X,MAT,PM,ARIN,ISTF);
      
   disp('.... Placing in Global Locations');
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
               if (NC > 0)
                  S(NR, NC) = S(NR, NC) + SE(I, J);
               end
            end
         end
      end
   end
end

%------------------------  function AddLoads  ---------------------------
function []=AddLoads();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT ARIN UDL S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ
%----- Loads due to uniformly distributed load on element
for N = 1:NE
   if (abs(UDL(N)) > 0)
      ISTF = 1;
      [ALMBDA,SE,SEP,EL]=elstiff(N,NOC,X,MAT,PM,ARIN,ISTF);
      I1 = NOC(N, 1);
      I2 = NOC(N, 2);
      ED(1) = 0;
      ED(4) = 0;
      ED(2) = UDL(N) * EL / 2;
      ED(5) = ED(2);
      ED(3) = UDL(N) * EL ^ 2 / 12;
      ED(6) = -ED(3);
      
      EDP = ALMBDA' * ED';
      
      for I = 1:3
         F(3 * I1 - 3 + I) = F(3 * I1 - 3 + I) + EDP(I);
         F(3 * I2 - 3 + I) = F(3 * I2 - 3 + I) + EDP(I + 3);
      end
   end
end

%------------------------  function ModifyForBC  ---------------------------
function []=ModifyForBC();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT ARIN UDL S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- Decide Penalty Parameter CNST -----
CNST = 0;
for I = 1:NQ
   if CNST < S(I, 1); CNST = S(I, 1); end
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
   if IR > I2; IR = I2; end
   IC = abs(I2 - I1) + 1;
   S(IR, IC) = S(IR, IC) + CNST * BT(I, 1) * BT(I, 2);
   F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3);
   F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3);
end

%------------------------  function BandSolver  ---------------------------
function []=BandSolver();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT ARIN UDL S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- Equation Solving using Band Solver -----
disp('Solving using Band Solver(bansol.m)');
[F] = bansol(NQ,NBW,S,F);

%------------------------  function EndActions  ---------------------------
function []=EndActions();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT ARIN UDL S
global PM NU U MPC BT EF REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ
%---- Member End-Actions
for N = 1:NE
   ISTF = 1;
   [ALMBDA,SE,SEP,EL]=elstiff(N,NOC,X,MAT,PM,ARIN,ISTF);
   I1 = NOC(N, 1);
   I2 = NOC(N, 2);
   for I = 1:3
      ED(I) = F(3 * I1 - 3 + I);
      ED(I + 3) = F(3 * I2 - 3 + I);
   end       
   EDP = ALMBDA * ED';
%  END FORCES DUE TO DISTRIBUTED LOADS
   if (abs(UDL(N)) > 0)
      ED(1) = 0;
      ED(4) = 0;
      ED(2) = -UDL(N) * EL / 2;
      ED(5) = ED(2);
      ED(3) = -UDL(N) * EL ^ 2 / 12;
      ED(6) = -ED(3);
   else
      ED = zeros(1,6);
   end
   for I = 1:6
      EF(N,I) = ED(I);
      for K = 1:6
         EF(N,I) = EF(N,I) + SEP(I, K) * EDP(K);
      end
   end
end

%------------------------  function ReactionCalc  ---------------------------
function []=ReactionCalc();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT TH DT S
global PM NU U MPC BT REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
%----- Reaction Calculation -----
for I = 1:ND
   N = NU(I);
   REACT(I) = CNST * (U(I) - F(N));
end

%------------------------  function Output  ---------------------------
function []=Output();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT TH DT S
global PM NU U MPC BT EF REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT LOUT2

disp(sprintf('Output for Input Data from file %s\n',FILE1));
fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);

disp('NODE#     X-Displ.        Y-Displ.        Z-Rot.');
fprintf(LOUT,'NODE#     X-Displ.        Y-Displ.        Z-Rot.\n');
for I=1:NN
   II = (I - 1) * 3;
	% print a matrix
	disp(sprintf(' %4d %15.4E %15.4E %15.4E',[I,F(II+1),F(II+2),F(II+3)]'));
	fprintf(LOUT,' %4d %15.4E %15.4E %15.4E\n',[I,F(II+1),F(II+2),F(II+3)]');
end
%----- Reaction Calculation -----
disp(sprintf('  DOF#     Reaction'));
fprintf(LOUT,'  DOF#     Reaction\n');
for I = 1:ND
   N = NU(I);
   disp(sprintf(' %4d %15.4E',N,REACT(I)));
   fprintf(LOUT,' %4d %15.4E\n',N,REACT(I));
end

%---- Member End-Actions
disp(sprintf(' Member End-Forces'));
fprintf(LOUT,' Member End-Forces\n');
for N = 1:NE
   disp(sprintf(' Member # %d', N));
   fprintf(LOUT,' Member # %d\n', N);
   for I = 1:2
      II = (I - 1) * 3;
      disp(sprintf('%15.4E %15.4E %15.4E', EF(N,II+1),EF(N,II+2),EF(N,II+3)));
      fprintf(LOUT,'%15.4E %15.4E %15.4E\n', EF(N,II+1),EF(N,II+2),EF(N,II+3));
   end
end

fclose(LOUT);
disp(sprintf('The Results are available in the text file %s', FILE2));

%------------------------ function elstiff ---------------------------
function [ALMBDA,SE,SEP,EL] = elstiff(N,NOC,X,MAT,PM,ARIN,ISTF)
%===== SUBROUTINE ELEMENT STIFFNESS =====
%----- Element Stiffness Matrix -----
I1 = NOC(N, 1);
I2 = NOC(N, 2);
M = MAT(N);
X21 = X(I2, 1) - X(I1, 1);
Y21 = X(I2, 2) - X(I1, 2);
EL = sqrt(X21 * X21 + Y21 * Y21);
EAL = PM(M, 1) * ARIN(N, 1) / EL;
EIZL = PM(M, 1) * ARIN(N, 2) / EL;
      
SEP = zeros(6);
SEP(1, 1) = EAL;
SEP(1, 4) = -EAL;
SEP(4, 4) = EAL;
SEP(2, 2) = 12 * EIZL / EL ^ 2;
SEP(2, 3) = 6 * EIZL / EL;
SEP(2, 5) = -SEP(2, 2);
SEP(2, 6) = SEP(2, 3);
SEP(3, 3) = 4 * EIZL;
SEP(3, 5) = -6 * EIZL / EL;
SEP(3, 6) = 2 * EIZL;
SEP(5, 5) = 12 * EIZL / EL ^ 2;
SEP(5, 6) = -6 * EIZL / EL;
SEP(6, 6) = 4 * EIZL;
for I = 1:6
   for J = I:6
      SEP(J, I) = SEP(I, J);
   end
end
%   CONVERT ELEMENT STIFFNESS MATRIX TO GLOBAL SYSTEM
DCOS(1, 1) = X21 / EL;
DCOS(1, 2) = Y21 / EL;
DCOS(1, 3) = 0;
DCOS(2, 1) = -DCOS(1, 2);
DCOS(2, 2) = DCOS(1, 1);
DCOS(2, 3) = 0;
DCOS(3, 1) = 0;
DCOS(3, 2) = 0;
DCOS(3, 3) = 1;
      
ALMBDA = zeros(6);
             
for K = 1:2
   IK = 3 * (K - 1);
   for I = 1:3
      for J = 1:3
         ALMBDA(I + IK, J + IK) = DCOS(I, J);
      end
   end
end

if ISTF ~=1
    SE = SEP * ALMBDA;
    SEP = SE;
    SE = ALMBDA' * SEP;
else
    SE = zeros(6);  % dummy
end