function [] = frame3dn()
clear all
close all
%------------------------ FRAME3D  ---------------------------
disp('========================================');
disp('            PROGRAM FRAME3D              ');
disp('        3-D Frame Analysis              ');
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
[NN, NE, NM, NDIM, NEN, NDN, NNREF] = ...
   deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5),TMP(6),TMP(7));

NNT = NN + NNREF;
NQ = NDN * NN;

DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));

NPR=2;  %  E, G

%----- Coordinates -----
DUMMY = fgets(LINP);
for I=1:NNT
   TMP = str2num(fgets(LINP));
   [N, X(N,:)]=deal(TMP(1),TMP(2:1+NDIM));
end
%----- Connectivity -----
DUMMY = fgets(LINP);
for I=1:NE
   TMP = str2num(fgets(LINP));
   [N,NOC(N,:), MAT(N,:), ARIN(N,:),UDL(N,:)] = ...
      deal(TMP(1),TMP(2:4), TMP(5), ...
           TMP(6:9), TMP(10:11));
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
   if (abs(UDL(N,1)) > 0)|(abs(UDL(N,2)) > 0)
      ISTF = 1;
      [ALMBDA,SE,SEP,EL]=elstiff(N,NOC,X,MAT,PM,ARIN,ISTF);
      I1 = NOC(N, 1);
      I2 = NOC(N, 2);
      ED(1) = 0;
      ED(4) = 0;
      ED(7) = 0;
      ED(10) = 0;
      ED(2) = UDL(N, 1) * EL / 2;
      ED(8) = ED(2);
      ED(6) = UDL(N, 1) * EL^2 / 12;
      ED(12) = -ED(6);
      ED(3) = UDL(N, 2) * EL / 2;
      ED(9) = ED(3);
      ED(5) = -UDL(N, 2) * EL^2 / 12;
      ED(11) = -ED(5);
      
      EDP = ALMBDA' * ED';
      
      for I = 1:6
         F(6 * I1 - 6 + I) = F(6 * I1 - 6 + I) + EDP(I);
         F(6 * I2 - 6 + I) = F(6 * I2 - 6 + I) + EDP(I + 6);
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
   for I = 1:6
      ED(I) = F(6 * I1 - 6 + I);
      ED(I + 6) = F(6 * I2 - 6 + I);
   end       
   EDP = ALMBDA * ED';
%  END FORCES DUE TO DISTRIBUTED LOADS
   if (abs(UDL(N,1)) > 0)|(abs(UDL(N,2)) > 0)
      ED(1) = 0;
      ED(4) = 0;
      ED(7) = 0;
      ED(10) = 0;
      ED(2) = UDL(N, 1) * EL / 2;
      ED(8) = ED(2);
      ED(6) = UDL(N, 1) * EL^2 / 12;
      ED(12) = -ED(6);
      ED(3) = UDL(N, 2) * EL / 2;
      ED(9) = ED(3);
      ED(5) = -UDL(N, 2) * EL^2 / 12;
      ED(11) = -ED(5);
   else
      ED = zeros(1,12);
   end
   for I = 1:12
      EF(N,I) = ED(I);
      for K = 1:12
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

%----------------------------------------------
disp(sprintf('Output for Input Data from file %s\n',FILE1));
fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);

I=[1:NN]';
I1 = 3 * I - 2;
I2 = I1 + 1;
I3 = I1 + 2;
I4 = I1 + 3;
I5 = I1 + 4;
I6 = I1 + 5;  
disp('NODE#     X-Displ.        Y-Displ.        Z-Displ.');
fprintf(LOUT,'NODE#     X-Displ.        Y-Displ.        Z-Rot.\n');
% print a matrix
disp(sprintf(' %4d %15.4E %15.4E %15.4E\n',[I,F(I1),F(I2),F(I3)]'));
fprintf(LOUT,' %4d %15.4E %15.4E %15.4E\n',[I,F(I1),F(I2),F(I3)]');

disp('NODE#     X-Rot.          Y-Rot.          Z-Rot.');
fprintf(LOUT,'NODE#     X-Displ.        Y-Displ.        Z-Rot.\n');
% print a matrix
disp(sprintf(' %4d %15.4E %15.4E %15.4E\n',[I,F(I4),F(I5),F(I6)]'));
fprintf(LOUT,' %4d %15.4E %15.4E %15.4E\n',[I,F(I4),F(I5),F(I6)]');

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
      II = (I - 1) * 6;
      disp(sprintf('%15.4E %15.4E %15.4E %15.4E %15.4E %15.4E', ...
         EF(N,II+1), EF(N,II+2), EF(N,II+3), EF(N,II+4), EF(N,II+5), EF(N,II+6)));
      fprintf(LOUT,'%15.4E %15.4E %15.4E %15.4E %15.4E %15.4E\n', ...
         EF(N,II+1), EF(N,II+2), EF(N,II+3), EF(N,II+4), EF(N,II+5), EF(N,II+6));
   end
end

fclose(LOUT);
disp(sprintf('The Results are available in the text file %s', FILE2));

%------------------------ elstiff ---------------------------
function [ALMBDA,SE,SEP,EL] = elstiff(N,NOC,X,MAT,PM,ARIN,ISTF)
%===== SUBROUTINE ELEMENT STIFFNESS =====
%----- Element Stiffness Matrix -----
I1 = NOC(N, 1);
I2 = NOC(N, 2);
I3 = NOC(N, 3);

M = MAT(N);
X21 = X(I2, 1) - X(I1, 1);
Y21 = X(I2, 2) - X(I1, 2);
Z21 = X(I2, 3) - X(I1, 3);
EL = sqrt(X21 * X21 + Y21 * Y21 + Z21 * Z21);
EAL = PM(M, 1) * ARIN(N, 1) / EL;
EIYL = PM(M, 1) * ARIN(N, 2) / EL;
EIZL = PM(M, 1) * ARIN(N, 3) / EL;
GJL = PM(M, 2) * ARIN(N, 4) / EL;
      
SEP = zeros(12);
SEP(1, 1) = EAL;
SEP(1, 7) = -EAL;
SEP(7, 7) = EAL;
SEP(4, 4) = GJL;
SEP(4, 10) = -GJL;
SEP(10, 10) = GJL;
SEP(2, 2) = 12 * EIZL / EL^2;
SEP(2, 6) = 6 * EIZL / EL;
SEP(2, 8) = -SEP(2, 2);
SEP(2, 12) = SEP(2, 6);
SEP(3, 3) = 12 * EIYL / EL^2;
SEP(3, 5) = -6 * EIYL / EL;
SEP(3, 9) = -SEP(3, 3);
SEP(3, 11) = SEP(3, 5);
SEP(5, 5) = 4 * EIYL;
SEP(5, 9) = 6 * EIYL / EL;
SEP(5, 11) = 2 * EIYL;
SEP(6, 6) = 4 * EIZL;
SEP(6, 8) = -6 * EIZL / EL;
SEP(6, 12) = 2 * EIZL;
SEP(8, 8) = 12 * EIZL / EL^2;
SEP(8, 12) = -6 * EIZL / EL;
SEP(9, 9) = 12 * EIYL / EL^2;
SEP(9, 11) = 6 * EIYL / EL;
SEP(11, 11) = 4 * EIYL;
SEP(12, 12) = 4 * EIZL;
for I = 1:12
   for J = I:12
      SEP(J, I) = SEP(I, J);
   end
end
%   CONVERT ELEMENT STIFFNESS MATRIX TO GLOBAL SYSTEM
DCOS(1, 1) = X21 / EL;
DCOS(1, 2) = Y21 / EL;
DCOS(1, 3) = Z21 / EL;
EIP1 = X(I3, 1) - X(I1, 1);
EIP2 = X(I3, 2) - X(I1, 2);
EIP3 = X(I3, 3) - X(I1, 3);
C1 = DCOS(1, 2) * EIP3 - DCOS(1, 3) * EIP2;
C2 = DCOS(1, 3) * EIP1 - DCOS(1, 1) * EIP3;
C3 = DCOS(1, 1) * EIP2 - DCOS(1, 2) * EIP1;
CC = sqrt(C1 * C1 + C2 * C2 + C3 * C3);
DCOS(3, 1) = C1 / CC;
DCOS(3, 2) = C2 / CC;
DCOS(3, 3) = C3 / CC;
DCOS(2, 1) = DCOS(3, 2) * DCOS(1, 3) - DCOS(1, 2) * DCOS(3, 3);
DCOS(2, 2) = DCOS(1, 1) * DCOS(3, 3) - DCOS(3, 1) * DCOS(1, 3);
DCOS(2, 3) = DCOS(3, 1) * DCOS(1, 2) - DCOS(1, 1) * DCOS(3, 2);
      
ALMBDA = zeros(12);
             
for K = 1:4
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
    SE = zeros(12);  % dummy
end