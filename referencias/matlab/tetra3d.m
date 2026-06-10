function [] = tetra3d()
clear all
close all
%------------------------ TETRA2  ---------------------------
disp('==========================================');
disp('           PROGRAM TETRA2                 ');
disp('    Three Dimensional Stress Analysis     ');
disp('         Tetrahedral Elements             ');
disp('   T.R.Chandrupatla and A.D.Belegundu     ');
disp('==========================================');

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
global X NOC F AREA MAT TH DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2 FILE3
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

%----- Total dof is  NQ
NQ = NDN * NN;

DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));

NPR=3;   %  E, NU, ALPHA

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
   [N,NOC(N,:), MAT(N,:), DT(N,:)] = ...
      deal(TMP(1),TMP(2:1+NEN), TMP(2+NEN), TMP(3+NEN));
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
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- Global Stiffness Matrix
S = zeros(NQ,NBW);

for N = 1:NE
   disp(sprintf('Forming Stiffness Matrix of Element %d', N));
   
   [DJ, D, B, DB] = dbmat(N, MAT, PM, NOC, X);
   
%--------  Element Stiffness -----
   SE = B' * DB * abs(DJ) / 6;
%-------- Temperature Load Vector QT() -----
   AL = PM(MAT(N), 3);
   C = AL * DT(N);
   for I = 1:12
      DSUM = DB(1, I) + DB(2, I) + DB(3, I);
      QT(I) = C * abs(DJ) * DSUM / 6;
   end 	  
   
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
         F(NR) = F(NR) + QT(I);
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
global PM NU U MPC BT STRESS PSTRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
%-----  Stress Calculations -----
%--- Stresses at Integration Points
for N = 1:NE
   [DJ, D, B, DB] = dbmat(N, MAT, PM, NOC, X);
   AL = PM(MAT(N), 3);

%  --- Stress Evaluation (Element Nodal Displacements stored in QT() )
   for I = 1:4
      IN = 3 * (NOC(N, I) - 1);
      II = 3 * (I - 1);
      for J = 1:3
         QT(II + J) = F(IN + J);
      end
   end
   C1 = AL * DT(N);
   for I = 1:6
      STR(I) = 0;
      for K = 1:12
         STR(I) = STR(I) + DB(I, K) * QT(K);
      end
      STR(I) = STR(I) - C1 * (D(I, 1) + D(I, 2) + D(I, 3));
   end

%  --- Principal Stress Calculations
   AI1 = STR(1) + STR(2) + STR(3);
   AI21 = STR(1) * STR(2) + STR(2) * STR(3) + STR(3) * STR(1);
   AI22 = STR(4) * STR(4) + STR(5) * STR(5) + STR(6) * STR(6);
   AI2 = AI21 - AI22;
   AI31 = STR(1) * STR(2) * STR(3) + 2 * STR(4) * STR(5) * STR(6);
   AI32 = STR(1)*STR(4)^2 + STR(2)*STR(5)^2 + STR(3)*STR(6)^2;
   AI3 = AI31 - AI32;
   C1 = AI2 - AI1^2 / 3;
   C2 = -2 * (AI1 / 3)^3 + AI1 * AI2 / 3 - AI3;
   C3 = 2 * sqrt(-C1 / 3);
   TH = -3 * C2 / (C1 * C3);
   TH2 = abs(1 - TH * TH);
   if (TH == 0); TH = pi / 2; end
   if (TH > 0) ; TH = atan(sqrt(TH2) / TH); end
   if (TH < 0) ; TH = pi - atan(sqrt(TH2) / TH); end
   TH = TH / 3;
%  --- Principal Stresses
   P1 = AI1 / 3 + C3 * cos(TH);
   P2 = AI1 / 3 + C3 * cos(TH + 2 * pi / 3);
   P3 = AI1 / 3 + C3 * cos(TH + 4 * pi / 3);
   
   for I = 1:6
      STRESS(N, I) = STR(I);
   end
   PSTRESS(N, 1) = P1;
   PSTRESS(N, 2) = P2;
   PSTRESS(N, 3) = P3;
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
global PM NU U MPC BT STRESS PSTRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT

%----------------------------------------------
disp(sprintf('Output for Input Data from file %s\n',FILE1));
fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);

disp(' Node#    X-Displ         Y-Displ         Z-Displ');
fprintf(LOUT,' Node#    X-Displ         Y-Displ         Z-Displ\n');
I=[1:NN]';
% print a matrix
disp(sprintf(' %4d %15.4E %15.4E %15.4E\n',[I,F(3*I-2),F(3*I-1),F(3*I)]'));
fprintf(LOUT,' %4d %15.4E %15.4E %15.4E\n',[I,F(3*I-2),F(3*I-1),F(3*I)]');


%-----  Stress Calculations -----
for N = 1:NE   
   fprintf(LOUT,' Stresses in element No. %d\n',N);
   fprintf(LOUT,'     Normal Stresses    SX,  SY,  SZ  :');
   fprintf(LOUT,' %15.4E %15.4E %15.4E \n',STRESS(1), STRESS(2), STRESS(3));
   fprintf(LOUT,'     Shear Stresses    TYZ, TXZ, TXY  :');
   fprintf(LOUT,' %15.4E %15.4E %15.4E \n',STRESS(4), STRESS(5), STRESS(6));
   fprintf(LOUT,'     Principal Stresses P1,  P2,  P3  :');
   fprintf(LOUT,' %15.4E %15.4E %15.4E \n',PSTRESS(N,1),PSTRESS(N,2),PSTRESS(N,3));
end

%----- Reaction Calculation -----
disp(sprintf('  DOF#     Reaction'));
fprintf(LOUT,'  DOF#     Reaction\n');
for I = 1:ND
   N = NU(I);
   disp(sprintf(' %4d %15.4E',N,REACT(I)));
   fprintf(LOUT,' %4d %15.4E\n',N,REACT(I));
end
fclose(LOUT);
disp(sprintf('The Results are available in the text file %s', FILE2));


%------------------------  dbmat  ---------------------------
function [DJ, D, B, DB] = dbmat(N, MAT, PM, NOC, X);

%  --- Material Properties
M = MAT(N);
E = PM(M, 1);
PNU = PM(M, 2);
AL = PM(M, 3);
C4 = E / ((1 + PNU) * (1 - 2 * PNU));
C1 = C4 * (1 - PNU);
C2 = C4 * PNU;
C3 = .5 * E / (1 + PNU);
      
D = zeros(6);
D(1, 1) = C1;
D(1, 2) = C2;
D(1, 3) = C2;
D(2, 1) = C2;
D(2, 2) = C1;
D(2, 3) = C2;
D(3, 1) = C2;
D(3, 2) = C2;
D(3, 3) = C1;
D(4, 4) = C3;
D(5, 5) = C3;
D(6, 6) = C3;

%   --- Strain-Displacement Matrix B()
I1 = NOC(N, 1);
I2 = NOC(N, 2);
I3 = NOC(N, 3);
I4 = NOC(N, 4);
X14 = X(I1, 1) - X(I4, 1);
X24 = X(I2, 1) - X(I4, 1);
X34 = X(I3, 1) - X(I4, 1);
Y14 = X(I1, 2) - X(I4, 2);
Y24 = X(I2, 2) - X(I4, 2);
Y34 = X(I3, 2) - X(I4, 2);
Z14 = X(I1, 3) - X(I4, 3);
Z24 = X(I2, 3) - X(I4, 3);
Z34 = X(I3, 3) - X(I4, 3);
DJ1 = X14 * (Y24 * Z34 - Z24 * Y34);
DJ2 = Y14 * (Z24 * X34 - X24 * Z34);
DJ3 = Z14 * (X24 * Y34 - Y24 * X34);
DJ = DJ1 + DJ2 + DJ3;
A11 = (Y24 * Z34 - Z24 * Y34) / DJ;
A21 = (Z24 * X34 - X24 * Z34) / DJ;
A31 = (X24 * Y34 - Y24 * X34) / DJ;
A12 = (Y34 * Z14 - Z34 * Y14) / DJ;
A22 = (Z34 * X14 - X34 * Z14) / DJ;
A32 = (X34 * Y14 - Y34 * X14) / DJ;
A13 = (Y14 * Z24 - Z14 * Y24) / DJ;
A23 = (Z14 * X24 - X14 * Z24) / DJ;
A33 = (X14 * Y24 - Y14 * X24) / DJ;

%   ---  B Matrix
B = zeros(6, 12);      
B(1, 1) = A11;
B(1, 4) = A12;
B(1, 7) = A13;
B(1, 10) = -A11 - A12 - A13;
B(2, 2) = A21;
B(2, 5) = A22;
B(2, 8) = A23;
B(2, 11) = -A21 - A22 - A23;
B(3, 3) = A31;
B(3, 6) = A32;
B(3, 9) = A33;
B(3, 12) = -A31 - A32 - A33;
B(4, 2) = A31;
B(4, 3) = A21;
B(4, 5) = A32;
B(4, 6) = A22;
B(4, 8) = A33;
B(4, 9) = A23;
B(4, 11) = B(3, 12);
B(4, 12) = B(2, 11);
B(5, 1) = A31;
B(5, 3) = A11;
B(5, 4) = A32;
B(5, 6) = A12;
B(5, 7) = A33;
B(5, 9) = A13;
B(5, 10) = B(3, 12);
B(5, 12) = B(1, 10);
B(6, 1) = A21;
B(6, 2) = A11;
B(6, 4) = A22;
B(6, 5) = A12;
B(6, 7) = A23;
B(6, 8) = A13;
B(6, 10) = B(2, 11);
B(6, 11) = B(1, 10);

%   --- DB Matrix DB = D*B
DB = D*B;
 
