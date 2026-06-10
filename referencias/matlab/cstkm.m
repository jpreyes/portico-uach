function []=cstkm()
clear all
close all
%------------------------  CSTKM2  ---------------------------
disp('=======================================');
disp('         PROGRAM CSTKM2                ');
disp('    2-D  CONSTANT STRAIN TRIANGLE      ');
disp('   T.R.Chandrupatla and A.D.Belegundu  ');
disp('=======================================');

InputData;
Bandwidth;
Stiffness;
ModifyForBC;
AddSprMass;
Output;

%------------------------  function InputData  ---------------------------
function [] = InputData();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT TH DT S GM
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT LOUT2
global NQ
global LC  

disp(' Problem Type');
disp('  1) Plane Stress Analysis');
disp('  2) Plane Strain Analysis');
LC = input('  Choose 1(default) or 2 :');
if isempty(LC) | LC<1 | LC>2
   LC = 1;
end

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

NPR = 4;  %  E,  NU,  ALPHA,  RHO

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
   [N,NOC(N,:), MAT(N,:), TH(N,:), DT(N,:)] = ...
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
global X NOC F AREA MAT TH DT S GM
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
global X NOC F AREA MAT TH DT S GM
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ
global LC 

%----- Global Stiffness Matrix
S = zeros(NQ,NBW);
GM = zeros(NQ,NBW);

for N = 1:NE
   disp(sprintf('Forming Stiffness Matrix of Element %d', N));
   %--------  Element Mass and Stiffness  -----
   [SE,EM] = elkm(LC,N,MAT,PM,X,NOC,TH);
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
                  GM(NR, NC) = GM(NR, NC) + EM(I, J);
               end
            end
         end
      end
   end
end

%------------------------  function ModifyForBC  ---------------------------
function []=ModifyForBC();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT TH DT S
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

%------------------------  function AddSprMass  ---------------------------
function []=AddSprMass();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT TH DT S GM
global PM NU U MPC BT  REACT
global CNST
global NQ

%-----  Additional Springs and Lumped Masses  -----
while 1    
	disp('SPRING SUPPORTS  < dof# = 0 Exits this mode >');
	N = input('dof# = ');
   if N == 0;break ;end;
	C = input('Spring Const');
   S(N, 1) = S(N, 1) + C;
end
while 1    
 	disp('LUMPED MASSES  < dof# = 0 Exits this mode >');
	N = input('dof# = ');
   if N == 0;break ;end
	C = input('Lumped Mass');
	GM(N, 1) = GM(N, 1) + C;
end

%------------------------  function Output  ---------------------------
function []=Output();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT TH DT S GM
global PM NU U MPC BT STRESS PRINSTRESS PLTSTRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT LOUT2
global LC IPL
global NQ

fprintf(LOUT,'Stiffness and Mass for Data in File %s\n',FILE1);
fprintf(LOUT,'Num. of DOF    Bandwidth\n');
fprintf(LOUT,'%4d %15d\n',NQ, NBW');
fprintf(LOUT,'BANDED STIFFNESS MATRIX\n');
for i=1:NQ
   fprintf(LOUT,'%15.7E  ',S(i,:));
   fprintf(LOUT,'\n');
end
fprintf(LOUT,'BANDED MASS MATRIX\n');
for i=1:NQ
   fprintf(LOUT,'%15.7E  ',GM(i,:));
   fprintf(LOUT,'\n');
end
fclose(LOUT);

disp(sprintf('Global Stiffness and Mass Matrices are in file %s', FILE2));
disp('Run INVITR  or  JACOBI program to get Eigenvalues and Eigenvectors');

%------------------------  function elkm  ---------------------------
function [SE,EM] = elkm(LC,N,MAT,PM,X,NOC,TH);
%----- D(), B() and DB() matrices
%--- Material Properties
M = MAT(N);
E = PM(M, 1);
PNU = PM(M, 2);
AL = PM(M, 3);
%--- D() Matrix
if (LC == 1)
   %--- Plane Stress
   C1 = E / (1 - PNU ^ 2);
   C2 = C1 * PNU;
else
   %--- Plane Strain
   C = E / ((1 + PNU) * (1 - 2 * PNU));
   C1 = C * (1 - PNU);
   C2 = C * PNU;
end
C3 = .5 * E / (1 + PNU);
      
      D(1, 1) = C1;
      D(1, 2) = C2;
      D(1, 3) = 0;
      D(2, 1) = C2;
      D(2, 2) = C1;
      D(2, 3) = 0;
      D(3, 1) = 0;
      D(3, 2) = 0;
      D(3, 3) = C3;
%--- Strain-Displacement Matrix B()
      I1 = NOC(N, 1);
      I2 = NOC(N, 2);
      I3 = NOC(N, 3);
      X1 = X(I1, 1);
      Y1 = X(I1, 2);
      X2 = X(I2, 1);
      Y2 = X(I2, 2);
      X3 = X(I3, 1);
      Y3 = X(I3, 2);
      X21 = X2 - X1;
      X32 = X3 - X2;
      X13 = X1 - X3;
      Y12 = Y1 - Y2;
      Y23 = Y2 - Y3;
      Y31 = Y3 - Y1;
      DJ = X13 * Y23 - X32 * Y31;
%--- Definition of B() Matrix
      B(1, 1) = Y23 / DJ;
      B(2, 1) = 0;
      B(3, 1) = X32 / DJ;
      B(1, 2) = 0;
      B(2, 2) = X32 / DJ;
      B(3, 2) = Y23 / DJ;
      B(1, 3) = Y31 / DJ;
      B(2, 3) = 0;
      B(3, 3) = X13 / DJ;
      B(1, 4) = 0;
      B(2, 4) = X13 / DJ;
      B(3, 4) = Y31 / DJ;
      B(1, 5) = Y12 / DJ;
      B(2, 5) = 0;
      B(3, 5) = X21 / DJ;
      B(1, 6) = 0;
      B(2, 6) = X21 / DJ;
      B(3, 6) = Y12 / DJ;
      
%--- DB Matrix DB = D*B
for I = 1:3
   for J = 1:6
      C = 0;
      for K = 1:3
         C = C + D(I, K) * B(K, J);
      end
      DB(I, J) = C;
   end
end
%--- Element Stiffness SE()
for I = 1:6
   for J = 1:6
      C = 0;
      for K = 1:3
              C = C + .5 * abs(DJ) * B(K, I) * DB(K, J) * TH(N);
 		end
      SE(I, J) = C;
   end
end
%--- Element Mass  EM()
RHO = PM(M, 4);
CM = RHO * TH(N) * .5 * abs(DJ) / 12;
EM = zeros(6);

%--- Non-zero elements of mass matrix are defined
      EM(1, 1) = 2 * CM;
      EM(1, 3) = CM;
      EM(1, 5) = CM;
      EM(2, 2) = 2 * CM;
      EM(2, 4) = CM;
      EM(2, 6) = CM;
      EM(3, 1) = CM;
      EM(3, 3) = 2 * CM;
      EM(3, 5) = CM;
      EM(4, 2) = CM;
      EM(4, 4) = 2 * CM;
      EM(4, 6) = CM;
      EM(5, 1) = CM;
      EM(5, 3) = CM;
      EM(5, 5) = 2 * CM;
      EM(6, 2) = CM;
      EM(6, 4) = CM;
      EM(6, 6) = 2 * CM;
