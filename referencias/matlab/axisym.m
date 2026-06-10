function []=axisym2n()
clear all
close all
%------------------------  AXISYM2  ---------------------------
disp('=======================================');
disp('         PROGRAM AXISYM2               ');
disp('    AXISYMMETRIC STRESS ANALYSIS       ');
disp('    		WITH TEMPERATURE 			      ');
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
global X NOC F AREA MAT TH DT S
global PM NU U MPC BT STRESS REACT PLTSTRESS
global CNST
global TITLE FILE1 FILE2
global LINP LOUT LOUT2
global NQ IPL

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

%----- TOTAL DOF IS "NQ"
NQ = NDN * NN;

DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));


disp(blanks(1));
disp('PLOT CHOICE');
disp('  1) No Plot Data');
disp('  2) Create Data File for Radial Stress');
disp('  3) Create Data File for Z- Stress');
disp('  4) Create Data File for Hoop Stress');
disp('  5) Create Data File for von Mises Stress');
IPL = input('  Choose Number <1 to 5>: ');
%     --- default is no data
if isempty(IPL) | IPL<1 | IPL>5
   IPL = 1;
end
if IPL > 1 
    disp(blanks(1));
    FILE3 = input('Give Data File Name for Element Stresses ','s');
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
   [N, NOC(N,:), MAT(N), DT(N)] = ...
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
   [N, PM(N,:)] = deal(TMP(1), TMP(2:4));
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
global X NOC F AREA MAT TH DT S
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
global X NOC F AREA MAT TH DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- Global Stiffness Matrix
S = zeros(NQ,NBW);
for N = 1:NE
   disp(sprintf('Forming Stiffness Matrix of Element %d', N));
   [DJ, D, B, DB, RBAR] = dbmat(N, MAT, PM, NOC, X);
%  --- Element Stiffness
   for I = 1:6
      for J = 1:6
         C = 0;
         for K = 1:4
            C = C + abs(DJ) * B(K, I) * DB(K, J) * pi * RBAR;
         end
         SE(I, J) = C;
      end
   end
%  --- Temperature Load Vector
   M=MAT(N);
   PNU = PM(M, 2);
   AL = PM(M, 3);
   C = AL * DT(N) * pi * RBAR * abs(DJ);
   for I = 1:6
      TL(I) = C * (DB(1, I) + DB(2, I) + DB(4, I));
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
         F(NR) = F(NR) + TL(I);
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

%------------------------  function BandSolver  ---------------------------
function []=BandSolver();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT TH DT S
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
global X NOC F AREA MAT TH DT S
global PM NU U MPC BT STRESS PRINSTRESS REACT PLTSTRESS
global CNST
global TITLE FILE1 FILE2
global LINP LOUT IPL

for N = 1:NE
   [DJ, D, B, DB, RBAR] = dbmat(N, MAT, PM, NOC, X);
%   ----- Stress Evaluation
   M = MAT(N);
   PNU = PM(M, 2);
   AL = PM(M, 3);
   I1 = NOC(N, 1);
   I2 = NOC(N, 2);
   I3 = NOC(N, 3);
   Q(1) = F(2 * I1 - 1);
   Q(2) = F(2 * I1);
   Q(3) = F(2 * I2 - 1);
   Q(4) = F(2 * I2);
   Q(5) = F(2 * I3 - 1);
   Q(6) = F(2 * I3);
   C1 = AL * DT(N);
   for I = 1:4
      C = 0;
      for K = 1:6
         C = C + DB(I, K) * Q(K);
      end
      STR(I) = C - C1 * (D(I, 1) + D(I, 2) + D(I, 4));
   end
%  --- Principal Stress Calculations
   if STR(3) == 0
      S1 = STR(1);
      S2 = STR(2);
      ANG = 0;
      if S2 > S1
         S1 = STR(2);
         S2 = STR(1);
         ANG = 90;
      end
   else
      C = .5 * (STR(1) + STR(2));
      R = sqrt(.25 * (STR(1) - STR(2))^ 2 + (STR(3))^ 2);
      S1 = C + R;
      S2 = C - R;
      if C > STR(1)
         ANG = 57.2957795 * atan(STR(3) / (S1 - STR(1)));
         if STR(3) > 0; ANG = 90 - ANG; end
         if STR(3) > 0; ANG = -90 - ANG; end
      else
         ANG = 57.29577951 * atan(STR(3) / (STR(1) - S2));
      end
   end
   STRESS(N, 1) = STR(1);
   STRESS(N, 2) = STR(2);
   STRESS(N, 3) = STR(3);
   STRESS(N, 4) = STR(4);
   PRINSTRESS(N, 1) = S1;
   PRINSTRESS(N, 2) = S2;
   PRINSTRESS(N, 3) = ANG;
   if IPL == 2 
      PLTSTRESS(N) = STR(1);
   elseif IPL == 3
      PLTSTRESS(N) = STR(2);
   elseif IPL == 4
       PLTSTRESS(N) = STR(4);
   elseif IPL == 5    
      S3 = STR(4);
      C = (S1 - S2) ^ 2 + (S2 - S3)^ 2 + (S3 - S1) ^ 2;
      PLTSTRESS(N) = sqrt(.5 * C);
   end
end


%------------------------  function ReactionCalc  ---------------------------
function []=ReactionCalc();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT TH DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
for I = 1:ND
   N = NU(I);
   REACT(I) = CNST * (U(I) - F(N));
end

%------------------------  function Output  ---------------------------
function []=Output();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT TH DT S
global PM NU U MPC BT STRESS PRINSTRESS  REACT PLTSTRESS
global CNST
global TITLE FILE1 FILE2
global LINP LOUT LOUT2 IPL

disp(sprintf('Output for Input Data from file %s\n',FILE1));
fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);

disp(' Node#    R-Displ         Z-Displ');
fprintf(LOUT,' Node#    R-Displ         Z-Displ\n');
for I=1:NN;
	% print a matrix
	disp(sprintf(' %4d %15.4E %15.4E',[I,F(2*I-1),F(2*I)]'));
	fprintf(LOUT,' %4d %15.4E %15.4E\n',[I,F(2*I-1),F(2*I)]');
end
%----- Reaction Calculation -----
disp(sprintf('  DOF#     Reaction'));
fprintf(LOUT,'  DOF#     Reaction\n');
for I = 1:ND
   N = NU(I);
   disp(sprintf(' %4d %15.4E',N,REACT(I)));
   fprintf(LOUT,' %4d %15.4E\n',N,REACT(I));
end

%-----  Stress Calculations
if IPL ==2
   fprintf(LOUT2,'Radial Stress\n');
elseif IPL ==3
	fprintf(LOUT2,'Z- Stress\n');
elseif IPL ==4
	fprintf(LOUT2,'Hoop Stress\n');
elseif IPL ==5
	fprintf(LOUT2,'von Mises Stress\n');
end

disp(sprintf('%5s %14s %14s %14s %14s %14s %14s %14s',...
   			 'ELEM#','SR','SZ','TRZ','ST','S1','S2','ANGLE SR-->S1'));
fprintf(LOUT,'%5s %14s %14s %14s %14s %14s %14s %14s\n',...
   			 'ELEM#','SR','SZ','TRZ','ST','S1','S2','ANGLE SR-->S1');
for N = 1:NE
   disp(sprintf('%5d %14.3E %14.3E %14.3E %14.3E %14.3E %14.3E %14.3E', ...
      				N, STRESS(N,1),STRESS(N,2),STRESS(N,3),STRESS(N,4), ...
      				PRINSTRESS(N,1),PRINSTRESS(N,2),PRINSTRESS(N,3)));
   fprintf(LOUT,'%5d %14.5E %14.5E %14.5E %14.5E %14.3E %14.5E %14.5E\n', ...
      				N, STRESS(N,1),STRESS(N,2),STRESS(N,3),STRESS(N,4), ...
      				PRINSTRESS(N,1),PRINSTRESS(N,2),PRINSTRESS(N,3));

   if IPL > 1 
      fprintf(LOUT2,'%14.5E\n',PLTSTRESS(N));
   end
end
if IPL>1
   fclose(LOUT2);
end

fclose(LOUT);
disp(sprintf('Results are now available in the text file %s', FILE2));


%------------------------  dbmat  ---------------------------
function [DJ, D, B, DB, RBAR] = dbmat(N, MAT, PM, NOC, X);
%----- D(), B() and DB() matrices

%--- First the D-Matrix
      M = MAT(N);
      E = PM(M, 1);
      PNU = PM(M, 2);
      AL = PM(M, 3);
      C1 = E * (1 - PNU) / ((1 + PNU) * (1 - 2 * PNU));
      C2 = PNU / (1 - PNU);
      
      D = zeros(4);
      
      D(1, 1) = C1;
      D(1, 2) = C1 * C2;
      D(1, 4) = C1 * C2;
      D(2, 1) = D(1, 2);
      D(2, 2) = C1;
      D(2, 4) = C1 * C2;
      D(3, 3) = .5 * E / (1 + PNU);
      D(4, 1) = D(1, 4);
      D(4, 2) = D(2, 4);
      D(4, 4) = C1;
      
%--- Strain-Displacement Matrix B()
      I1 = NOC(N, 1);
      I2 = NOC(N, 2);
      I3 = NOC(N, 3);
      R1 = X(I1, 1);
      Z1 = X(I1, 2);
      R2 = X(I2, 1);
      Z2 = X(I2, 2);
      R3 = X(I3, 1);
      Z3 = X(I3, 2);
      R21 = R2 - R1;
      R32 = R3 - R2;
      R13 = R1 - R3;
      Z12 = Z1 - Z2;
      Z23 = Z2 - Z3;
      Z31 = Z3 - Z1;
      DJ = R13 * Z23 - R32 * Z31;
      RBAR = (R1 + R2 + R3) / 3;
%--- Definition of B() Matrix
      B(1, 1) = Z23 / DJ;
      B(2, 1) = 0;
      B(3, 1) = R32 / DJ;
      B(4, 1) = 1 / (3 * RBAR);
      B(1, 2) = 0;
      B(2, 2) = R32 / DJ;
      B(3, 2) = Z23 / DJ;
      B(4, 2) = 0;
      B(1, 3) = Z31 / DJ;
      B(2, 3) = 0;
      B(3, 3) = R13 / DJ;
      B(4, 3) = 1 / (3 * RBAR);
      B(1, 4) = 0;
      B(2, 4) = R13 / DJ;
      B(3, 4) = Z31 / DJ;
      B(4, 4) = 0;
      B(1, 5) = Z12 / DJ;
      B(2, 5) = 0;
      B(3, 5) = R21 / DJ;
      B(4, 5) = 1 / (3 * RBAR);
      B(1, 6) = 0;
      B(2, 6) = R21 / DJ;
      B(3, 6) = Z12 / DJ;
      B(4, 6) = 0;
      
%--- DB Matrix DB = D*B
for I = 1:4
   for J = 1:6
      C = 0;
      for K = 1:4
         C = C + D(I, K) * B(K, J);
      end
      DB(I, J) = C;
   end
end