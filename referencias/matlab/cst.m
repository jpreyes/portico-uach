function []=cst2n()
clear all
close all
%------------------------  CST2  ---------------------------
disp('=======================================');
disp('         PROGRAM CST2                  ');
disp('    2-D  CONSTANT STRAIN TRIANGLE      ');
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
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT LOUT2
global NQ
global LC IPL

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

NPR=3;  %E, NU,  ALPHA

disp(blanks(1));
disp('PLOT CHOICE');
disp('  1) No Plot Data');
disp('  2) Create Data File for in-plane Shear Stress');
disp('  3) Create Data File for Von Mises Stress');
IPL = input('  Choose 1(defalut), 2, or 3 :');
%     --- default is no data
if isempty(IPL) | IPL<1 | IPL>3
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
   [N,NOC(N,:), MAT(N,:), TH(N,:), DT(N,:)] = ...
      deal(TMP(1),TMP(2:1+NEN), TMP(2+NEN), TMP(3+NEN), TMP(4+NEN));
end

%----- Specified Displacements -----
DUMMY = fgets(LINP);
for I=1:ND
   TMP = str2num(fgets(LINP));
   [NU(I),U(I)] = deal(TMP(1), TMP(2));
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
global LC IPL

%----- Global Stiffness Matrix
S = zeros(NQ,NBW);
for N = 1:NE
   disp(sprintf('Forming Stiffness Matrix of Element %d', N));
   [DJ, D, B, DB] = dbmat(N, LC, MAT, PM, NOC, X);
%  --- Element Stiffness
   for I = 1:6
      for J = 1:6
         C = 0;
         for K = 1:3
            C = C + .5 * abs(DJ) * B(K, I) * DB(K, J) * TH(N);
         end
         SE(I, J) = C;
      end
   end
%  --- Temperature Load Vector
   M=MAT(N);
   PNU = PM(M, 2);
   AL = PM(M, 3);
   C = AL * DT(N);
   if (LC == 2); C = C * (1 + PNU); end
   for I = 1:6
      TL(I) = .5 * C * TH(N) * abs(DJ) * (DB(1, I) + DB(2, I));
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
global PM NU U MPC BT STRESS PRINSTRESS PLTSTRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global LC IPL

%-----  Stress Calculations
for N = 1:NE
   [DJ, D, B, DB] = dbmat(N, LC, MAT, PM, NOC, X);
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
   if LC == 2; C1 = C1 * (1 + PNU);end
   for I = 1:3
      C = 0;
      for K = 1:6
         C = C + DB(I, K) * Q(K);
      end
      STR(I) = C - C1 * (D(I, 1) + D(I, 2));
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
   PRINSTRESS(N, 1) = S1;
   PRINSTRESS(N, 2) = S2;
   PRINSTRESS(N, 3) = ANG;
   if IPL == 2 
      PLTSTRESS(N) = .5 * (S1 - S2);
   elseif IPL == 3
      S3 = 0;
      if LC == 2; S3 = PNU * (S1 + S2); end
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
global PM NU U MPC BT STRESS PRINSTRESS PLTSTRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT LOUT2
global LC IPL

disp(sprintf('Output for Input Data from file %s\n',FILE1));
fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);
if LC == 1; fprintf(LOUT,'Plane Stress Analysis\n'); end
if LC == 2; fprintf(LOUT,'Plane Strain Analysis\n'); end


%----- Displacements -----
disp(' Node#    X-Displ         Y-Displ');
fprintf(LOUT,' Node#    X-Displ         Y-Displ\n');
for I=1:NN;
	% print a matrix
	disp(sprintf(' %4d %15.4E %15.4E',[I,F(2*I-1),F(2*I)]'));
	fprintf(LOUT,' %4d %15.4E %15.4E\n',[I,F(2*I-1),F(2*I)]');
end
%----- Reaction  -----
disp(sprintf('  DOF#     Reaction'));
fprintf(LOUT,'  DOF#     Reaction\n');
for I = 1:ND
   N = NU(I);
   disp(sprintf(' %4d %15.4E',N,REACT(I)));
   fprintf(LOUT,' %4d %15.4E\n',N,REACT(I));
end

if IPL ==2
   fprintf(LOUT2,'Max. in-plane Shear Stress\n');
elseif IPL ==3
	fprintf(LOUT2,'Von Mises Stress\n');
end
%-----  Stress Calculations
disp(sprintf('%5s %14s %14s %14s %14s %14s %14s',...
   			 'ELEM#','SX','SY','TXY','S1','S2','ANGLE SX-->S1'));
fprintf(LOUT,'%5s %14s %14s %14s %14s %14s %14s\n',...
   			 'ELEM#','SX','SY','TXY','S1','S2','ANGLE SX-->S1');
for N = 1:NE
   disp(sprintf('%5d %14.3E %14.3E %14.3E %14.3E %14.3E %14.3E', ...
      				N, STRESS(N,1),STRESS(N,2),STRESS(N,3), ...
      				PRINSTRESS(N,1),PRINSTRESS(N,2),PRINSTRESS(N,3)));
   fprintf(LOUT,'%5d %14.5E %14.5E %14.5E %14.5E %14.5E %14.5E\n', ...
      				N, STRESS(N,1),STRESS(N,2),STRESS(N,3), ...
      				PRINSTRESS(N,1),PRINSTRESS(N,2),PRINSTRESS(N,3));
   if IPL == 2 
      fprintf(LOUT2,'%14.5E\n',PLTSTRESS(N));
   elseif IPL == 3
      fprintf(LOUT2,'%14.5E\n',PLTSTRESS(N));
   end
end
if IPL>1
   fclose(LOUT2);
end
fclose(LOUT);
disp(sprintf('Results are now available in the text file %s', FILE2));


%------------------------  function dbmat  ---------------------------
function [DJ, D, B, DB] = dbmat(N, LC, MAT, PM, NOC, X);
%----- D(), B() and DB() matrices

%  --- Material Properties
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