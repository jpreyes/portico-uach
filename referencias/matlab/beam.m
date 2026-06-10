function []=beam2n()
clear all
close all
clear
close
%------------------------ BEAM2  ---------------------------
disp('========================================');
disp('            PROGRAM BEAM2               ');
disp('        Beam Bending Analysis           ');
disp('   T.R.Chandrupatla and A.D.Belegundu   ');
disp('========================================');

InputData;
Bandwidth;
Stiffness;
ModifyForBC;
BandSolver;
ReactionCalc;
Output;

%------------------------  function InputData  ---------------------------
function [] = InputData();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT SMI S
global PM NU U MPC BT REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
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

NPR=1;   % E

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
   [N,NOC(N,:), MAT(N), SMI(N)] = ...
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
global X NOC F AREA MAT SMI S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
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
global X NOC F AREA MAT SMI S
global PM NU U MPC BT  REACT
global CNST
global NQ
%----- Global Stiffness Matrix
S = zeros(NQ,NBW);

for N = 1:NE
   disp(sprintf('Forming Stiffness Matrix of Element %d', N));
   
%--------  Element Stiffness and Temperature Load  -----
   N1 = NOC(N, 1);
   N2 = NOC(N, 2);
	M = MAT(N);
   EL = abs(X(N1) - X(N2));
   EIL = PM(M, 1) * SMI(N) / EL^3;
   SE(1, 1) = 12 * EIL;
   SE(1, 2) = EIL * 6 * EL;
   SE(1, 3) = -12 * EIL;
   SE(1, 4) = EIL * 6 * EL;
   SE(2, 1) = SE(1, 2);
   SE(2, 2) = EIL * 4 * EL * EL;
   SE(2, 3) = -EIL * 6 * EL;
   SE(2, 4) = EIL * 2 * EL * EL;
   SE(3, 1) = SE(1, 3);
   SE(3, 2) = SE(2, 3);
   SE(3, 3) = EIL * 12;
   SE(3, 4) = -EIL * 6 * EL;
   SE(4, 1) = SE(1, 4);
   SE(4, 2) = SE(2, 4);
   SE(4, 3) = SE(3, 4);
   SE(4, 4) = EIL * 4 * EL * EL;
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

%------------------------  function ModifyForBC  ---------------------------
function []=ModifyForBC();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT SMI S
global PM NU U MPC BT  REACT
global CNST
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
global X NOC F AREA MAT SMI S
global PM NU U MPC BT  REACT
global CNST
global NQ
%----- Equation Solving using Band Solver -----
disp('Solving using Band Solver(bansol.m)');
[F] = bansol(NQ,NBW,S,F);


%------------------------  function ReactionCalc  ---------------------------
function []=ReactionCalc();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT SMI S
global PM NU U MPC BT  REACT
global CNST
for I = 1:ND
   N = NU(I);
   REACT(I) = CNST * (U(I) - F(N));
end

%------------------------  function Output  ---------------------------
function []=Output();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT SMI S
global PM NU U MPC BT REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT 

disp(sprintf('Output for Input Data from file %s\n',FILE1));
fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);

disp(' Node#    X-Displ         Rotation');
fprintf(LOUT,' Node#    X-Displ         Rotation\n');
I=[1:NN]';
% print a matrix
disp(sprintf(' %4d %15.4E %15.4E\n',[I,F(2*I-1),F(2*I)]'));
fprintf(LOUT,' %4d %15.4E %15.4E\n',[I,F(2*I-1),F(2*I)]');

%----- Reaction Calculation -----
disp(sprintf('  DOF#     Reaction'));
fprintf(LOUT,'  DOF#     Reaction\n');
for I = 1:ND
   N = NU(I);
   R = CNST * (U(I) - F(N));
   disp(sprintf(' %4d %15.4E',N,REACT(I)));
   fprintf(LOUT,' %4d %15.4E\n',N,REACT(I));
end

fclose(LOUT);
disp(sprintf('The Results are available in the text file %s', FILE2));

