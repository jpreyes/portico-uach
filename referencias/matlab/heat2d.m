function []=heat2d()
clear all
close all
%-----------------------  HEAT1D2  ------------------------
disp('***************************************');
disp('*          PROGRAM HEAT2D2            *');
disp('*   HEAT 2-D  WITH 3-NODED TRIANGLES  *');
disp('* T.R.Chandrupatla and A.D.Belegundu  *');
disp('***************************************');      

InputData;
Bandwidth;
Stiffness;
ModifyForBC;
BandSolver;
HeatFlowCalc;
Output;

%------------------------  function InputData  ---------------------------
function [] = InputData();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT EHS F S
global PM NU U MPC BT Q
global CNST
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2
global NQ
global LC IPL

disp(blanks(1));
FILE1 = input('Input Data File Name ','s');
LINP  = fopen(FILE1,'r');
FILE2 = input('Output Data File Name ','s');
LOUT  = fopen(FILE2,'w');

DUMMY = fgets(LINP);
TITLE = fgets(LINP); % Read line without carrige return
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[NN, NE, NM, NDIM, NEN, NDN] = deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5),TMP(6));
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));

NPR = 1;
NMPC = 0;
NDN = 1;
NDIM=2;
NEN=3;
%--- NCH = 1 = ELEMENT HEAT SOURCE, EHS(I),I=1,...,NE
%--- ND = NO. OF SPECIFIED TEMPERATURES
%--- NL = NO. OF NODAL HEAT SOURCES
%---  NPR =1 (THERMAL CONDUCTIVITY) AND NMPC = 0

disp(blanks(1));
disp('PLOT CHOICE');
disp('  1) No Plot Data');
disp('  2) Create Data File for Nodal Temperatures');
IPL = input('  Choose 1(defalut) or 2 :');
%     --- default is no data
if isempty(IPL) | IPL<1 | IPL>2
   IPL = 1;
end
if IPL > 1 
    disp(blanks(1));
    FILE3 = input('Give Data File Name for Nodal Temperatures ','s');
    LOUT2  = fopen(FILE3,'w');
end


%----- Coordinates -----
DUMMY = fgets(LINP);
for I=1:NN
   TMP = str2num(fgets(LINP));
   [N, X(N,:)]=deal(TMP(1),TMP(2:1+NDIM));
%	N = fscanf(LINP,'%d',[1 1]);
%   X(N,1:NDIM) = fscanf(LINP,'%g\n',[NDIM 1])';
end
%----- Connectivity -----
DUMMY = fgets(LINP);
for I=1:NE
   TMP = str2num(fgets(LINP));
   [N,NOC(N,:), MAT(N,:), EHS(N,:)] = ...
      deal(TMP(1),TMP(2:1+NEN), TMP(2+NEN), TMP(3+NEN));
end

%----- Temperature BC -----
DUMMY = fgets(LINP);
for I=1:ND
   TMP = str2num(fgets(LINP));
   [NU(I,:),U(I,:)] = deal(TMP(1), TMP(2));
end

%----- Nodal Heat Sources -----
DUMMY = fgets(LINP);
F = zeros(NN,1);
for I=1:NL
   TMP = str2num(fgets(LINP));
   [N,F(N)]=deal(TMP(1),TMP(2));
end

%----- Thermal Conductivity -----
DUMMY = fgets(LINP);
for I=1:NM
   TMP = str2num(fgets(LINP));
   [N, PM(N,:)] = deal(TMP(1), TMP(2:NPR+1));
end

%------------------------  function Bandwidth  ---------------------------
function []=Bandwidth();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT EHS F S
global PM NU U MPC BT Q
global CNST
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2
global NQ
global LC IPL

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
disp(sprintf('Bandwidth = %d', NBW));

%------------------------  function Stiffness  ---------------------------
function [] = Stiffness();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT EHS F S
global PM NU U MPC BT Q
global CNST
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2
global NQ
global LC IPL
%----- INITIALIZATION OF CONDUCTIVITY MATRIX
S = zeros(NN,NBW);
DUMMY = fgets(LINP);
NHF = fscanf(LINP,'%d\n',[1 1]);
if (NHF > 0)
   for I=1:NHF
      TMP = str2num(fgets(LINP));
      [N1, N2, V] = deal(TMP(1), TMP(2), TMP(3));
      ELEN = sqrt((X(N1,1)-X(N2,1))^2 +(X(N1,2)-X(N2,2))^2);
      F(N1) = F(N1) - ELEN * V / 2;
      F(N2) = F(N2) - ELEN * V / 2;
   end
end
DUMMY = fgets(LINP);
NCONV = fscanf(LINP,'%d\n',[1 1]);
if (NCONV > 0)
   for I = 1:NCONV
      TMP = str2num(fgets(LINP));
      [N1, N2, H, TINF] = deal(TMP(1), TMP(2), TMP(3), TMP(4));
      ELEN = sqrt((X(N1,1)-X(N2,1))^2+(X(N1,2)-X(N2,2))^2);
      F(N1) = F(N1) + ELEN * H * TINF / 2;
      F(N2) = F(N2) + ELEN * H * TINF / 2;
      S(N1, 1) = S(N1, 1) + H * ELEN / 3;
      S(N2, 1) = S(N2, 1) + H * ELEN / 3;
      if (N1 >= N2)
         N3 = N1;
         N1 = N2;
         N2 = N3;
      end
      S(N1,N2-N1+1) = S(N1,N2-N1+1) + H * ELEN / 6;
   end
end
fclose(LINP);

%--- CONDUCTIVITY MATRIX
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
%  --- ELEMENT HEAT SOURCES
   if EHS(I) ~= 0
      C = EHS(I) * AREA / 3;
      F(I1) = F(I1) + C;
      F(I2) = F(I2) + C;
      F(I3) = F(I3) + C;
   end
   BT(1, 1) = Y23;
   BT(1, 2) = Y31;
   BT(1, 3) = Y12;
   BT(2, 1) = X32;
   BT(2, 2) = X13;
   BT(2, 3) = X21;
        
   BT = BT ./ DETJ;

   for II = 1:3
      for JJ = 1:3
         II1 = NOC(I, II);
         II2 = NOC(I, JJ);
         if II1 <= II2
            SUM = 0.;
            for J = 1:2
               SUM = SUM + BT(J, II) * BT(J, JJ);
            end        
            IC = II2 - II1 + 1;
            S(II1, IC) = S(II1, IC) + SUM * AREA * PM(MAT(I), 1);
         end
		end
	end
end

%------------------------  function ModifyForBC  ---------------------------
function []=ModifyForBC();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT EHS F S
global PM NU U MPC BT Q
global CNST
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2
global NQ
global LC IPL

if ND >0
%----- Decide Penalty Parameter CNST -----
   CNST = 0;
   for I = 1:NN
      if CNST < S(I, 1); CNST = S(I, 1); end
   end
   CNST = CNST * 10000;

%----- Modify for Boundary Conditions -----
%    --- Temperature BC ---
   for I = 1:ND
      N = NU(I);
      S(N, 1) = S(N, 1) + CNST;
      F(N) = F(N) + CNST * U(I);
   end
end

%------------------------  function BandSolver  ---------------------------
function []=BandSolver();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT EHS F S
global PM NU U MPC BT Q
global CNST
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2
global NQ
global LC IPL
%----- Equation Solving using Band Solver -----
disp('Solving using Band Solver(bansol.m)');
[F] = bansol(NN,NBW,S,F);

%------------------------  function HeatFlowCalc  ---------------------------
function []=HeatFlowCalc();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT EHS F S
global PM NU U MPC BT Q
global CNST
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2
global NQ
global LC IPL

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
   BT(1, 1) = Y23;
   BT(1, 2) = Y31;
   BT(1, 3) = Y12;
   BT(2, 1) = X32;
   BT(2, 2) = X13;
   BT(2, 3) = X21;
        
   BT = BT ./ DETJ;

   QX = BT(1, 1) * F(I1) + BT(1, 2) * F(I2) + BT(1, 3) * F(I3);
   QX = -QX * PM(MAT(I), 1);
   QY = BT(2, 1) * F(I1) + BT(2, 2) * F(I2) + BT(2, 3) * F(I3);
   QY = -QY * PM(MAT(I), 1);
   Q(I,1) = QX;
   Q(I,2) = QY;
end
%------------------------  function Output  ---------------------------
function []=Output();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC MAT EHS F S
global PM NU U MPC BT Q
global CNST
global TITLE FILE1 FILE2 FILE3
global LINP LOUT LOUT2
global NQ
global LC IPL

disp(sprintf('Output for Input Data from file %s\n',FILE1));
fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);

disp(' Node#    Temperature');
fprintf(LOUT,' Node#    Temperature\n');

I=[1:NN]';
disp(sprintf(' %4d %15.4f\n',[I,F(I)]'));
fprintf(LOUT,' %4d %15.4f\n',[I,F(I)]');
if IPL > 1
   fprintf(LOUT2,'(Node Temp.) for Data in File %s\n', FILE1);
   fprintf(LOUT2,'%15.4f\n',[F(I)]');
   fclose(LOUT2);
end

disp(' -- CONDUCTION HEAT FLOW PER UNIT AREA IN EACH ELEMENT -- ');
disp(' ELEMENT#    QX= -K*DT/DX    QY= -K*DT/DY ');
fprintf(LOUT,' -- CONDUCTION HEAT FLOW PER UNIT AREA IN EACH ELEMENT -- \n');
fprintf(LOUT,'ELEMENT#   QX= -K*DT/DX   QY= -K*DT/DY \n');
for I = 1:NE  
	disp(sprintf('	%d %15.5G %15.5G',I, Q(I,1), Q(I,2)));
	fprintf(LOUT, '	%d %15.5G %15.5G\n',I, Q(I,1), Q(I,2));
end
fclose(LOUT);
disp(sprintf('Complete results are in file : %s',FILE2));
if IPL > 1 
   disp(sprintf('Element Stress Data in file : %s',FILE3));
   disp('Run CONTOUR1/2 to plot ISOTHERMS');
end
      
