function []=heat1d2n()
clear all
close all
%-----------------------  HEAT1D2  ------------------------
disp('***************************************');
disp('*          PROGRAM HEAT1D2            *');
disp('* T.R.Chandrupatla and A.D.Belegundu  *');
disp('***************************************');      

InputData;
Stiffness;
ModifyForBC;
BandSolver;
Output;

%------------------------  function InputData  ---------------------------
function [] = InputData();
global NN NE NBC NQ NBW TC
global X NB BC V H F S
global TITLE FILE1 FILE2
global LINP LOUT

disp(blanks(1));
FILE1 = input('Input Data File Name ','s');
LINP  = fopen(FILE1,'r');
FILE2 = input('Output Data File Name ','s');
LOUT  = fopen(FILE2,'w');

DUMMY = fgets(LINP);
TITLE = fgets(LINP);
DUMMY = fgets(LINP);

TMP = str2num(fgets(LINP));
[NE, NBC, NQ] = deal(TMP(1),TMP(2),TMP(3));
NN = NE + 1;
NBW = 2;
% NBW IS THE HALF-BAND-WIDTH

%----- Thermal Conductivity -----
DUMMY = fgets(LINP);
for I=1:NE
   TMP = str2num(fgets(LINP));
   [N, TC(N,1)]=deal(TMP(1),TMP(2));
end
%----- Coordinates -----
DUMMY = fgets(LINP);
for I=1:NN
%   TMP = str2num(fgets(LINP));
%   [N, X(N,1)]=deal(TMP(1),TMP(2));
	N = fscanf(LINP,'%d',[1 1]);
   X(N,1) = fscanf(LINP,'%g\n',[1 1]);
end

DUMMY = fgets(LINP);

%----- Boundary Conditions -----
for I = 1:NBC
   NB(I) = fscanf(LINP,'%d',[1 1]);
   BC(I,:) = upper(fscanf(LINP,'%s\n',[1 1]));
   if ~isempty(findstr(BC(I,:),'TEMP'))
      V(I) = fscanf(LINP,'%g\n',[1 1]);
   elseif ~isempty(findstr(BC(I,:),'FLUX'))
      V(I) = fscanf(LINP,'%g\n',[1 1]);
   elseif ~isempty(findstr(BC(I,:),'CONV'))
      H(I)= fscanf(LINP,'%g',[1 1]);
      V(I)= fscanf(LINP,'%g\n',[1 1]);
   end
end
%--- CALCULATE AND INPUT NODAL HEAT SOURCE VECTOR ---
F = zeros(NN,1);
DUMMY = fgets(LINP);
if NQ > 0
   for I = 1:NQ
      N = fscanf(LINP,'%d',[1 1]);
      F(N) = fscanf(LINP,'%g\n',[1 1]);
   end
end
fclose(LINP);

%------------------------  function Stiffness  ---------------------------
function []=Stiffness();
global NN NE NBC NQ NBW TC
global X NB BC V H F S
global TITLE FILE1 FILE2
global LINP LOUT

%--- STIFFNESS MATRIX ---
S = zeros(NN,NBW);
for I = 1:NE
   I1 = I;
   I2 = I + 1;
   ELL = abs(X(I2) - X(I1));
   EKL = TC(I) / ELL;
   S(I1, 1) = S(I1, 1) + EKL;
   S(I2, 1) = S(I2, 1) + EKL;
   S(I1, 2) = S(I1, 2) - EKL;
end

%------------------------  function ModifyForBC  ---------------------------
function []=ModifyForBC();
global NN NE NBC NQ NBW TC
global X NB BC V H F S
global TITLE FILE1 FILE2
global LINP LOUT

%--- ACCOUNT FOR B.C.'S ---
AMAX = 0;
for I = 1:NN
   if S(I, 1) > AMAX; AMAX = S(I, 1); end
end
CNST = AMAX * 10000.;
for I = 1:NBC
   N = NB(I);
   if ~isempty(findstr(BC(I,:),'CONV'))
      S(N, 1) = S(N, 1) + H(I);
      F(N) = F(N) + H(I) * V(I);
   elseif ~isempty(findstr(BC(I,:),'FLUX'))
      F(N) = F(N) - V(I);
   elseif ~isempty(findstr(BC(I,:),'TEMP'))
      S(N, 1) = S(N, 1) + CNST;
      F(N) = F(N) + CNST * V(I);
   end
end

%------------------------  function BandSolver  ---------------------------
function []=BandSolver();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global NN NE NBC NQ NBW TC
global X NB BC V H F S
global TITLE FILE1 FILE2
global LINP LOUT

[F] = bansol(NN,NBW,S,F);

%------------------------  function Output  ---------------------------
function []=Output();
global NN NE NBC NQ NBW TC
global X NB BC V H F S
global TITLE FILE1 FILE2
global LINP LOUT

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);
disp('NODE#   TEMPERATURE');
fprintf(LOUT,'NODE#   TEMPERATURE\n');
for I=1:NN
	disp(sprintf('%d %15.5G', I, F(I)));
	fprintf(LOUT,'%d %15.5G\n', I, F(I));
end
fclose(LOUT);
disp(sprintf('RESULTS ARE IN FILE : %s', FILE2));
