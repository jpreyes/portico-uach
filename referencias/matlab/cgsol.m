function []=cgsol();
disp('*********   PROGRAM CGSOL   **********')
disp('*     CONJUGATE GRADIENT METHOD      *')
disp('*   FOR SOLVING AX=B, A Symmetric    *')
disp('* T.R.Chandrupatla and A.D.Belegundu *')
disp('**************************************')

InputData;
CgSol;
Output;

%------------------------  function InputData  ---------------------------
function []=InputData()
global A B N LOUT FILE1

disp(blanks(1));
FILE1 = input('Input Data File Name ','s');
LINP  = fopen(FILE1,'r');
FILE2 = input('Output Data File Name ','s');
LOUT  = fopen(FILE2,'w');

DUMMY = fgets(LINP);
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[N]= deal(TMP(1));
DUMMY = fgets(LINP);
for I=1:N
   TMP = str2num(fgets(LINP));
   [A(I,:)] = deal(TMP(1:N));
end
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[B(1:N)] = deal(TMP(1:N));
fclose(LINP);      
    
%------------------------  function CgSol  ---------------------------
function [] = CgSol()
global A B N LOUT FILE1
global X ITER      

% function CG call
ITER = 0;
[X, ITER] = CG(A, B, N, ITER);

%------------------------  function Output  ---------------------------
function [] = Output()
global A B N LOUT FILE1
global X ITER      

fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);
disp(sprintf('Iterations = %i', ITER));
fprintf('Iterations = %i', ITER);
disp('Solution Vector');
disp(sprintf('%12.4E\n',X));
for I=1:N
  fprintf(LOUT,'%15.4E\n',[X(I)]);
end
fclose(LOUT);
      
%------------------------  function CG  ---------------------------
function [X,ITER] = CG (A, B, N, ITER)
% Conjugate Gradiant
for I = 1:N
   X(I) = 0.;
   G(I) = -B(I);
   D(I) = B(I);
end
GG1 = 0.;
for I = 1:N
   GG1 = GG1 + G(I) * G(I);
end

while GG1 >= 1e-6
	ITER = ITER + 1;
	DAD = 0.;
   for I = 1:N
	   C = 0.;
      for J = 1:N
          C = C + A(I, J) * D(J);
      end
      AD(I) = C;
      DAD = DAD + C * D(I);
   end
   AL = GG1 / DAD;
	GG2 = 0.;
   for I = 1:N
      X(I) = X(I) + AL * D(I);
      G(I) = G(I) + AL * AD(I);
      GG2 = GG2 + G(I) * G(I);
   end
   BT = GG2 / GG1;
   for I = 1:N
      D(I) = -G(I) + BT * D(I);
   end
   GG1 = GG2;
end
