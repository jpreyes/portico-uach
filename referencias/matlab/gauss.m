function []=gauss();
close all;clear all;

%*****  Gauss Elimination Method (General Matrix) *****
%                   (no pivoting)
%*****           Row-wise Elimination             *****
disp('=======================================');
disp('          PROGRAM GAUSS                ');
disp('  T.R.Chandrupatla and A.D.Belegundu   ');
disp('=======================================');

InputData;
GaussRow;
Output;

%------------------------  function InputData  ---------------------------
function []=InputData();
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
%------------------------  function GaussRow  ---------------------------

function []=GaussRow();
global A B N LOUT FILE1
%----- Forward Elimination -----
for K = 1:N-1
   K1 = K + 1;
   for I = K1:N
      C = A(I, K) / A(K, K);
      for J = K1:N
         A(I, J) = A(I, J) - C * A(K, J);
      end
   B(I) = B(I) - C * B(K);
   end
end
%----- Back Substitution -----
B(N) = B(N) / A(N, N);
for II = 1:N-1
   I = N - II;
   I1 = I + 1;
   C = 1 / A(I, I);
   B(I) = C * B(I);
   for K = I1:N
      B(I) = B(I) - C * A(I, K) * B(K);
   end
end
  
%------------------------  function Output  ---------------------------
function []=Output();
global A B N LOUT FILE1

fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);
disp('Solution Vector');
disp(sprintf('%12.4E\n',B'));
for I=1:N
  fprintf(LOUT,'%15.4E\n',[B(I)]);
end
fclose(LOUT);
