%------------------------  function InputData  ---------------------------
function []=skyline2n();
%   ****  Sky Line Method  ****
%   **** Gauss Elimination ****
disp(' =======================================');
disp('            PROGRAM SKYLINE             ');
disp('   T.R.Chandrupatla and A.D.Belegundu   ');
disp(' =======================================');

InputData;
SkyLine;
Output;

%------------------------  function InputData  ---------------------------
function []=InputData();
global A B ID N
global LOUT  FILE1
 
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
TMP = str2num(fgets(LINP));
[ID(1:N)] = deal(TMP(1:N));
NSUM = 0;
for I = 2:N
   ID(I) = ID(I) + ID(I - 1);
end
NSUM = ID(N); 
DUMMY = fgets(LINP);
NT = 0; NI = 1;
for I = 1:N
   if I ~= 1
      NT = ID(I - 1);
      NI = ID(I) - ID(I - 1);
   end
   TMP = str2num(fgets(LINP));
   [A(NT+1:NT+NI)] = deal(TMP(1:NI));
end
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[B(1:N)] = deal(TMP(1:N));
fclose(LINP);      

%------------------------  function SkyLine  ---------------------------
function []=SkyLine();
global A B ID N

B = skylin(A, B, ID, N);

%------------------------  function Output  ---------------------------
function []=Output();
global A B ID N 
global LOUT  FILE1
fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);
disp('Solution Vector');
disp(sprintf('%12.4E\n',B'));
for I=1:N
  fprintf(LOUT,'%15.4E\n',[B(I)]);
end
fclose(LOUT);

%------------------------  function skylin  ---------------------------
function F = skylin(A, F, ID, NQ);
% Skyline Solver Function

%--- Forward Elimination ---
for J = 2:NQ
   NJ = ID(J) - ID(J - 1);
   if NJ ~= 1
   	K1 = 0;
  	 	NJ = J - NJ + 1;
   	for K = NJ:J - 1
    		K1 = K1 + 1;
     		KJ = ID(J - 1) + K1;
     		KK = ID(K);
   	  	C = A(KJ) / A(KK);
   	  	for I = K + 1:J
     			NI = ID(I) - ID(I - 1);
     			if (I - K + 1) <= NI
     				IJ = ID(J) - J + I;
     				KI = ID(I) - I + K;
               A(IJ) = A(IJ) - C * A(KI);
            end
  			end
      end
   end
end

for K = 1:NQ-1
   KK = ID(K);
   C = F(K) / A(KK);
   for I = K+1:NQ
      NI = ID(I) - ID(I - 1);
      if (I - K + 1) <= NI
      	KI = ID(I) - I + K;
      	F(I) = F(I) - C * A(KI);
      end
   end
end

%---- Back substitution ----
NS = ID(NQ);
F(NQ) = F(NQ) / A(NS);
for I1 = 1:NQ-1
   I = NQ - I1;
   II = ID(I);
   C = 1 / A(II);
   F(I) = C * F(I);
   for J = I+1:NQ
      J1 = J - I + 1;
      NJ = ID(J) - ID(J - 1);
      if J1 <= NJ
         IJ = ID(J) - J + I;
      	F(I) = F(I) - C * A(IJ) * F(J);
      end
   end
end

