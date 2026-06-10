function []=invitr()
clear all
close all
%------------------------ INVITR2  ---------------------------
disp('*****      PROGRAM INVITR2        *****');
disp('*      Inverse Iteration Method       *');
disp('*  for eigenvalues and eigenvectors   *');
disp('*        Searching in Subspace        *');
disp('*         for Banded Matrices         *');
disp('* T.R.Chandrupatla and A.D.Belegundu  *');
disp('***************************************');

InputData;
BanSolve1; 		%<----Stiffness to Upper Triangle
InverseIter;
Output;

%------------------------  function InputData  ---------------------------
function [] = InputData();
global NQ NBW S GM F ST
global TOL SH NEV
global TITLE FILE1 FILE2
global LINP LOUT

FILE1 = input('Input Data File Name <DOS file name> ','s');
%FILE1 = 'inv_jac1.km';
LINP  = fopen(FILE1,'r');
FILE2 = input('Output Data File Name <DOS file name> ','s');
%FILE2 = 'inv_jac1.inv';
LOUT  = fopen(FILE2,'w');

%----- Tolerance
TMP = input(' Tolerance (Default TOL=1E-6) ','s');
if isempty(TMP)
	TOL = .000001;
else
   TOL = str2num(TMP);
end

%----- Read Parameters
DUMMY = fgets(LINP);
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[NQ, NBW] = deal(TMP(1),TMP(2));

%----- Number of eigenvalues
TMP = input(' Number of Eigenvalues Desired (Default NEV=1) ','s');
if isempty(TMP)
	NEV = 1;
else
   NEV = str2num(TMP);
end

%--- Read in Stiffness Matrix ---
DUMMY = fgets(LINP);
for I=1:NQ
   S(I, :) = str2num(fgets(LINP));
end
%--- Read in Mass Matrix ---
DUMMY = fgets(LINP);
for I=1:NQ
   GM(I, :) = str2num(fgets(LINP));
end

%--- Shifting value ---
TMP = input(' Shifting value (Defalut SH=0) ','s');
if isempty(TMP)
   SH = 0;
else
   SH = str2num(TMP);
   S = S - SH * GM;
end

%--- Starting Eigenvector ---
TMP = input(' Starting Vector (Default ST=ones(NQ,1))');
if isempty(TMP)
   for I=1:NQ
      ST(I) = 1.;
   end
else
   for I=1:NQ
      ST(I) = input(sprintf('ST(%d) = ',I));
   end
end
disp(' ');
fclose(LINP);

%------------------------  function BandSolve1  ---------------------------
function []=BanSolve1();
global NQ NBW S GM F ST
global TOL SH NEV
global TITLE FILE1 FILE2
global LINP LOUT
%---- Gauss Elimination LDU Approach ----
%--- REDUCTION TO UPPER TRIANGULAR MATRIX ---
for K = 1:NQ - 1
   NK = NQ - K + 1;
   if NK > NBW; NK = NBW; end;
   for I = 2:NK
      C1 = S(K, I) / S(K, 1);
      I1 = K + I - 1;
      for J = I:NK
         J1 = J - I + 1;
         S(I1, J1) = S(I1, J1) - C1 * S(K, J);
      end
   end
end
%------------------------  function InverseIter  ---------------------------
function []=InverseIter();
global NQ NBW S GM F ST
global TOL SH NEV NEV1 NITER
global EVL EVC
global TITLE FILE1 FILE2
global LINP LOUT

ITMAX = 50;
NEV1 = NEV;
for NV = 1:NEV
%--- Starting Value for Eigenvector
   for I = 1:NQ
      EV1(I) = ST(I);
   end
   ITER = 0;
   EL2 = TOL/10;
   EL1 = 0;
   while abs(EL2 - EL1)/abs(EL2) > TOL
   	EL1 = EL2;
      ITER = ITER + 1;
      if ITER > ITMAX
         disp(sprintf('No Convergence for Eigenvalue# %d ', NV));
         NEV1 = NV - 1;
         return;
      end
      if NV > 1
      	%----  Starting Vector Orthogonal to Evaluated Vectors
         for I = 1 : NV - 1
            CV = 0;
            for K = 1 : NQ
               KA = K - NBW + 1;
               KZ = K + NBW - 1;
               if KA < 1;KA = 1; end;
               if KZ > NQ; KZ = NQ; end;
               for L = KA : KZ
                  if L < K 
                     K1 = L;
                     L1 = K - L + 1;
                  else
                     K1 = K;
                     L1 = L - K + 1;
                  end
                  CV = CV + EVS(K) * GM(K1, L1) * EVC(L, I);
               end
            end
            for K = 1:NQ
                EV1(K) = EV1(K) - CV * EVC(K, I);
            end
         end
      end
      for I = 1 : NQ
         IA = I - NBW + 1;
         IZ = I + NBW - 1;
         EVT(I) = 0;
         if IA < 1; IA = 1; end;
         if IZ > NQ; IZ = NQ; end;
         for K = IA : IZ
            if K < I 
               I1 = K;
               K1 = I - K + 1;
            else
               I1 = I;
               K1 = K - I + 1;
            end
            EVT(I) = EVT(I) + GM(I1, K1) * EV1(K);
         end
         EV2(I) = EVT(I);
      end
      %----- Reduction of the right hand side
     	for K = 1 : NQ - 1
        	NK = NQ - K + 1;
        	if NK > NBW; NK = NBW; end
         for I = 2:NK
            I1 = K + I - 1;
           	C1 = 1 / S(K, 1);
           	EV2(I1) = EV2(I1) - C1 * S(K, I) * EV2(K);
         end
     	end  
     	%----- Back Substitution
     	EV2(NQ) = EV2(NQ) / S(NQ, 1);
     	for II = 1 : NQ - 1
      	I = NQ - II;
        	C1 = 1 / S(I, 1);
        	NI = NQ - I + 1;
        	if NI > NBW; NI = NBW; end
        	EV2(I) = C1 * EV2(I);
        	for K = 2:NI
         	EV2(I) = EV2(I) - C1 * S(I, K) * EV2(I + K - 1);
        	end
     	end
     	%----- Reduce Right Side and Solve
     	C1 = 0;
     	C2 = 0;
     	for I = 1:NQ
      	C1 = C1 + EV2(I) * EVT(I);
     	end
     	for I = 1:NQ
        	IA = I - NBW + 1;
        	IZ = I + NBW - 1;
        	EVT(I) = 0;
        	if IA < 1; IA = 1; end;
        	if IZ > NQ; IZ = NQ; end;
        	for K = IA:IZ
           	if K < I
               I1 = K;
               K1 = I - K + 1;
           	else
              	I1 = I;
              	K1 = K - I + 1;
           	end
           	EVT(I) = EVT(I) + GM(I1, K1) * EV2(K);
        	end
    	end
      for I = 1:NQ
         C2 = C2 + EV2(I) * EVT(I);
      end
      EL2 = C1 / C2;
      C2 = sqrt(C2);
      for I = 1:NQ
         EV1(I) = EV2(I) / C2;
         EVS(I) = EV1(I);
      end
   end 
   for I = 1:NQ
      EVC(I, NV) = EV1(I);
   end
   NITER(NV) = ITER;
   EL2 = EL2 + SH;
	EVL(NV) = EL2;
end

%------------------------  function Output  ---------------------------
function []=Output();
global NQ NBW S GM F ST
global TOL SH NEV NEV1 NITER
global EVL EVC
global TITLE FILE1 FILE2
global LINP LOUT

%----- Print Eigenvalues and Eigenvectors
disp(sprintf('Eigenvalue Number'));
fprintf(LOUT,'Eigenvalue Number\n');
disp(sprintf('%17d', [1:NEV1]));
fprintf(LOUT,'%17d', [1:NEV1]);
fprintf(LOUT,'\n');
disp(sprintf('Iteration Number'));
fprintf(LOUT,'Iteration Number\n');
disp(sprintf('%17d', NITER));
fprintf(LOUT,'%17d', NITER);
fprintf(LOUT,'\n');
disp(sprintf('Eigenvalues'));
fprintf(LOUT,'Eigenvalues\n');
disp(sprintf('%17.7e', EVL));
fprintf(LOUT,'%17.7e', EVL);
fprintf(LOUT,'\n');
disp(sprintf('Eigenvectors(each column)'));
fprintf(LOUT,'Eigenvectors(each column)\n');
for I = 1:NQ
		disp(sprintf('%17.7e', EVC(I,:)));
		fprintf(LOUT,'%17.7e', EVC(I,:));
		fprintf(LOUT,'\n');
end

fclose(LOUT);
