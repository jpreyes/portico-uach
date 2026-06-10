function []=jacobi()
clear all
close all
%------------------------ JACOBI2  ---------------------------
disp(blanks(2));
disp('*****      PROGRAM JACOBI2        *****');
disp('*      Generalised Jacobi Method      *');
disp('*  for eigenvalues and eigenvectors   *');
disp('*        for Symmetric Matrices       *');
disp('* T.R.Chandrupatla and A.D.Belegundu  *');
disp('***************************************');
disp(blanks(2));
FILE1 = input('Input Data File Name ','s');
LINP  = fopen(FILE1,'r');
FILE2 = input('Output Data File Name ','s');
LOUT  = fopen(FILE2,'w');

PI = 3.1415927;
%----- Tolerance
TMP = input(' Tolerance (Default TOL=1E-6) ','s');
if isempty(TMP)
	TOL = .000001;
else
   TOL = str2num(TMP);
end
TMP = input(' Maximum Number of Sweeps < default = 50 > ','s');
if isempty(TMP)
	NSWMAX = 50;
else
   NSWMAX = str2num(TMP);
end

DUMMY = fgets(LINP);
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[NQ, NBW] = deal(TMP(1),TMP(2));
%--- Read in Banded Stiffness Matrix ---
DUMMY = fgets(LINP);
for I=1:NQ
   S1(I, :) = str2num(fgets(LINP));
end
%--- Read in Banded Mass Matrix ---
DUMMY = fgets(LINP);
for I=1:NQ
   GM1(I, :) = str2num(fgets(LINP));
end

% -- Banded Stiffness Matrix into Square, Symmetric S(NQ,NQ) -----
S=zeros(NQ,NQ);
for I = 1: NQ; 
   for JN = 1:NBW;
     STIFF = S1(I,JN); J = I + JN - 1;
        if J <= NQ 
          S(I, J) = STIFF; S(J, I) = STIFF;
        end;
   end;
end;
GM=zeros(NQ,NQ);
for I = 1: NQ; 
   for JN = 1:NBW;
     EMASS = GM1(I,JN); J = I + JN - 1;
        if J <= NQ 
          GM(I, J) = EMASS; GM(J, I) = EMASS;
        end;
   end;
end;
disp(' ');
fclose(LINP);

% -- Initialize Eigenvector Matrix -----
     EVC=eye(NQ); 
% -- Define Tolerances     
     diag(S)
     C1=min(diag(S));C2=max(diag(GM));
     TOLS = TOL * C1; TOLM = TOL * C2; 
     K1 = 1; I1 = 1; NSW = 0;
     IFL=1;
     while IFL==1
     NSW = NSW + 1;
     if NSW > NSWMAX
        disp(sprintf('No convergence in specified # of sweeps'));
        fprintf(LOUT,'No convergence in specified # of sweeps\n');
	     fclose(LOUT);
	  end
     disp(' ----  SWEEP NUMBER  '); NSW
     for K = K1: NQ - 1
     for I = I1: K
     J = NQ - K + I;
     if (abs(S(I, J)) > TOLS | abs(GM(I, J)) > TOLM) 
       AA = S(I, I) * GM(I, J) - GM(I, I) * S(I, J);
       BB = S(J, J) * GM(I, J) - GM(J, J) * S(I, J);
       CC = S(I, I) * GM(J, J) - GM(I, I) * S(J, J);
     CAB = .25 * CC * CC + AA * BB;
     if CAB < 0 
         disp( 'Square Root of Negative Term -- Check Matrices')
         fprintf(LOUT, 'Square Root of Negative Term -- Check Matrices\n')
         fclose(LOUT);
     end
     if AA == 0 
        BET = 0; ALP = -S(I, J) / S(I, I);
     elseif BB == 0 
        ALP = 0; BET = -S(I, J) / S(J, J);
     else
        SQC = sqrt(CAB); 
        if CC < 0 
           SQC = -SQC;
        end 
        ALP = (-.5 * CC + SQC) / AA;
        BET = -AA * ALP / BB;
     end
     % -- Only Upper Triangular Part is used in Diagonalization
     if I > 1 
       for N = 1: I - 1
         SI = S(N, I); SJ = S(N, J); EMI = GM(N, I); EMJ = GM(N, J);
         S(N, I) = SI + BET * SJ; S(N, J) = SJ + ALP * SI;
         GM(N, I) = EMI + BET * EMJ; GM(N, J) = EMJ + ALP * EMI;
       end
     end
     if J < NQ 
       for N = J + 1: NQ
       SI = S(I, N); SJ = S(J, N); EMI = GM(I, N); EMJ = GM(J, N);
       S(I, N) = SI + BET * SJ; S(J, N) = SJ + ALP * SI;
       GM(I, N) = EMI + BET * EMJ; GM(J, N) = EMJ + ALP * EMI;
       end
     end
     if I < J 
       for N = I + 1: J - 1
       SI = S(I, N); SJ = S(N, J); EMI = GM(I, N); EMJ = GM(N, J);
       S(I, N) = SI + BET * SJ; S(N, J) = SJ + ALP * SI;
       GM(I, N) = EMI + BET * EMJ; GM(N, J) = EMJ + ALP * EMI;
       end
     end
     SII = S(I, I); SIJ = S(I, J); SJJ = S(J, J);
     S(I, J) = 0; S(I, I) = SII + 2 * BET * SIJ + BET * BET * SJJ;
     S(J, J) = SJJ + 2 * ALP * SIJ + ALP * ALP * SII;
     EII = GM(I, I); EIJ = GM(I, J); EJJ = GM(J, J);
     GM(I, J) = 0; GM(I, I) = EII + 2 * BET * EIJ + BET * BET * EJJ;
     GM(J, J) = EJJ + 2 * ALP * EIJ + ALP * ALP * EII;
     % *** EIGENVECTORS ***
     for N = 1: NQ
       EVI = EVC(N, I); EVJ = EVC(N, J);
       EVC(N, I) = EVI + BET * EVJ; EVC(N, J) = EVJ + ALP * EVI;
     end
     end   %endif
     end  %NEXT I
     end  %NEXT K
     for K = 1: NQ - 1
     for I = 1: K
     J = NQ - K + I;
     IFL = 0;
     if (abs(S(I, J)) > TOLS | abs(GM(I, J)) > TOLM )
        K1 = K; I1 = I; IFL = 1;
     end
     if IFL == 1 
        break;
     end
     end  %NEXT I
     if IFL == 1 
        break;
     end
     end  %NEXT K
     end  %end while loop
     %-----  Calculation of Eigenvalues -----
     for I = 1:NQ
        if abs(GM(I, I)) < TOLM 
           GM(I, I) = TOLM;
        end
        EVL(I) = S(I, I) / GM(I, I);
     end
     %----- Scaling of Eigenvectors
     for I = 1:NQ
        GM2 = sqrt(abs(GM(I, I)));
        for J = 1:NQ
           EVC(J, I) = EVC(J, I) / GM2;
        end
     end
     %-----   RESULTS   -----
     %--- Ascending Order of Eigenvalues
     [EVL,NORD]=sort(EVL);
     EVC=EVC(:,NORD);
     
     disp('Eigenvalues and Eigenvectors for Data in File:'); FILE1
     fprintf(LOUT,'Eigenvalues and Eigenvectors for Data in File\n'); FILE1
     
%----- Print Eigenvalues and Eigenvectors
disp(sprintf('Eigenvalue Number'));
fprintf(LOUT,'Eigenvalue Number\n');
disp(sprintf('%17d', [1:NQ]));
fprintf(LOUT,'%17d', [1:NQ]);
fprintf(LOUT,'\n');
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
