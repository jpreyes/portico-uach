function geneigen()
%------------------------  function geneigen  ---------------------------
disp('****** GENERALIZED EIGENVALUE PROBLEM ******');
disp('*              Kx = lambda Mx              *');
disp('*       T.R.Chandrupatla, A.D. Belegundu   *');
disp('********************************************');

%----- Define Input and Output File Names
FILE1 = input('Input Data File Name <DOS file name>','s');
%FILE1 = 'inv_jac1.km';
LINP  = fopen(FILE1,'r');
FILE2 = input('Output Data File Name <DOS file name>','s');
%FILE2 = 'inv_jac1.gen';
LOUT  = fopen(FILE2,'w');

%----- Read in Number of Equations
DUMMY = fgets(LINP);
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[NQ, NBW] = deal(TMP(1),TMP(2));

%----- Banded Stiffness Matrix into S(NQ,NQ) -----
DUMMY = fgets(LINP);
for I = 1:NQ
   TMP = str2num(fgets(LINP));
   for JN = 1:NBW
      STIFF = TMP(JN);
      J = I + JN - 1;
      if J <= NQ
         S(I, J) = STIFF;
         S(J, I) = STIFF;
      end
   end
end
%----- Banded Mass Matrix into GM(NQ,NQ) -----
DUMMY = fgets(LINP);
for I = 1:NQ
   TMP = str2num(fgets(LINP));
   for JN = 1:NBW
      EMASS = TMP(JN);
      J = I + JN - 1;
      if J <= NQ
         GM(I, J) = EMASS;
         GM(J, I) = EMASS;
      end
   end
end

%----- Cholesky Factorization of Mass Matrix
GM = CHOLESKY(GM, NQ);
%----- Update of Stiffness Matrix - Standard Form Ax=(lambda)x
S = UPDTSTIFF(S, GM, NQ);
%----- Tri-diagonalize D() Diagonal, B() Sub-diagonal
%----- S() has the rotation matrix
[S, D, B] = TRIDIAG(S, NQ);
%----- Find Eigenvalues and Eigenvectors
[D, S, ITER] = EIGENTD(D, B, S, NQ);
%----- Determine Eigenvectors
S = EIGENVEC(S, GM, NQ);

%----- Print Eigenvalues and Eigenvectors
disp(sprintf('Eigenvalues'));
fprintf(LOUT,'Eigenvalues\n');
disp(sprintf('%17.7e', D));
fprintf(LOUT,'%17.7e', D);
fprintf(LOUT,'\n');
disp(sprintf('Eigenvectors(each column)'));
fprintf(LOUT,'Eigenvectors(each column)\n');
for I = 1:NQ
		disp(sprintf('%17.7e', S(I,:)));
		fprintf(LOUT,'%17.7e', S(I,:));
		fprintf(LOUT,'\n');
end

fclose(LINP);
fclose(LOUT);

%------------------------  function CHOLESKY  ---------------------------
function A = CHOLESKY(A, N)
%----- L into lower left triangle of A
for K = 1:N
	A(K, K) = sqrt(A(K, K));
   for I = K + 1:N
      A(I, K) = A(I, K) / A(K, K);
   end
   for J = K + 1:N
      for I = J:N
         A(I, J) = A(I, J) - A(I, K) * A(J, K);
      end
   end
end
%----- Inverse placed in upper triangle
%	    Use 1/diag for inverse diagonal
for K = 1:N - 1
   for I = 1:N - K
      IK = I + K;
      A(I, IK) = 0;
      for J = I:IK - 1
         C = A(I, J);
         if I == J; C = 1 / C ;end;
         A(I, IK) = A(I, IK) - A(IK, J) * C;
      end
      A(I, IK) = A(I, IK) / A(IK, IK);
   end
end

%------------------------  function UPDTSTIFF  ---------------------------
function A = UPDTSTIFF(A, B, N)
%
for J = N:-1:1
   for I = 1:N
      T = 0;
      for K = 1:J
         C = B(K, J);
         if K == J; C = 1 / C; end;
         T = T + C * A(I, K);
      end
      A(I, J) = T;
   end
end

for I = N:-1:1
   for J = 1:N
      T = 0;
      for K = 1:I
         C = B(K, I);
         if K == I; C = 1 / C; end
         T = T + C * A(K, J);
      end
      A(I, J) = T;
   end
end

%------------------------  function TRIDIAG  ---------------------------
function [A, D, B] = TRIDIAG(A, N)
%
for I = 1:N - 2
   AA = 0;
   for J = I + 1:N
      AA = AA + A(J, I) * A(J, I);
   end
   AA = sqrt(AA);
   WW = 2 * AA * (AA + abs(A(I + 1, I)));
   WW = sqrt(WW);
   IA = sign(A(I + 1, I));
%----- Diagonal and Next to Diagonal Term
   D(I) = A(I, I);
   B(I) = -IA * AA;
%----- Unit Vector W() in Column I from Row I+1 to N
   for J = I + 1:N
      A(J, I) = A(J, I) / WW;
   end
   A(I + 1, I) = A(I + 1, I) + IA * AA / WW;
%----- W'A in Row I from Col I+1 to N
   BET = 0;
   for J = I + 1:N
      A(I, J) = 0;
      for K = I + 1:N
         A(I, J) = A(I, J) + A(K, I) * A(K, J);
      end
      BET = BET + A(I, J) * A(J, I);
   end
%----- Modified A()                                              
   for J = I + 1:N
      for K = I + 1:N
         A(J, K) = A(J, K) - 2 * A(J, I) * A(I, K) - 2 * A(I, J) * A(K, I);
         A(J, K) = A(J, K) + 4 * BET * A(K, I) * A(J, I);
      end
   end
end

D(N - 1) = A(N - 1, N - 1);
B(N - 1) = A(N - 1, N);
D(N) = A(N, N);
A(N - 1, N - 1) = 1;
A(N, N) = 1;
A(N - 1, N) = 0;
A(N, N - 1) = 0;

%----- Now Create the Q matrix in A()
for I = 1:N - 2
   II = N - I - 1;
   A(II, II) = 1;
   for J = 1:I + 1
      IJ = II + J;
      A(II, IJ) = 0;
      for K = 1:I + 1
         IK = II + K;
         A(II, IJ) = A(II, IJ) + A(IK, IJ) * A(IK, II);
      end
      for K = 1:I + 1
         IK = II + K;
         A(IK, IJ) = A(IK, IJ) - 2 * A(IK, II) * A(II, IJ);
      end
   end
   for J = 1:I + 1
      IJ = II + J;
      A(II, IJ) = 0;
      A(IJ, II) = 0;
   end
end

%------------------------  function EIGENTD  ---------------------------
function [A, Q, ITER] = EIGENTD(A, B, Q, N)
%
ITER = 0;
M = N;
while M > 1
   ITER = ITER + 1;
   D = .5 * (A(M - 1) - A(M));
   BB = B(M - 1) * B(M - 1);
   BOT = D + sign(D) * sqrt(D * D + BB);
   P = A(1) - A(M) + BB / BOT;
   X = B(1);
   for I = 1:M - 1
      PP = sqrt(P * P + X * X);
      C = -P / PP;
      S = X / PP;
      if I > 1; B(I - 1) = C * B(I - 1) - S * X; end;
      A1 = A(I);
      A2 = A(I + 1);
      B1 = B(I);
      A(I) = A1 * C * C - 2 * B1 * C * S + A2 * S * S;
      B(I) = (A1 - A2) * C * S + B1 * (C * C - S * S);
      A(I + 1) = A1 * S * S + 2 * B1 * C * S + A2 * C * C;
%----- Update Q()
      for K = 1:N
         A1 = Q(K, I);
         A2 = Q(K, I + 1);
         Q(K, I) = C * A1 - S * A2;
         Q(K, I + 1) = S * A1 + C * A2;
      end
      if I == M - 1; break ; end
      X = -B(I + 1) * S;
      B(I + 1) = B(I + 1) * C;
      P = B(I);
   end   
   while (M > 1) & (abs(B(M - 1)) < .000001)
      M = M - 1;
   end
end

%------------------------  function EIGENVEC  ---------------------------
function A = EIGENVEC(A, B, N)
%
for I = 1:N
   for J = 1:N
      T = 0;
      for K = I:N
         C = B(I, K);
         if I == K; C = 1 / C; end;
         T = T + C * A(K, J);
      end 
      A(I, J) = T;
   end
end
