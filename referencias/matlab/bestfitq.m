clear all
close all
disp('=======================================');
disp('         PROGRAM BESTFITQ              ');
disp(' T.R. Chandrupatla and A.D.Belegundu   ');
disp('=======================================');

FILE1 = input('Input Mesh Data File Name: ','s');
LINP  = fopen(FILE1,'r');disp(blanks(1));
FILE2 = input('File Name Containing Element Values: ','s');
LINP2  = fopen(FILE2,'r');disp(blanks(1));
FILE3 = input('Output Data File Name For Nodal Values:','s');
LOUT  = fopen(FILE3,'w');

%---------- Read FE data from FILE1 ----------
DUMMY = fgets(LINP);
TITLE = fgets(LINP);
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[NN, NE, NM, NDIM, NEN, NDN] = deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5),TMP(6));

%NQ = NDN * NN;
NQ = NN;

DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));

if (NDIM ~= 2)|(NEN ~= 4)
    disp('This program is for 4 noded quadrilaterals only')
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
   [N,NOC(N,:)] = deal(TMP(1),TMP(2:1+NEN));
end

fclose(LINP);

%---------- Read FE data from FILE2 ----------
DUMMY = fgets(LINP2);
for I = 1:NE
   TMP = str2num(fgets(LINP2));
   for J = 1:NEN
      [V(I, :)] = deal(TMP(1:4));
   end
end
fclose(LINP2);

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

AL = .5773502692;
%----- Shape Function Values
SH(1, 1) = .25 * (1 + AL) ^ 2;
SH(1, 2) = .25 * (1 - AL * AL);
SH(1, 3) = .25 * (1 - AL) ^ 2;
SH(1, 4) = SH(1, 2);
SH(2, 1) = SH(1, 2);
SH(2, 2) = SH(1, 1);
SH(2, 3) = SH(1, 2);
SH(2, 4) = SH(1, 3);
SH(3, 1) = SH(1, 3);
SH(3, 2) = SH(1, 2);
SH(3, 3) = SH(1, 1);
SH(3, 4) = SH(1, 2);
SH(4, 1) = SH(1, 2);
SH(4, 2) = SH(1, 3);
SH(4, 3) = SH(1, 2);
SH(4, 4) = SH(1, 1);
%----- Element Stiffness
for I = 1:4
   for J = 1:4
      C = 0;
      for K = 1:4
         C = C + SH(I, K) * SH(K, J);
      end
      SE(I, J) = C;
   end
end
%----- Global Stiffness Matrix
S = zeros(NQ,NBW);
F = zeros(NQ,1);

%----- Stiffness and Loads
for N = 1:NE
   for I = 1:4
      C = 0;
      for J = 1:4
         C = C + SH(I, J) * V(N, J);
      end
      FE(I) = C;
   end
%disp('.... Placing in Banded Locations');
   for I = 1:4
      NR = NOC(N, I);
      for J = 1:4
         NC = NOC(N, J) - NR + 1;
         if NC > 0
            S(NR, NC) = S(NR, NC) + SE(I, J);
         end
      end
      F(NR) = F(NR) + FE(I);
   end
end
F
S(:,1)
%----- Equation Solving
disp('.... Solving Equations');
[F] = bansol(NN,NBW,S,F);
F
fprintf(LOUT,'Nodal Values for Data in Files %s and %s\n', FILE1, FILE2);
for I = 1:NN
   fprintf(LOUT,'%f\n', F(I));
end
fclose(LOUT);
