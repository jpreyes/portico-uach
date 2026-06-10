clear all
close all
disp('=======================================');
disp('         PROGRAM BESTFIT               ');
disp('        FOR 3-NODED TRIANGLES          ');
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

if (NDIM ~= 2)|(NEN ~= 3)
    disp('This program is for 3 noded triangles  only')
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
   [N,NOC(N,:)] = ...
      deal(TMP(1),TMP(2:1+NEN));
end

fclose(LINP);

%---------- Read FE data from FILE2 ----------
DUMMY = fgets(LINP2);
for I = 1:NE
   FS(I) = str2num(fgets(LINP2));
end
fclose(LINP2);

%----- Bandwidth NBW from Connectivity NOC()
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

%----- Global Stiffness Matrix
S = zeros(NQ,NBW);
F = zeros(NQ,1);

for N = 1:NE
   %disp(sprintf('Forming Stiffness Matrix of Element %d'));
   %--- Element Stiffness Formation
   I1 = NOC(N, 1);   I2 = NOC(N, 2);   I3 = NOC(N, 3);
   X1 = X(I1, 1); 	Y1 = X(I1, 2);		
   X2 = X(I2, 1); 	Y2 = X(I2, 2);
   X3 = X(I3, 1); 	Y3 = X(I3, 2);
   X21 = X2 - X1;    X32 = X3 - X2;	 	X13 = X1 - X3;
   Y12 = Y1 - Y2;    Y23 = Y2 - Y3;		Y31 = Y3 - Y1;
   DJ = X13 * Y23 - X32 * Y31; %DETERMINANT OF JACOBIAN
   AE = abs(DJ) / 24;
   SE(1, 1) = 2 * AE; SE(1, 2) = AE; SE(1, 3) = AE;
   SE(2, 1) = AE; SE(2, 2) = 2 * AE; SE(2, 3) = AE;
   SE(3, 1) = AE; SE(3, 2) = AE; SE(3, 3) = 2 * AE;
   A1 = FS(N) * abs(DJ) / 6;
   FE(1) = A1; FE(2) = A1; FE(3) = A1;
   % disp(sprintf('Placing in Global Locations'));
   for II = 1:3
   	NR = NOC(N, II);
   	F(NR) = F(NR) + FE(II);
      for JJ = 1:3
         NC = NOC(N, JJ) - NR + 1;
      	if NC > 0
            S(NR, NC) = S(NR, NC) + SE(II, JJ);
         end
      end
   end
end

%----- Equation Solving
disp('.... Solving Equations');
[F] = bansol(NN,NBW,S,F);

%-----PRINT FINAL DENSITIES
fprintf(LOUT,'Nodal Values for Data in Files %s and %s\n', FILE1, FILE2);
for I = 1:NN
   fprintf(LOUT,'%f\n', F(I));
end
fclose(LOUT);
