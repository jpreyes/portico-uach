function []=hexafron()
clear all;
close all;
disp(blanks(2));
disp('*****        PROGRAM HEXAFRON         *')
disp('*  3-D STRESS ANALYSIS USING  8-NODE  *')
disp('*   ISOPARAMETRIC HEXAHEDRAL ELEMENT  *')
disp('*         USING FRONTAL SOLVER        *')
disp('* T.R.Chandrupatla and A.D.Belegundu  *');
disp('***************************************');
disp(blanks(2));

FILE1 = input('Input Data File Name ','s');
LINP  = fopen(FILE1,'r');
FILE2 = input('Output Data File Name ','s');
LOUT  = fopen(FILE2,'w');
LSCR = fopen('SCRATCH.DAT','w+');

DUMMY = fgets(LINP);
TITLE = fgets(LINP);
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[NN, NE, NM, NDIM, NEN, NDN] = deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5),TMP(6));
NQ = NDN * NN;
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));

NPR=3;  %  E, NU, ALPHA

[X,NOC,F,MAT,DT,PM,NU,U,MPC,BT] = GETDATA(LINP,LOUT,NN,NE,NM, ...
                                NDIM,NEN,NDN,ND,NL,NPR,NMPC,NQ);
NEDF = NEN * NDN; 
[NFRON,NTOGO,NDCNT,IBL,INDX,NOC] = PREFRONT(NN,NDN,NEN,NE,NMPC,NOC,MPC);
      S=zeros(IBL,IBL); ISBL=zeros(IBL);
      ICOUNT = 0;
     %=====  FRONTAL ASSEMBLY & ELIMINATON ETC.  =====
     %----- Corner Nodes and Integration Points
[XNI,XI] = INTEG;

     MTN1 = 0;
     for N = 1:NE
        disp('... Forming Stiffness Matrix of Element '), N
        MTN = MAT(N);
        if MTN ~= MTN1 
            [MTN1,AL,D] = DMAT(MTN,PM);
        end
       
        [SE,QT] = ELSTIF(N,AL,DT,XI,XNI,X,NOC,D);
        
        if N == 1 
           CNST = 0;
           for I = 1:NEDF
               CNST = CNST + SE(I, I);
           end
           CNST = 1E+11 * CNST;
           if NMPC > 0
             [NFRON,ISBL,S,F] = MPCFRON(NMPC,NFRON,CNST,INDX,MPC,BT,S,F);
           end
        end
 %----- Account for temperature loads QT()
        for I = 1:NEN
          IL = 3 * (I - 1); IG = 3 * (abs(NOC(N, I)) - 1);
          for J = 1:3
             IL = IL + 1; IG = IG + 1;
             F(IG) = F(IG) + QT(IL);
          end
        end
        %Frontal assembly  and Forward Elimination
       [NFRON,ICOUNT,S,F,ISBL,INDX] = FRONT(N,ND,NEDF,NEN,NDN,NE,NFRON,NTOGO,CNST,INDX, ...
                                SE,S,F,U,NU,NOC,ICOUNT,LSCR,NDCNT,ISBL);

     end
     %----- Assembly and reduction are complete
     %----- Now Backsubstitute
     [F] = BACKSUB(LSCR,ICOUNT);
     fclose(LSCR);

   
disp(sprintf('Output for Input Data from file %s\n',FILE1));
fprintf(LOUT,'Output for Input Data from file %s\n',FILE1);
   
disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);

disp(' Node#    X-Displ         Y-Displ         Z-Displ');
fprintf(LOUT,' Node#    X-Displ         Y-Displ         Z-Displ\n');
% print a matrix
for I = 1:NN
  disp(sprintf(' %4d %15.4E %15.4E %15.4E\n',[I,F(3*I-2),F(3*I-1),F(3*I)]'));
  fprintf(LOUT,' %4d %15.4E %15.4E %15.4E\n',[I,F(3*I-2),F(3*I-1),F(3*I)]');
end

%----- Reaction Calculation -----
disp(sprintf('  DOF#     Reaction'));
fprintf(LOUT,'  DOF#     Reaction\n');
for I = 1:ND
   N = NU(I);
   R = CNST * (U(I) - F(N));
   disp(sprintf(' %4d %15.4E',N,R));
   fprintf(LOUT,' %4d %15.4E\n',N,R);
end

%-----  Stress Calculations -----
disp(sprintf('ELEM#	 von Mises Stresses at Integ_points'));
fprintf(LOUT,'ELEM#	 von Mises Stresses at Integ_points\n');

     MTN1 = 0;
     for N = 1:NE
        MTN = MAT(N); if MTN ~= MTN1 
            [MTN1,AL,D] = DMAT(MTN,PM);
        end
        CAL = AL * DT(N);
        for IP = 1:8
           %--- Von Mises Stress at Integration Points
           [DJ,B,DB] = DBMAT(N,IP,XI,XNI,X,NOC,D);
           %--- Element Nodal Displacements stored in QT()
           for I = 1:8
              IN = 3 * (abs(NOC(N, I)) - 1); II = 3 * (I - 1);
              for J = 1:3
                 QT(II + J) = F(IN + J);
             end
           end
           %--- Stress Calculation STR = DB * Q
           for I = 1:6
              STR(I) = 0;
              for J = 1:24
                 STR(I) = STR(I) + DB(I, J) * QT(J);
             end
             STR(I) = STR(I) - CAL * (D(I, 1) + D(I, 2) + D(I, 3));
           end
           %--- Calculation of Von Mises Stress at IP
           SIV1 = STR(1) + STR(2) + STR(3);
           SIV2 = STR(1) * STR(2) + STR(2) * STR(3) + STR(3) * STR(1);
           SIV2 = SIV2 - STR(4) ^ 2 - STR(5) ^ 2 - STR(6) ^ 2;
           VM(IP) = sqrt(SIV1 * SIV1 - 3 * SIV2);
         end
        disp(sprintf('%5d  %14.4E %14.4E %14.4E %14.4E',N,VM(1:4)));
	    fprintf(LOUT,'%5d  %14.4E %14.4E %14.4E %14.4E\n',N,VM(1:4));
        disp(sprintf('%5d  %14.4E %14.4E %14.4E %14.4E',N,VM(5:8)));
	    fprintf(LOUT,'%5d  %14.4E %14.4E %14.4E %14.4E\n',N,VM(5:8));    
    end
fclose(LOUT);
disp(blanks(1)); 
disp(sprintf('The Results are available in the text file %s', FILE2));

     
%===============  READ DATA  ====================
function [X,NOC,F,MAT,DT,PM,NU,U,MPC,BT] = GETDATA(LINP,LOUT,NN,NE,NM,NDIM, ...
                                            NEN,NDN,ND,NL,NPR,NMPC,NQ)

DUMMY = fgets(LINP);
for I=1:NN
   TMP = str2num(fgets(LINP));
   [N, X(N,:)]=deal(TMP(1),TMP(2:1+NDIM));
end
%----- Connectivity -----
DUMMY = fgets(LINP);
for I=1:NE
   TMP = str2num(fgets(LINP));
   [N,NOC(N,:), MAT(N,:), DT(N,:)] = ...
      deal(TMP(1),TMP(2:1+NEN), TMP(2+NEN), TMP(3+NEN));
end
%----- Specified Displacements -----
DUMMY = fgets(LINP);
for I=1:ND
   TMP = str2num(fgets(LINP));
   [NU(I,:),U(I,:)] = deal(TMP(1), TMP(2));
end
%----- Component Loads -----
DUMMY = fgets(LINP);
F = zeros(NQ,1);
for I=1:NL
   TMP = str2num(fgets(LINP));
   [N,F(N)]=deal(TMP(1),TMP(2));
end
%----- Material Properties -----
DUMMY = fgets(LINP);
for I=1:NM
   TMP = str2num(fgets(LINP));
   [N, PM(N,:)] = deal(TMP(1), TMP(2:NPR+1));
end
%----- Multi-point Constraints B1*Qi+B2*Qj=B0
if NMPC > 0
   DUMMY = fgets(LINP);
   for I=1:NMPC
   	TMP = str2num(fgets(LINP));
      [BT(I,1), MPC(I,1), BT(I,2), MPC(I,2), BT(I,3)] = ...
         			deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5));
   end
end
fclose(LINP);

function [XNI,XI] = INTEG();
%------- Integration Points XNI() --------
     C = .57735026919;
     XI(1, 1) = -1; XI(2, 1) = -1; XI(3, 1) = -1;
     XI(1, 2) = 1; XI(2, 2) = -1; XI(3, 2) = -1;
     XI(1, 3) = 1; XI(2, 3) = 1; XI(3, 3) = -1;
     XI(1, 4) = -1; XI(2, 4) = 1; XI(3, 4) = -1;
     XI(1, 5) = -1; XI(2, 5) = -1; XI(3, 5) = 1;
     XI(1, 6) = 1; XI(2, 6) = -1; XI(3, 6) = 1;
     XI(1, 7) = 1; XI(2, 7) = 1; XI(3, 7) = 1;
     XI(1, 8) = -1; XI(2, 8) = 1; XI(3, 8) = 1;
     for I = 1: 8
        XNI(1, I) = C * XI(1, I);
        XNI(2, I) = C * XI(2, I);
        XNI(3, I) = C * XI(3, I);
    end
    
function [MTN1,AL,D] = DMAT(MTN,PM);
     %--- D() Matrix relating Stresses to Strains
     E = PM(MTN, 1); PNU = PM(MTN, 2); AL = PM(MTN, 3);
     C1 = E / ((1 + PNU) * (1 - 2 * PNU)); C2 = .5 * E / (1 + PNU);
     D=zeros(6,6);
     D(1, 1) = C1 * (1 - PNU); D(1, 2) = C1 * PNU; D(1, 3) = D(1, 2);
     D(2, 1) = D(1, 2); D(2, 2) = D(1, 1); D(2, 3) = D(1, 2);
     D(3, 1) = D(1, 3); D(3, 2) = D(2, 3); D(3, 3) = D(1, 1);
     D(4, 4) = C2; D(5, 5) = C2; D(6, 6) = C2;
     MTN1 = MTN;

function [SE,QT] = ELSTIF(N,AL,DT,XI,XNI,X,NOC,D);
%--------  Element Stiffness  -----
     SE = zeros(24,24); QT = zeros(24,1);     
     DTE = DT(N);
     %--- Weight Factor is ONE
     %--- Loop on Integration Points
     for IP = 1: 8
        %disp('Integration Point = '), IP
        [DJ,B,DB] = DBMAT(N,IP,XI,XNI,X,NOC,D);
        %--- Element Stiffness Matrix  SE
        SE = SE + DJ*B'*D*B;
        %--- Determine Temperature Load QT()
        C = AL * DTE;
        for I = 1:24
           DSUM = DB(1, I) + DB(2, I) + DB(3, I);
           QT(I) = QT(I) + C * abs(DJ) * DSUM / 6;
        end
     end
    
   
     
function [DJ,B,DB] = DBMAT(N,IP,XI,XNI,X,NOC,D);
%-------  DB()  MATRIX  ------
     %--- Gradient of Shape Functions - The GN() Matrix
     for I = 1:3
        for J = 1:8
           C = 1;
           for K = 1: 3
              if K ~= I 
                 C = C * (1 + XI(K, J) * XNI(K, IP));
             end
           end
           GN(I, J) = .125 * XI(I, J) * C;
        end
     end
     %--- Formation of Jacobian  TJ
     
     for I = 1: 3
        for J = 1: 3
           TJ(I, J) = 0;
           for K = 1: 8
              KN = abs(NOC(N, K));
              TJ(I, J) = TJ(I, J) + GN(I, K) * X(KN, J);
           end
        end
     end
     %--- Determinant of the JACOBIAN
     DJ1 = TJ(1, 1) * (TJ(2, 2) * TJ(3, 3) - TJ(3, 2) * TJ(2, 3));
     DJ2 = TJ(1, 2) * (TJ(2, 3) * TJ(3, 1) - TJ(3, 3) * TJ(2, 1));
     DJ3 = TJ(1, 3) * (TJ(2, 1) * TJ(3, 2) - TJ(3, 1) * TJ(2, 2));
     DJ = DJ1 + DJ2 + DJ3;
     %--- Inverse of the Jacobian AJ()
     AJ(1, 1) = (TJ(2, 2) * TJ(3, 3) - TJ(2, 3) * TJ(3, 2)) / DJ;
     AJ(1, 2) = (TJ(3, 2) * TJ(1, 3) - TJ(3, 3) * TJ(1, 2)) / DJ;
     AJ(1, 3) = (TJ(1, 2) * TJ(2, 3) - TJ(1, 3) * TJ(2, 2)) / DJ;
     AJ(2, 1) = (TJ(2, 3) * TJ(3, 1) - TJ(2, 1) * TJ(3, 3)) / DJ;
     AJ(2, 2) = (TJ(1, 1) * TJ(3, 3) - TJ(1, 3) * TJ(3, 1)) / DJ;
     AJ(2, 3) = (TJ(1, 3) * TJ(2, 1) - TJ(1, 1) * TJ(2, 3)) / DJ;
     AJ(3, 1) = (TJ(2, 1) * TJ(3, 2) - TJ(2, 2) * TJ(3, 1)) / DJ;
     AJ(3, 2) = (TJ(1, 2) * TJ(3, 1) - TJ(1, 1) * TJ(3, 2)) / DJ;
     AJ(3, 3) = (TJ(1, 1) * TJ(2, 2) - TJ(1, 2) * TJ(2, 1)) / DJ;
     %--- H() Matrix relates local derivatives of  u  to local
     %    displacements  q
     H=zeros(9,24);
     for I = 1: 3
        for J = 1: 3
           IR = 3 * (I - 1) + J;
           for K = 1: 8
              IC = 3 * (K - 1) + I;
              H(IR, IC) = GN(J, K);
          end
        end
    end
     %--- G() Matrix relates strains to local derivatives of  u
     G=zeros(6,9);
     G(1, 1) = AJ(1, 1); G(1, 2) = AJ(1, 2); G(1, 3) = AJ(1, 3);
     G(2, 4) = AJ(2, 1); G(2, 5) = AJ(2, 2); G(2, 6) = AJ(2, 3);
     G(3, 7) = AJ(3, 1); G(3, 8) = AJ(3, 2); G(3, 9) = AJ(3, 3);
     G(4, 4) = AJ(3, 1); G(4, 5) = AJ(3, 2); G(4, 6) = AJ(3, 3);
     G(4, 7) = AJ(2, 1); G(4, 8) = AJ(2, 2); G(4, 9) = AJ(2, 3);
     G(5, 1) = AJ(3, 1); G(5, 2) = AJ(3, 2); G(5, 3) = AJ(3, 3);
     G(5, 7) = AJ(1, 1); G(5, 8) = AJ(1, 2); G(5, 9) = AJ(1, 3);
     G(6, 1) = AJ(2, 1); G(6, 2) = AJ(2, 2); G(6, 3) = AJ(2, 3);
     G(6, 4) = AJ(1, 1); G(6, 5) = AJ(1, 2); G(6, 6) = AJ(1, 3);
     %--- B() Matrix relates strains to  q
     B = G*H;
     %--- DB() Matrix relates stresses to  q
     DB = D*B;
     
     
function [NFRON,NTOGO,NDCNT,IBL,INDX,NOC] = PREFRONT(NN,NDN,NEN,NE,NMPC,NOC,MPC);
        %----- Mark Last Appearance of Node / Make it negative in NOC()
        % Last appearance is first appearance for reverse element order
        for I = 1:NN
           IFLG=0; 
           for J = NE:-1:1
              for K = 1:NEN
                 if I == NOC(J, K) 
                    NOC(J, K) = -I; IFLG=1;
                 end
                 if IFLG==1,break,end
              end
              if IFLG==1,break,end
           end
        end
    
        %===== Block Size Determination
        NQ = NN * NDN;
        IDE=zeros(NQ,1);
        for I = 1:NMPC
           for J = 1:2
              IDE(MPC(I, J)) = 1;
           end
        end
        IFRON = 0; 
        for I = 1:NQ
           IFRON = IFRON + IDE(I);
        end
        IBL = IFRON;
        for N = 1:NE
           INEG = 0;
           for I = 1:NEN
              I1 = NOC(N, I); IA = NDN * (abs(I1) - 1);
              for J = 1:NDN
                 IA = IA + 1;
                 if IDE(IA) == 0 
                    IFRON = IFRON + 1; IDE(IA) = 1;
                 end
              end
              if I1 < 0
                 INEG = INEG + 1;
              end
           end
           if IBL < IFRON
              IBL = IFRON;
           end
           IFRON = IFRON - NDN * INEG;
        end
        disp(blanks(1));
        disp('Block size = '),IBL
        NFRON = 0; NTOGO = 0; NDCNT = 0;
        for I = 1:IBL
            INDX(I) = I;
        end
        
        
function [NFRON,ISBL,S,F] = MPCFRON(NMPC,NFRON,CNST,INDX,MPC,BT,S,F);
        %----- Modifications for Multipoint Constraints by Penalty Method
        for I = 1:NMPC
           I1 = MPC(I, 1);
           IFL = 0;
           for J = 1: NFRON
              J1 = INDX(J);
              if I1 == ISBL(J1)
                 IFL = 1; break
             end
           end
           if IFL == 0 
              NFRON = NFRON + 1; J1 = INDX(NFRON); ISBL(J1) = I1;
           end
           I2 = MPC(I, 2);
           IFL = 0;
           for K = 1:NFRON
              K1 = INDX(K);
              if K1 == ISBL(K1) 
                 IFL = 1; break
              end
           end
           if IFL == 0 
              NFRON = NFRON + 1; K1 = INDX(NFRON); ISBL(K1) = I2;
           end

           %----- Stiffness Modification
           S(J1, J1) = S(J1, J1) + CNST * BT(I, 1) ^ 2;
           S(K1, K1) = S(K1, K1) + CNST * BT(I, 2) ^ 2;
           S(J1, K1) = S(J1, K1) + CNST * BT(I, 1) * BT(I, 2);
           S(K1, J1) = S(J1, K1);
           %----- Force Modification
           F(I1) = F(I1) + CNST * BT(I, 3) * BT(I, 1);
           F(I2) = F(I2) + CNST * BT(I, 3) * BT(I, 2);
        end
                
function [NFRON,ICOUNT,S,F,ISBL,INDX] = FRONT(N,ND,NEDF,NEN,NDN,NE,NFRON,NTOGO,CNST,INDX, ...
                                SE,S,F,U,NU,NOC,ICOUNT,LSCR,NDCNT,ISBL);
        %----- Frontal Method Assembly and Elimination
%----------------  Assembly of Element N  --------------------

        for I = 1:NEN
           I1 = NOC(N, I); IA = abs(I1); IS1 = sign(I1);
           IDF = NDN * (IA - 1); IE1 = NDN * (I - 1);
           for J = 1:NDN
              IDF = IDF + 1; IE1 = IE1 + 1; IFL = 0;
              if NFRON > NTOGO 
                  for II = NTOGO + 1:NFRON
                    IX = INDX(II);
                    if IDF == ISBL(IX)
                       IFL = 1; break;
                    end
                  end
              end
      
              if IFL == 0 
                 NFRON = NFRON + 1; II = NFRON; IX = INDX(II);
              end
              ISBL(IX) = IDF; IEBL(IE1) = IX;
              if IS1 == -1 
                 NTOGO = NTOGO + 1;
                 ITEMP = INDX(NTOGO);
                 INDX(NTOGO) = INDX(II);
                 INDX(II) = ITEMP;
              end
            end
        end
          
        for I = 1:NEDF
           I1 = IEBL(I);
           for J = 1: NEDF
              J1 = IEBL(J);
              S(I1, J1) = S(I1, J1) + SE(I, J);
           end
        end
     
%------------------------------------------------------------------
     if NDCNT < ND 
%-----  Modification for displacement BCs / Penalty Approach  -----
        for I = 1:NTOGO
           I1 = INDX(I);
           IG = ISBL(I1);
              for J = 1: ND
                 if IG == NU(J) 
                    S(I1, I1) = S(I1, I1) + CNST;
                    F(IG) = F(IG) + CNST * U(J);
                    NDCNT = NDCNT + 1;       %Counter for check
                    break;
                 end
             end
         end
     end

     %------------   Elimination of completed variables   ---------------
        NTG1 = NTOGO;
        for II = 1:NTG1
           IPV = INDX(1); IPG = ISBL(IPV);
           PIVOT = S(IPV, IPV);
        
        %-----  Write separator "0" and PIVOT value to disk  -----
           ICOUNT = ICOUNT + 1;
           IBA = 0;
           fwrite(LSCR,IBA,'int16');
           fwrite(LSCR,PIVOT,'float32');   
           
           S(IPV, IPV) = 0;
           for I = 2:NFRON
              I1 = INDX(I); IG = ISBL(I1);
              if S(I1, IPV) ~= 0 
                  C = S(I1, IPV) / PIVOT; S(I1, IPV) = 0;
                  for J = 2:NFRON
                     J1 = INDX(J);
                     if S(IPV, J1) ~= 0 
                        S(I1, J1) = S(I1, J1) - C * S(IPV, J1);
                     end
                  end
                  F(IG) = F(IG) - C * F(IPG);
              end
          end
          for J = 2:NFRON
        %-----  Write Variable# and Reduced Coeff/PIVOT to disk  -----
              J1 = INDX(J);
              if S(IPV, J1) ~= 0 
                 ICOUNT = ICOUNT + 1; 
                 IBA = ISBL(J1); CC = S(IPV,J1)/PIVOT;
                 fwrite(LSCR,IBA,'int16'); fwrite(LSCR,CC,'float32'); 
                 
                 S(IPV, J1) = 0;
              end
          end
          ICOUNT = ICOUNT + 1;
        %-----  Write Eliminated Variable# and RHS/PIVOT to disk  -----
          CC = F(IPG)/PIVOT;
          fwrite(LSCR,IPG,'int16'); fwrite(LSCR,CC,'float32');      
          F(IPG) = 0;
        %----- (NTOGO) into (1); (NFRON) into (NTOGO)
        %----- IPV into (NFRON) and reduce front & NTOGO sizes by 1
           if NTOGO > 1 
              INDX(1) = INDX(NTOGO);
          end
           INDX(NTOGO) = INDX(NFRON); INDX(NFRON) = IPV;
           NFRON = NFRON - 1; NTOGO = NTOGO - 1;
        end
       

function [F] = BACKSUB(LSCR,ICOUNT);
        %===== Backsubstitution
        IBEG = 6*(ICOUNT-1); fseek(LSCR,IBEG,'bof');
        N1=fread(LSCR,1,'int16'); F(N1)=fread(LSCR,1,'float32');
       
        while ICOUNT > 0
          IBEG = 6*(ICOUNT-1); fseek(LSCR,IBEG,'bof');
          N1=fread(LSCR,1,'int16'); F(N1)=fread(LSCR,1,'float32');
          
          ICOUNT = ICOUNT - 1;
          N2 = 1;  %dummy positive value
          while N2 > 0
            IBEG = 6*(ICOUNT-1); fseek(LSCR,IBEG,'bof');
            N2=fread(LSCR,1,'int16'); C=fread(LSCR,1,'float32');
            

            ICOUNT = ICOUNT - 1;
            if N2 ~= 0 
              F(N1) = F(N1) - C * F(N2);
            end
          end
        end