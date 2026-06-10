function [] = contour1()
global hh
%clear all;
disp('***************************************');
disp('*          PROGRAM CONTOUR1           *');
disp('*    	   LINE CONTOUR PLOT	     	*');
disp('* T.R.Chandrupatla and A.D.Belegundu 	*');
disp('***************************************');      

FILE1 = input('Input Mesh Data File Name: ','s');
LINP  = fopen(FILE1,'r');
DUMMY = fgets(LINP);
TITLE = fgets(LINP);
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[NN, NE, NM, NDIM, NEN, NDN] = deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5),TMP(6));
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));

%IC=['y','m','c','r','g','b','w','k'];
%----- Coordinates -----
DUMMY = fgets(LINP);
for I=1:NN
   TMP = str2num(fgets(LINP));
   [j,X(j,1),X(j,2)] = deal(TMP(1),TMP(2), TMP(3));
end

%----- Connectivity -----
DUMMY = fgets(LINP);
for J=1:NE
   TMP = str2num(fgets(LINP));
   [I,NOC(I,:)] = ...
      deal(TMP(1),TMP(2:NEN+1));
end

%read nodal variables
disp(blanks(1));
FILE2 = input('Input File Name Containing Nodal Values: ','s');
LINP2  = fopen(FILE2,'r');
     %----- Nodal Values
DUMMY = fgets(LINP2);
    for I = 1 : NN
       FF(I) = str2num(fgets(LINP2));
    end
    
    FMAX = max(FF); FMIN = min(FF);
    
  NCL = 10; %no. of colours
  STP = (FMAX - FMIN) / NCL;
  % Red-Orange-Green-Blue-Magenta colors for 10 contour lines [MM(1,:) to MM(10,:)]
  for i=1:10
    a=((10-i)/12); b=1; c=1; 
    H=[a,b,c];MM(i,:)=hsv2rgb(H);
  end
    
%find boundary lines
     %Edges defined by nodes in NOC to nodes in NCON
    for IE = 1 : NE
       for I = 1 : NEN
          I1 = I + 1; 
          if I1 > NEN 
             I1 = 1;
          end
          NCON(IE, I) = NOC(IE, I1);
       end
    end
    for IE = 1 : NE
     for I = 1 : NEN
       I1 = NCON(IE, I); I2 = NOC(IE, I);
     	 INDX = 0;
     	 for JE = IE + 1 : NE
     		for J = 1 : NEN
           flow=2;
           if (NCON(JE, J) == 0)
                 flow=1;
           end
           if ((I1 ~= NCON(JE, J)) & (I1 ~= NOC(JE, J))) 
              flow=1;
           end
           if ((I2 ~= NCON(JE, J)) & (I2 ~= NOC(JE, J)))
              flow=1;
           end
           if (flow==2)
              NCON(JE, J) = 0; INDX = INDX + 1;
           end
           
			end
       end
       if INDX > 0
          NCON(IE, I) = 0;
       end
    end
   end
   axis off
   %============  Draw Boundary  ==============
  for IE = 1 : NE
     for I = 1 : NEN
        if NCON(IE, I) > 0
           I1 = NCON(IE, I); I2 = NOC(IE, I);
           xe(1)=X(I1,1);xe(2)=X(I2,1);
           ye(1)=X(I1,2);ye(2)=X(I2,2);
           line (xe,ye, 'LineWidth', 2, 'color', 'b')
        end
     end
  end

  
 %===========  Contour Plotting  ===========
     for IE = 1 : NE
        if (NEN == 3)
           for IEN = 1 : NEN
              IEE = NOC(IE, IEN);
              U(IEN) = FF(IEE);
              XX(IEN) = X(IEE, 1);
              YY(IEN) = X(IEE, 2);
           end
           LPLOT(U,XX,YY,FMIN,STP,NCL, MM);
        elseif (NEN == 4)
           XB = 0; YB = 0; UB = 0;
           for IT = 1 : NEN
              NIT = NOC(IE, IT);
              XB = XB + .25 * X(NIT, 1);
              YB = YB + .25 * X(NIT, 2);
              UB = UB + .25 * FF(NIT);
           end
           for IT = 1 : NEN
              IT1 = IT + 1;
              if (IT1 > 4)
                 IT1 = 1;
              end              
              XX(1) = XB; YY(1) = YB; U(1) = UB;              
              NIE = NOC(IE, IT);
              XX(2) = X(NIE, 1); YY(2) = X(NIE, 2); U(2) = FF(NIE);
              NIE = NOC(IE, IT1);
              XX(3) = X(NIE, 1); YY(3) = X(NIE, 2); U(3) = FF(NIE);
              LPLOT(U,XX,YY,FMIN,STP,NCL, MM);
           end
        else
           disp( 'NUMBER OF ELEMENT NODES > 4 IS NOT SUPPORTED')
           stop
        end
     end
     
      %draw legend
        for i=1 : NCL
           cc(i,1)=FMIN+STP*i; 
        end
        cc1=num2str(cc,3);
        legend(hh,cc1,-1)

             
function [] = LPLOT(U,XX,YY,FMIN,STP,NCL, MM)
global hh
%THREE POINTS IN ASCENDING ORDER
		     for I = 1 : 2
           C = U(I); II = I;
           for J = I + 1 : 3
              if C > U(J)
                 C = U(J); II = J;
              end
           end
           U(II) = U(I); U(I) = C;
           C1 = XX(II); XX(II) = XX(I); XX(I) = C1;
           C1 = YY(II); YY(II) = YY(I); YY(I) = C1;
        end
        
      SU = (U(1) - FMIN) / STP;
      II = fix(SU+1.e-10)+1;
      if II>NCL
         II=NCL;
      end
      
       UT = FMIN + II * STP;
        while UT <= U(3)
           X1 = ((U(3) - UT) * XX(1) + (UT - U(1)) * XX(3)) / (U(3) - U(1));
           Y1 = ((U(3) - UT) * YY(1) + (UT - U(1)) * YY(3)) / (U(3) - U(1));
           L = 1;
           if UT > U(2) 
             L = 3;
           end
           X2 = ((U(L) - UT) * XX(2) + (UT - U(2)) * XX(L)) / (U(L) - U(2));
           Y2 = ((U(L) - UT) * YY(2) + (UT - U(2)) * YY(L)) / (U(L) - U(2));
           xe(1)=X1;xe(2)=X2;ye(1)=Y1;ye(2)=Y2;           
           hh(II,1)=line(xe,ye,'LineWidth',2,'color',MM(II,:));
           UT = UT + STP;
           II = II + 1;
        end
