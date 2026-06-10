function [] = contour2()
clear all;global hh cc
disp('***************************************');
disp('*          PROGRAM CONTOUR2 (SHADED)  *');
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
  
 
   %============  Draw Contour  ==============
   
   for N=1 : NE
      hold on;
      if (NEN == 3)
        for IEN = 1 : NEN
         IEE = NOC(N, IEN);
         UU(IEN) = FF(IEE);
         XX(IEN) = X(IEE, 1);
         YY(IEN) = X(IEE, 2);
        end
      LPLOT(UU,XX,YY,FMIN,STP,NCL, MM);      
      elseif (NEN == 4)
           XB = 0; YB = 0; UB = 0;
           for IT = 1 : NEN
              NIT = NOC(N, IT);
              XB = XB + .25 * X(NIT, 1);
              YB = YB + .25 * X(NIT, 2);
              UB = UB + .25 * FF(NIT);
           end
           for IT = 1 : NEN
              IT1 = IT + 1;
              if (IT1 > 4)
                 IT1 = 1;
              end              
              XX(1) = XB; YY(1) = YB; UU(1) = UB;              
              NIE = NOC(N, IT);
              XX(2) = X(NIE, 1); YY(2) = X(NIE, 2); UU(2) = FF(NIE);
              NIE = NOC(N, IT1);
              XX(3) = X(NIE, 1); YY(3) = X(NIE, 2); UU(3) = FF(NIE);
              LPLOT(UU,XX,YY,FMIN,STP,NCL, MM);              
           end
        else
           disp( 'NUMBER OF ELEMENT NODES > 4 IS NOT SUPPORTED')
           stop
        end
    end
    
    %draw legend     
     dd = nonzeros(cc);
     cc1=num2str(dd,3);
     legend(hh,cc1,-1)

          
        function [] = LPLOT(UU,XX,YY,FMIN,STP,NCL, MM)
        global hh cc
      [U, INDEX] = sort(UU);     
      SU = (U(2) - FMIN) / STP;
      II = fix(SU+1.e-10)+1;
      if II>NCL
         II=NCL;
      end
      UT = FMIN + II*STP;IT=II;U2= U(2);
      ULEVEL2 = UT;ILEVEL2=II;
      
      %lower triangle
      if U(3)~=U(1)
        while UT > U(1)
          psi = (U2-U(1))/(U(3)-U(1));
          xpsi=XX(INDEX(1))+psi*(XX(INDEX(3))-XX(INDEX(1)));
          ypsi=YY(INDEX(1))+psi*(YY(INDEX(3))-YY(INDEX(1)));
          if U(2)~=U(1)
             psi2=(U2-U(1))/(U(2)-U(1));
          else
             psi2=1;
          end
          xc = XX(INDEX(1))+psi2*(XX(INDEX(2))-XX(INDEX(1)));
          yc=YY(INDEX(1))+psi2*(YY(INDEX(2))-YY(INDEX(1)));
          xe(1)=XX(INDEX(1));xe(2)=xc;xe(3)=xpsi;
          ye(1)=YY(INDEX(1));ye(2)=yc;ye(3)=ypsi;
          hh(IT,:)=fill(xe,ye,MM(IT,:),'LineWidth',.0000001,'LineStyle','none');
          cc(IT,1)=UT; 
          UT = UT-STP; U2=UT; IT = IT-1;
          if IT<1 
            IT=1;
          end    
        end
     else
          xe(1)=XX(INDEX(1));xe(2)=XX(INDEX(2));xe(3)=XX(INDEX(3));
          ye(1)=YY(INDEX(1));ye(2)=YY(INDEX(2));ye(3)=YY(INDEX(3));
			 hh(IT,:)=fill(xe,ye,MM(IT,:),'LineWidth',.0000001,'LineStyle','none');
	  end
     
     %upper triangle
     UT=ULEVEL2-STP; IT=ILEVEL2; U2=U(2);
     if U(3)~=U(1)
      while UT < U(3)
        psi = (U2-U(1))/(U(3)-U(1));
        xpsi=XX(INDEX(1))+psi*(XX(INDEX(3))-XX(INDEX(1)));
        ypsi=YY(INDEX(1))+psi*(YY(INDEX(3))-YY(INDEX(1)));
        if U(3)~=U(2)
           psi2=(U2-U(2))/(U(3)-U(2));
        else
           psi2=0;
        end
        xc = XX(INDEX(2))+psi2*(XX(INDEX(3))-XX(INDEX(2)));
        yc=YY(INDEX(2))+psi2*(YY(INDEX(3))-YY(INDEX(2)));
        xe(1)=XX(INDEX(3));xe(2)=xc;xe(3)=xpsi;
        ye(1)=YY(INDEX(3));ye(2)=yc;ye(3)=ypsi;
        hh(IT,:)=fill(xe,ye,MM(IT,:),'LineWidth',.0000001,'LineStyle','none');
        cc(IT,1)=UT; 
        UT = UT+STP; U2=UT; IT = IT+1;
        if IT>NCL 
          IT=NCL;
        end
      end
     else
      xe(1)=XX(INDEX(1));xe(2)=XX(INDEX(2));xe(3)=XX(INDEX(3));
      ye(1)=YY(INDEX(1));ye(2)=YY(INDEX(2));ye(3)=YY(INDEX(3));
	   hh(IT,:)=fill(xe,ye,MM(IT,:),'LineWidth',.0000001,'LineStyle','none');
     end


   
             
