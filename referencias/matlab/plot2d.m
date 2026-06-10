function [] = plot2d()
clear all;
disp('***************************************');
disp('*          PROGRAM PLOT2D2            *');
disp('*        TWO DIMENSIONAL PLOT			*');
disp('* T.R.Chandrupatla and A.D.Belegundu 	*');
disp('***************************************');      

FILE1 = input('Input Data File Name <DOS file name>  ','s');
disp(blanks(1));
LINP  = fopen(FILE1,'r');
DUMMY = fgets(LINP);
TITLE = fgets(LINP);
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[NN, NE, NM, NDIM, NEN, NDN] = deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5),TMP(6));
DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NMPC]= deal(TMP(1),TMP(2),TMP(3));

%----- Coordinates -----
DUMMY = fgets(LINP);
for I=1:NN
   TMP = str2num(fgets(LINP));
   [j,X(j,1),X(j,2)] = deal(TMP(1),TMP(2), TMP(3));
end
%calculate offset, epsilon, for printing node numbers
LX = abs(max(X(:,1))-min(X(:,1)));
LY = abs(max(X(:,2))-min(X(:,2)));
LMAX = max (LX,LY);
epsilon = 0.01*LMAX;

%----- Connectivity -----
DUMMY = fgets(LINP);

for J=1:NE
   TMP = str2num(fgets(LINP));
   [I,NOC(I,:)] = ...
      deal(TMP(1),TMP(2:NEN+1));
end
close;axis off;
for N=1 : NE
  for j=1 : NEN
    xe(j)=X(NOC(N,j),1);
    ye(j)=X(NOC(N,j),2);
  end
  xe(j+1)=xe(1);ye(j+1)=ye(1);
  line(xe,ye)
end
disp('  Type 1 for on/off node labels,?n');
disp('  Type 2 for on/off element labels,?n');
disp('  Type 3 for on/off node & element labels,?n');
disp('  Type 0 to EXIT,?n');

status = 1; k=0; kold = -1;
while status ~= 0
     user_input=input(' Type <1,2,3 OR 0 > : ','s');   
     if k == 1     %on-->off
        if kold ~= -1
          ktemp=str2num(user_input);
          if ktemp ~= kold 
            user_input = num2str(kold);
          end
         end
        k=0;
     else			 %off-->on
        k=1;
     end
     figure(1);
   status = str2num(user_input);
   switch status
     case 1 
      if k==1 
        for i=1:NN
           xc = X(i,1)-epsilon;yc=X(i,2)-epsilon;
           hh1(i,:)=text(xc,yc,num2str(i),'FontSize',8);
        end 
        kold = str2num(user_input);        
      else
        delete(hh1);  
      end
     case 2
      if k==1
       	for i=1:NE
           xc=0;yc=0;
           for ii=1 : NEN
              xc=xc + X(NOC(i,ii),1);
              yc=yc + X(NOC(i,ii),2);
           end
           xc=xc/NEN;yc=yc/NEN;
		  	  hh2(i,:) = text(xc,yc,num2str(i),'FontSize',8);
          end
          kold = str2num(user_input);        
       else
          delete (hh2);
       end    
     case 3
      if k==1         
        for i=1:NN
         xc = X(i,1);yc=X(i,2);
        	hh3(i,:) = text(xc,yc,num2str(i),'FontSize',8);
        end
        for i=1:NE
           xc=0;yc=0;
           for ii=1 : NEN
              xc=xc + X(NOC(i,ii),1);
              yc=yc + X(NOC(i,ii),2);
           end
           xc=xc/NEN;yc=yc/NEN;
		  	  hh4(i,:) = text(xc,yc,num2str(i),'FontSize',8);
        end
        kold = str2num(user_input);        
      else
        delete(hh3);delete(hh4);
      end
    end
end

        