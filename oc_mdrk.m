function coneq = oc_mdrk(p,x,s);
%Order Conditions Found in Tsai Chan paper for TDRK Methods,  not taking
%advantage of stage order simplifications. Currently Methods are only defined up to Fifth Order. 
%We are using the order conditions for our equality constraints in this optimization procedure.
 
%%%%Representing methods in Butcher Form
[A,Ahat,b,bhat] =  unpackMSMDRK(x,s); 
e=ones(s,1);                          
b=b(:);                             
bhat=bhat(:);
c=sum(A,2);                           %Define Abscissas for A and Ahat
chat=sum(Ahat,2);       

  %%Order Conditions;
      coneq(1) = b'*e-1;

   if p==2
      coneq(2) = b'*c+bhat'*e-1/2;
   end
      
   if p==3
      coneq(2) = b'*c+bhat'*e-1/2;
      coneq(3) = b'*c.^2+ 2*bhat'*c-1/3;
      coneq(4) = b'*(A*c)+b'*Ahat*e+bhat'*c-1/6;
      
  
   end

   if p==4
	  coneq(2) = b'*c+bhat'*e-1/2;
      coneq(3) = b'*c.^2+ 2*bhat'*c-1/3;
      coneq(4) = b'*(A*c)+b'*Ahat*e+bhat'*c-1/6;
      
      coneq(5) = b'*c.^3+3*bhat'*c.^2-1/4;
      coneq(6) = b'*(c.*((A*c)))+b'*(c.*(Ahat*e))+bhat'*c.^2+bhat'*(A*c)+bhat'*Ahat*e-1/8;
      coneq(7) = b'*A*c.^2+2*b'*(Ahat*c)+bhat'*c.^2-1/12;
      coneq(8) = b'*A^2*c+ b'*A*Ahat*e+ b'*(Ahat*c)+bhat'*(A*c)+bhat'*Ahat*e-1/24;
      
   end
   
      if p==5
	  coneq(2) = b'*c+bhat'*e-1/2;
      coneq(3) = b'*c.^2+ 2*bhat'*c-1/3;
      coneq(4) = b'*(A*c)+b'*Ahat*e+bhat'*c-1/6;
      coneq(5) = b'*c.^3+3*bhat'*c.^2-1/4;
      coneq(6) = b'*(c.*(A*c))+b'*(c.*(Ahat*e))+bhat'*c.^2+bhat'*A*c+bhat'*Ahat*e-1/8;
      coneq(7) = b'*A*c.^2+2*b'*(Ahat*c)+bhat'*c.^2-1/12;
      coneq(8) = b'*A^2*c+ b'*A*Ahat*e+ b'*(Ahat*c)+bhat'*A*c+bhat'*Ahat*e-1/24;
      
      
      coneq(9)  = b'*c.^4 + 4*bhat'*c.^3 - 1/5; 
      coneq(10) = b'*(c.^2.*(A*c)) + b'*(c.^2.*chat)+bhat'*c.^3+2*bhat'*(c.*(A*c))+2*bhat'*(c.*chat)-1/10;
      coneq(11) = b'*(c.*(A*c.^2))+2*b'*(c.*(Ahat*c))+bhat'*c.^3+bhat'*A*c.^2+2*bhat'*(Ahat*c)-1/15;
      coneq(12) = b'*(c.*(A*A*c))+b'*(c.*(A*chat))+b'*(c.*(Ahat*c))+bhat'*(c.*(A*c))+bhat'*(c.*chat)+2*bhat'*A*A*c+bhat'*A*chat+bhat'*Ahat*c-1/30;      
      coneq(13) = b'*(A*c).^2+2*b'*(chat.*(A*c))+b'*chat.^2+ 2*bhat'*(c.*(A*c))+2*bhat'*(c.*chat)-1/20;
      coneq(14) = b'*A*(c.*(A*c))+b'*A*(c.*chat)+b'*Ahat*c.^2+b'*Ahat*A*c+b'*Ahat*chat +bhat'*(c.*(A*c))+bhat'*(c.*chat)-1/40;
      coneq(15) = b'*A*A*c.^2+2*b'*A*Ahat*c+b'*Ahat*c.^2+bhat'*A*c.^2+2*bhat'*Ahat*c-1/60;
      coneq(16) = b'*A*A*A*c+b'*A*A*chat+b'*A*Ahat*c+b'*Ahat*A*c+b'*Ahat*chat+bhat'*A*A*c+bhat'*A*chat+bhat'*Ahat*c-1/120;
      end
      
   
  
 end
