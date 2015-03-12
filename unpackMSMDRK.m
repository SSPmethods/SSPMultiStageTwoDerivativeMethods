%%%Unpacking Explicit 2 Derivative Runge Kutta Method
%&& We are taking our coeficient vector X which is being fed through
%%% the optimizer and unpacking the pieces back into the original form

%A      SxS matrix and lower triangular because its explicit;
%Ahat   SxS matrix and lower triangular because its explicit;
%b      Vector of Length S
%bhat   Vector of Length S

%Count defines the number of unknows in each matrix we are building
function [A,Ahat,b,bhat] =  unpackMSMDRK(x,s)

count1= (s^2-s)/2;
count2= count1 + (s^2-s)/2;
count3= count2+s;
count4= count3 + s;

A=zeros(s,s);
count=1;
for i = 1:s-1
A(i+1,1:i)=x(count:count+i-1);
count= count+i;
end



Ahat=zeros(s,s);
count=1+count1;
for i = 1:s-1
Ahat(i+1,1:i)=x(count:count+i-1);
count= count+i;
end

b=x(count2+1:count3);
bhat=x(count3+1:count4);




end

