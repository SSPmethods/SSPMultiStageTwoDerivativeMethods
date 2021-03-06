function [rvect,Krange] = revealR(A,Ahat,b,bhat,Kmin, Kmax);
%Finds the SSP coeficient of a given method dependent on K.  
%Where K=(dtVV/dtFE).  We use this function to help build our SSP curve of
%a given method.

if nargin<6;
    Kmax=Kmin;
end

Krange=linspace(Kmin,Kmax);
rvect=0*Krange;
for i =1:length(Krange)
	K=Krange(i);
	r=0;
	s=length(A);z=zeros(s+1,1);I=eye(s+1);
	S=[[A;b],z];
	Shat=[[Ahat;bhat],z];
	e=ones(s+1,1);
	xx=0;
	rincrement=.1;
	while rincrement>10^-14;                  %Refines r(K) to the 14th decimal place

              while min(min(xx))==0;          %Keep Running until negative Coeficient appears
                    r=r+rincrement;        
                    r2=(r)^2/K^2;
                    R=inv((I+r*S+r2*Shat));   %Converting butcher to MSO
                    P=R*r*S;
                    Q=R*r2*Shat;
                    xx=[R*e P Q];
              end
        r=r-rincrement;                       %Backs up r to before negativity is introduced
        rincrement=.1*rincrement;             %shrinks increment by a power of 10
        xx=0;                                 %Allows us to re enter internal while loop
    end
    
    rvect(i) = r;                             %Stores are values of r for a given K.

end

if nargin==6;
    plot(Krange,rvect,'-bo');                 %If given a range of Ks to check, plot the generated curve
end

end
