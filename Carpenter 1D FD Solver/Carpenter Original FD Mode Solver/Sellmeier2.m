function [n_out] = Sellmeier2(lambda_in,doping)
%http://books.google.co.uk/books?id=jeWB-6K5u3EC&pg=PA104&lpg=PA104&dq=GeO2+sellmeier&source=bl&ots=K22do02tVq&sig=7MxB-dpnijGnVzhqdAR4QDNvu0s&hl=en&ei=OPjDTuuWGZTe8QPZ47SKCw&sa=X&oi=book_result&ct=result&resnum=2&ved=0CB8Q6AEwAQ#v=onepage&q=GeO2%20sellmeier&f=false
X = doping./100;

%100% SiO2
SA = [0.696166300 0.407942600 0.897479400];
SL = [0.0684043,0.1162414,9.896161];
%C1 = 0.6961663; C2 = 0.0684043; C3 = 0.4079426; C4 = 0.1162414; C5 = 0.8974794; C6 = 9.896161
%100% BK7
%SA = [1.03961212	0.231792344	1.01046945]; 
%SL = [0.00600069867	0.0200179144	103.560653];

%100% GeO2
GA = [0.80686642,0.71815848,0.85416831];
GL = [0.068972606,0.15396605,11.841931];

%lambda = (600:1800).*1e-9;
lambda = [lambda_in].*1e6;
n2 = zeros(size(lambda_in));

for i=1:length(lambda)
    n2(i) = 1+sum(((SA+X.*(GA-SA)).*lambda.^2)./(lambda.^2-(SL+X.*(GL-SL)).^2));
end

n_out = sqrt(n2);
%ratio = n(2)./n(1);
%n_out = n_in.*ratio;
end
%Optical materials By Joseph Habib Simmons, Kelly S. Potter pg. 104