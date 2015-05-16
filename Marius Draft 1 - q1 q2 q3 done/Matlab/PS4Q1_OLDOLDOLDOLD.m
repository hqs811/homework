%% Problemset 4, Question 1
cd('/users/mariusring/Documents/Northwestern/520 FINC Todorov/PS 4/Matlab')
clear all;
clc;

%% Q1 (c)
phi= [0,0,0,0]';
omega=[1,1,0.01,0.01]';
alpha=[0,0,0.09,0.09]';
beta=[0,0,0.9,0.9]';
T=[100,1000,100,1000]';
Rejections=zeros(4,1)'

Reps=1000

phi_hat=zeros(4,Reps);
tstat=zeros(4,Reps);
D=0;

for i=1:4
    for rep=1:Reps
        %display([i, rep])
        clear z sigma2 epsilon y var 

        sigma2=zeros(D+T(i),1);
        epsilon=zeros(D+T(i),1);
        y=zeros(D+T(i),1);

        z(1:D+T(i))=normrnd(zeros(D+T(i),1),diag(ones(D+T(i))),D+T(i),1);

        % initial value for sigma and y
        sigma2(1)=omega(i);
        epsilon(1)=z(1)*sqrt(sigma2(1));
        y(1)=0;

        for t=2:D+T(i)
            
            sigma2(t)=omega(i)+alpha(i)*epsilon(t-1)^2+beta(i)*sigma2(t-1);
            
            epsilon(t)=z(t)*sqrt(sigma2(t));
            
            y(t)=phi(i)*y(t-1)+epsilon(t);
        end
        y_=y(D+2:D+T(i));
        x_=y(D+1:D+T(i)-1);
        phi_hat(i,rep)=inv(x_'*x_)*x_'*y_;
        var=x_'*x_/T(i);
        tstat(i,rep)=sqrt(T(i))*phi_hat(i,rep)/sqrt(var);
        %display([mean((y_-phi_hat(i,rep)*x_).*x_)])
       
    end
        
    Critical=1.96;
    Rejections(i)=length(tstat(i,abs(tstat(i,:))>Critical))/Reps;
end
