function [x,val,k,process]=revisenm(fun,gfun,Hess,x0)
%����: ������ţ�ٷ������Լ������:  min f(x)
%����: x0�ǳ�ʼ��, fun, gfun, Hess �ֱ�����
%         Ŀ�꺯��,�ݶ�,Hesse ��ĺ���
%���:  x, val�ֱ��ǽ������ŵ������ֵ,  k�ǵ�������,process�ǵ�������
n=length(x0); maxk=150;
rho=0.55;sigma=0.4; tau=0.0;
k=0;  epsilon=1e-10;
process = x0;% ��¼��������
while(k<maxk)
    gk=feval(gfun,x0); % �����ݶ�
    muk=norm(gk)^(1+tau);
    Gk=feval(Hess,x0);  % ����Hesse��
    Ak=Gk+muk*eye(n);
    dk=-Ak\gk; %�ⷽ����Gk*dk=-gk, ������������
    if(norm(gk)<epsilon), break; end  %������ֹ׼��
    m=0; mk=0;
    while(m<20)   %��Armijo�����󲽳� 
        if(feval(fun,x0+rho^m*dk)<feval(fun,x0)+sigma*rho^m*gk'*dk)
            mk=m; break;
        end
        m=m+1;
    end
    x0=x0+rho^mk*dk;
    k=k+1;
    process = [process,x0];
end
x=x0;
val=feval(fun,x); 
%gval=norm(gfun(x));
%%%% Ŀ�꺯�� %%%%%%%%%
%function f=fun(x)
%f=100*(x(1)^2-x(2))^2+(x(1)-1)^2;
%%%% �ݶ� %%%%%%%%%%%%%%%%%%%
%function g=gfun(x)
%g=[400*x(1)*(x(1)^2-x(2))+2*(x(1)-1), -200*(x(1)^2-x(2))]';
%%%% Hesse �� %%%%%%%%%%%%%%%%%%%
%function He=Hess(x)
%n=length(x);
%He=zeros(n,n);
%He=[1200*x(1)^2-400*x(2)+2, -400*x(1); 
 %        -400*x(1),                         200        ];
