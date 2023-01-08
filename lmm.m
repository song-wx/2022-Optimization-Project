function [x,val,k,process]=lmm(Fk,JFk,x0)
%����: ��L-M�����������Է�����: F(x)=0
%����: x0�ǳ�ʼ��, Fk, JFk �ֱ�����F(xk)��F'(xk)�ĺ���
%���:  x, val�ֱ��ǽ��ƽ⼰||F(xk)||��ֵ,  k�ǵ�������.
maxk=1000;   %��������������
rho=0.55;sigma=0.4; muk=norm(Fk(x0));
k=0;  epsilon=1e-6; n=length(x0);
process = x0;
while(k<maxk)
    fk=Fk(x0); %���㺯��ֵ
    jfk=JFk(x0); %����Jacobi��
    gk=jfk'*fk;
    dk=-(jfk'*jfk+muk*eye(n))\gk; %�ⷽ����Gk*dk=-gk, ������������
    if(norm(gk)<epsilon), break; end  %������ֹ׼��
    m=0; mk=0;
    while(m<20)   % ��Armijo�����󲽳� 
        newf=0.5*norm(Fk(x0+rho^m*dk))^2;
        oldf=0.5*norm(Fk(x0))^2;
        if(newf<oldf+sigma*rho^m*gk'*dk)
            mk=m; break;
        end
        m=m+1;
    end
    x0=x0+rho^mk*dk;
    muk=norm(Fk(x0));
    k=k+1;
    process = [process, x0];
end
x=x0;
val=0.5*muk^2; 
%gval=norm(gfun(x));
%%%% Ŀ�꺯�� %%%%%%%%%

% function y=Fk(x)
% y(1)=x(1)-0.7*sin(x(1))-0.2*cos(x(2));
% y(2)=x(2)-0.7*cos(x(1))+0.2*sin(x(2));
% y=y(:);
% %%%% Jacobi �� %%%%%%%%%%%%%%%%%%%
% function JF=JFk(x)
%    JF=[1-0.7*cos(x(1)), 0.2*sin(x(2));
%         0.7*sin(x(1)), 1+0.2*cos(x(2))];
