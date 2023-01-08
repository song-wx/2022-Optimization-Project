function [x,val,k,process]=grad(fun,gfun,x0)
% ����: �������½��������Լ������:  min f(x)
%����:  x0�ǳ�ʼ��, fun, gfun�ֱ���Ŀ�꺯�����ݶ�
%���:  x, val�ֱ��ǽ������ŵ������ֵ,  k�ǵ���������process�ǵ�������
maxk=50000;   %����������
rho=0.5;sigma=0.04;
k=0;  epsilon=1e-10;
process=x0;  %%��¼��������
while(k<maxk)
    g=feval(gfun,x0);  %����feval������ͨ��gfun�ݶȸ�ʽ�͵�x0�����ݶ�
    d=-g;    %������������
    if(norm(d)<epsilon), break; end
    m=0; mk=0;
    while(m<20)   %Armijo����
        if(feval(fun,x0+rho^m*d)<feval(fun,x0)+sigma*rho^m*g'*d)
            mk=m; break;
        end
        m=m+1;
    end
    x0=x0+rho^mk*d;
    process=[process,x0];
    k=k+1;
end
x=x0;
val=feval(fun,x0); 
