function [x,val,k,process]=frcg(fun,gfun,x0)
% ����: ��FR�����ݶȷ������Լ������:  min f(x)
%����:  x0�ǳ�ʼ��, fun, gfun�ֱ���Ŀ�꺯�����ݶ�
%���:  x, val�ֱ��ǽ������ŵ������ֵ,  k�ǵ�������, process��¼��������
maxk=50000;   %����������
rho=0.6;sigma=0.4;
k=0;  epsilon=1e-10; 
n=length(x0);
process = x0;  % ��¼�������� 
while(k<maxk)
    g=feval(gfun,x0);  %�����ݶ�
    itern=k-(n+1)*floor(k/(n+1));
    itern=itern+1;  %% FR�����ݶ���n+1Ϊһ��
    %������������
    if(itern==1)  
        d=-g;  
    else
        beta=(g'*g)/(g0'*g0);
        d=-g+beta*d0;  gd=g'*d;
        if(gd>=0.0)
            d=-g;  
        end
    end
    if(norm(g)<epsilon) 
        break; 
    end   %������ֹ����

    m=0; mk=0;

    while(m<200)   %Armijo����
        if(feval(fun,x0+rho^m*d)<feval(fun,x0)+sigma*rho^m*g'*d)
            mk=m; break;
        end
        m=m+1;
    end
    x0=x0+rho^mk*d;
    val=feval(fun,x0);
    g0=g;  d0=d; 
    k=k+1;
    process = [process, x0];
end
x=x0;
val=feval(fun,x); 