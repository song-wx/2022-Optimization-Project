function [x,val,k,process]=dfp(fun,gfun,Hess,x0)
%功能: 用DFP算法求解无约束问题:  min f(x)
%输入: x0是初始点, fun, gfun分别是目标函数及其梯度
%输出:  x, val分别是近似最优点和最优值,  k是迭代次, process是迭代过程
maxk=1e5;   %给出最大迭代次数
rho=0.55;sigma=0.4; epsilon=1e-10; 
k=0;   n=length(x0); 
Hk=inv(Hess(x0));   %Hk=eye(n);
process = x0;
while(k<maxk)
    gk=feval(gfun,x0); %计算梯度
    if(norm(gk)<epsilon), break; end  %检验终止准则
    dk=-Hk*gk;  %解方程组, 计算搜索方向
    m=0; mk=0;
    while(m<20)   % 用Armijo搜索求步长 
        if(feval(fun,x0+rho^m*dk)<feval(fun,x0)+sigma*rho^m*gk'*dk)
            mk=m; break;
        end
        m=m+1;
    end
    %DFP校正
    x=x0+rho^mk*dk;  
    sk=x-x0;  yk=feval(gfun,x)-gk;
    if(sk'*yk>0)
        Hk=Hk-(Hk*yk*yk'*Hk)/(yk'*Hk*yk)+(sk*sk')/(sk'*yk);
    end
    k=k+1;     x0=x;
    process = [process, x0];
end
val=feval(fun,x0); 

