function [xk,val,k,process]=trustm(fun, gfun, Hess, x0)
%功能: 牛顿型信赖域方法求解无约束优化问题 min f(x)
%输入: x0是初始迭代点
%输出: xk是近似极小点, val是近似极小值, k是迭代次数
n=length(x0);  x=x0; dta=1;
eta1=0.15; eta2=0.75;  dtabar=2.0;  %%eta为判别逼近是否成功的参数
tau1=0.5; tau2=2.0; epsilon=1e-10;  %%tau为控制半径增减的参数
k=0;  Bk=Hess(x);  %Bk=eye(n); 
process = x;%用于记录迭代过程
while(k<150000) %% k<150000
    gk=gfun(x);   
    if(norm(gk)<epsilon)
        break;
    end
    [d,~,~,~]=trustq(gk,Bk,dta);
    deltaq=-qk(gfun,Hess,x,d);
    deltaf=fun(x)-fun(x+d);
    rk=deltaf/deltaq;
    %%更改探测半径
    if(rk<=eta1)
        dta=tau1*dta; %%逼近效果差，减小探测半径
    else 
        if (rk>=eta2&&norm(d)==dta) %%逼近效果好，增大探测半径
            dta=min(tau2*dta,dtabar);
        else
            dta=dta;%%逼近效果一般，保持探测半径
        end
    end
    %%确定步长
    if(rk>eta1)
%       x0=x;     
        x=x+d;    
%       sk=x-x0;  yk=gfun(x)-gfun(x0);  
%       vk=sqrt(yk'*Bk*yk)*(sk/(sk'*yk)-Bk*yk/(yk'*Bk*yk));
%       Bk=Bk-Bk*yk*yk'*Bk/(yk'*Bk*yk)+sk*sk'/(sk'*yk)+vk*vk';
        Bk=Hess(x);%% 计算下一步需要更新的Hess
    end
    process = [process,x];
    k=k+1;
end
xk=x;
val=fun(xk);
% %%% 目标函数  %%%%%%%%%%%%%%%
% function f=fun(x)
%  f=100*(x(1)^2-x(2))^2+(x(1)-1)^2;
%  %%% 子问题目标函数 %%%%%%%%%%%%%
function qd=qk(gfun,Hess,x,d)
gk=gfun(x);  Bk=Hess(x);
qd=gk'*d+0.5*d'*Bk*d;
%%% 目标函数的梯度 %%%%%%%%%%%%%%
% function gf=gfun(x)
% gf=[400*x(1)*(x(1)^2-x(2))+2*(x(1)-1), -200*(x(1)^2-x(2))]';
% %%% 目标函数的Hesse阵 %%%%%%%%%%%%%%
% function He=Hess(x)
% n=length(x);
% He=zeros(n,n);
% He=[1200*x(1)^2-400*x(2)+2, -400*x(1); 
%          -400*x(1),                         200        ];
