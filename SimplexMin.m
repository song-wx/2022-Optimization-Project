function [x, z, table, flag] = SimplexMin(A, b, c, ind, method)

format rat;
if(nargin == 4)
    method  = 0;%%默认原始单纯形法
end

if method==0%%原始单纯形法
[m, n] = size(A);%初始化单纯形表
table0 = zeros(m+1, n+1);%大小为m+1 * n+1
un_ind = 1:n;
un_ind(ind) = [];%非基变量指标
x = zeros(1, n);%初始化解
B = A(:, ind);
xB = B\b;
table0(1:m, 1:n) = B\A;
table0(1:m, n+1) = xB;
cB = c(ind);
table0(m+1, n+1) = cB*xB;%目标函数值
w = cB/B;
table0(m+1, ind) = 0;%基变量判别数为零
table0(m+1, un_ind) = w*table0(1:m, un_ind) - c(un_ind);%非基变量判别数

table = table0;%%储存第一个单纯形表
while true
   if sum(table0(m+1,1:end-1)>0)==0%%判别数都小于零，则为最优解
       flag = 0;
       break
   end
   [~,in_ind] = max(table0(m+1,1:end-1));%入基向量索引
   pos_ind = find(table0(1:m, in_ind) > 0); %入基列大于零的行索引
   if isempty(pos_ind)%%若入基列行均小于零，则无解
       flag = -1;
       break
   end
   if size(pos_ind,1)~=1
   [~,out_line] = min(table0(pos_ind, n+1)./table0(pos_ind,in_ind)); %计算bi/yik，确定离基向量
   out_line = pos_ind(out_line);
   else
       out_line = pos_ind;
   end
   out_ind = ind(out_line); 
   ind(out_line) = in_ind;
   un_ind(in_ind) = out_ind;
   %%更新单纯形表
   table0(out_line,:) = table0(out_line,:)/table0(out_line, in_ind);
   un_out_ind = 1:m;
   un_out_ind(out_line) = [];%%非出基向量索引
   table0(un_out_ind,:) = table0(un_out_ind,:) - table0(out_line,:).*table0(un_out_ind, in_ind);
   table0(m+1,:) = table0(m+1,:) - table0(out_line,:).*table0(m+1, in_ind);
   table = [table, table0];
   x(ind) = table(1:m,end)';
   z = table(m+1,end);
end
    x = x(1,1:n-m);
    z = table(m+1,end);
end






if method==2%%两阶段法
    [m, n] = size(A);
    e = [zeros(1,n), ones(1,m)];
    A0 = [A,eye(m)];
    
    %%第一阶段求解
    [m, n] = size(A0);%初始化单纯形表
    table0 = zeros(m+1, n+1);%大小为m+1 * n+1
    ind = n-m+1:n;%第一阶段基变量指标
    un_ind = 1:n-m;%非基变量指标
    x1 = zeros(1, n);%初始化第一阶段解
    xB = b;
    table0(1:m, 1:n) = A0;
    table0(1:m, n+1) = xB;
    cB = e(1, ind);
    table0(m+1, n+1) = sum(b);%目标函数值
    w = cB;
    table0(m+1, ind) = 0;%基变量判别数为零
    table0(m+1, un_ind) = w*table0(1:m, un_ind) - e(un_ind);%非基变量判别数

    table = table0;%%储存第一阶段第一个单纯形表

while true
   if sum(table0(m+1,1:end-1)>0.001)==0%%判别数都小于零，则为最优解
       break
   end
   [~,in_ind] = max(table0(m+1,1:end-1));%入基向量索引
   pos_ind = find(table0(1:m, in_ind) > 0); %入基列大于零的行索引
   if size(pos_ind,1)~=1
    [~,out_line] = min(table0(pos_ind, n+1)./table0(pos_ind,in_ind)); %计算bi/yik，确定离基向量
    out_line = pos_ind(out_line);
   else
       out_line = pos_ind;
   end
   out_ind = ind(out_line); 
   ind(out_line) = in_ind;
   un_ind(in_ind) = out_ind;
   %%更新单纯形表
   table0(out_line,:) = table0(out_line,:)/table0(out_line, in_ind);
   un_out_ind = 1:m;
   un_out_ind(out_line) = [];%%非出基向量索引
   table0(un_out_ind,:) = table0(un_out_ind,:) - table0(out_line,:).*table0(un_out_ind, in_ind);
   table0(m+1,:) = table0(m+1,:) - table0(out_line,:).*table0(m+1, in_ind);
   table = [table, table0];
end
   x1(ind) = table(1:m,end)';
    if sum(x1(1,n-m+1:n)~=0)~=0%%若人工变量非零，则无解
        flag = -2;
        z = 0;
        return
    end
    %%若人工变量全零，且xa分量都是非基变量
k = [];
for i = n+1:n+m
    if(sum(ind==i)==1)
        k = [k,i];%记录xa中为基变量的下标
    end
end

if isempty(k)%%求解第二阶段
   
   x = zeros(1,n-m);
   table0 = [table(:,end-n:end-m-1),table(:,end)];
   table0(end,end) =dot(x1(ind),c(ind)); 
   un_ind = 1:n-m;
   un_ind(ind) = [];
   table0(end,ind) = 0;
   table0(end,un_ind) =c(ind)*table0(1:end-1,un_ind)-c(un_ind);
   table = [table, table0];
while true 
   if sum(table0(m+1,1:end-1)>0)==0%%判别数都小于零，则为最优解
       flag = 0;
       break
   end
   [~,in_ind] = max(table0(m+1,1:end-1));%入基向量索引
   pos_ind = find(table0(1:m, in_ind) > 0); %入基列大于零的行索引
   if isempty(pos_ind)%%若入基列行均小于零，则无解
       flag = -1;
       break
   end
   if size(pos_ind,1)~=1
   [~,out_line] = min(table0(pos_ind, end)./table0(pos_ind,in_ind)); %计算bi/yik，确定离基向量
   out_line = pos_ind(out_line);
   else
       out_line = pos_ind;
   end
   out_ind = ind(out_line); 
   ind(out_line) = in_ind;
   un_ind(in_ind) = out_ind;
   %%更新单纯形表
   table0(out_line,:) = table0(out_line,:)/table0(out_line, in_ind);
   un_out_ind = 1:m;
   un_out_ind(out_line) = [];%%非出基向量索引
   table0(un_out_ind,:) = table0(un_out_ind,:) - table0(out_line,:).*table0(un_out_ind, in_ind);
   table0(m+1,:) = table0(m+1,:) - table0(out_line,:).*table0(m+1, in_ind);
   table = [table, table0];
   x(ind) = table(1:m,end)';
   z = table(m+1,end);
end

   x(ind) = table(1:m,end)';
   z = table(m+1,end);

end

end






if(method ==1)%%大M法
[m, n] = size(A);%初始化单纯形表
A = [A,eye(m)];
M = 1000000;
c = [c, M*ones(1,m)];
ind = n+1:n+m;%基变量指标
un_ind = 1:n;%非基变量指标
[m, n] = size(A);
table0 = zeros(m+1, n+1);%大小为m+1 * n+1
x = zeros(1, n);%初始化解
B = A(:, ind);
xB = B\b;
table0(1:m, 1:n) = B\A;
table0(1:m, n+1) = xB;
cB = c(ind);
table0(m+1, n+1) = cB*xB;%目标函数值
w = cB/B;
table0(m+1, ind) = 0;%基变量判别数为零
table0(m+1, un_ind) = w*table0(1:m, un_ind) - c(un_ind);%非基变量判别数

table = table0;%%储存第一个单纯形表
while true
   if sum(table0(m+1,1:end-1)>0)==0%%判别数都小于零，则为最优解
       flag = 0;
       break
   end
   [~,in_ind] = max(table0(m+1,1:end-1));%入基向量索引
   pos_ind = find(table0(1:m, in_ind) > 0); %入基列大于零的行索引
   if isempty(pos_ind)%%若入基列行均小于零
       
       if sum(x(1,m+1:end)~=0)==0 %%若人工变量全零，则无界
       flag = 1;
       break
       else 
           flag = 2;%有部分人工变量为非零，则无可行解
       end
   end
   if size(pos_ind,1)~=1
   [~,out_line] = min(table0(pos_ind, n+1)./table0(pos_ind,in_ind)); %计算bi/yik，确定离基向量
   out_line = pos_ind(out_line);
   else
       out_line = pos_ind;
   end
   out_ind = ind(out_line); 
   ind(out_line) = in_ind;
   un_ind(in_ind) = out_ind;
   %%更新单纯形表
   table0(out_line,:) = table0(out_line,:)/table0(out_line, in_ind);
   un_out_ind = 1:m;
   un_out_ind(out_line) = [];%%非出基向量索引
   table0(un_out_ind,:) = table0(un_out_ind,:) - table0(out_line,:).*table0(un_out_ind, in_ind);
   table0(m+1,:) = table0(m+1,:) - table0(out_line,:).*table0(m+1, in_ind);
   table = [table, table0];
   x(ind) = table(1:m,end)';
   z = table(m+1,end);
end

   x = x(1,1:n-m);


end

end



    