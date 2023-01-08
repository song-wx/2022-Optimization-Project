function [x, w, z, table, flag] = DualSimplex(A, b, c)

format rat;
table = [];
[m, n] = size(A);%初始化单纯形表
table0 = zeros(m+1, n+1);%大小为m+1 * n+1
x = zeros(1, n);%初始化最优解

%%选择A中单位阵所在列为入基列
base = eye(m);
choose = nchoosek(1:n, m);%%每次从n列里选m列判断

for i=1:size(choose,1)
   if(A(:,choose(i,:))==eye(m))
   ind = choose(i, :); %%找到一组可行解
   break
   end
end

%%若无单位阵，则找一个非退化子矩阵作为基矩阵
if(isempty(ind))
for i=1:size(choose,1)
   if(rank(A(:, choose(i,:)) ~= 0))
   ind = choose(i, :); %%找到一组可行解
   break
   end
end
end

un_ind = 1:n;
un_ind(ind) = [];%非基变量指标
B = A(:, ind); 
b_head = B\b;
table0(1:m, 1:n) = B\A;
table0(1:m, n+1) = b_head;%%第一步单纯形表
table0(end,end) = dot(b_head, c(1,ind));
cB = c(1, ind);
w = cB/B;
table0(m+1, ind) = 0;%基变量判别数为零
table0(m+1, un_ind) = w*table0(1:m, un_ind) - c(un_ind);%非基变量判别数
table = table0;

while true

if sum(table(1:m,end) < 0)==0 %%若对偶可行的基本解也是原问题的可行解，则最优
    x(1,ind) = table(1:m,end);
    z = table(end, end);
    flag = 0;

    break
end

[~, out_line] = min(table0(1:m,end));
out_ind = ind(out_line);%%确定离基变量指标

neg_ind = find(table0(out_line,1:n)<0);

if(sum(neg_ind)==0)
   flag  = 2;%%若所有yij都大于等于零，则无可行解
end

[~,in_ind] = min(table0(end,neg_ind)./table0(out_line, neg_ind));
in_ind = neg_ind(in_ind)
 
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

w = -c(ind)/A(:, ind);%%对偶问题最优解
end


