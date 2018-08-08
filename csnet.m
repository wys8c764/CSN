function csn = csnet(data,c,alpha,boxsize,weighted)
%Construction of cell-specific network
%The function performs the transformation from gene expression matrix to
%cell-specific network (csn).
%data: Gene expression matrix, rows = genes, columns = cells
%c: Construct the CSNs for all cells, set c = [] (Default);
%   Construct the CSN for cell k, set  c = k
%alpha: Significant level (eg. 0.001, 0.01, 0.05 ...)
%       larger alpha leads to more edges, Default = 0.01
%boxsize: Size of neighborhood, Default = 0.1
%weighted: 1  edge is weighted
%          0  edge is not weighted (Default)
%csn: Cell-specific network, sparse matrix, rows = genes, columns = genes
%
%Note that too many cells or genes may lead to out of memory

if nargin < 5
    weighted = 0;
end
if nargin < 4 || isempty(boxsize)
    boxsize = 0.1;
end
if nargin <3 || isempty(alpha)
    alpha = 0.01;
end
if nargin <2 || isempty(c)
    c = [];
end

%Define the neighborhood of each plot
[n1,n2] = size(data);
upper = zeros(n1,n2);
lower = zeros(n1,n2);
for i = 1 : n1
    [s1,s2] = sort(data(i,:));
    n3 = n2-sum(sign(s1));
    h = round(boxsize*sum(sign(s1)));
    k = 1;
    while k <= n2
        s = 0;
        while k+s+1 <= n2 && s1(k+s+1) == s1(k)
            s = s+1;
        end
        if s >= h
            upper(i,s2(k:k+s)) = data(i,s2(k));
            lower(i,s2(k:k+s)) = data(i,s2(k));
        else
            upper(i,s2(k:k+s)) = data(i,s2(min(n2,k+s+h)));
            lower(i,s2(k:k+s)) = data(i,s2(max(n3*(n3>h)+1,k-h)));
        end
        k = k+s+1;
    end
end

%Construction of cell-specific network
if isempty(c)
    csn = cell(1,n2);
    B = zeros(n1,n2);
    for k = 1 : n2
        for j = 1 : n2
            B(:,j) = (data(:,j) <= upper(:,k) & data(:,j) >= lower(:,k));
        end
        a = sum(B,2);
        e = sign(data(:,k));
        d = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1));
        if weighted
            csn{k} = sparse(d.*(1-eye(n1)).*(e*e').*(d > 0));
        else
            csn{k} = sparse(d.*(1-eye(n1)).*(e*e') > ...
                -icdf('norm',alpha,0,1));
        end
        disp(['Cell ' num2str(k) ' is completed']);
    end
else
    B = zeros(n1,n2);
    k = c;
    for j = 1 : n2
        B(:,j) = (data(:,j) <= upper(:,k) & data(:,j) >= lower(:,k));
    end
    a = sum(B,2);
    e = sign(data(:,k));
    d = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1));
    if weighted
        csn = sparse(d.*(1-eye(n1)).*(e*e').*(d > 0));
    else
        csn = sparse(d.*(1-eye(n1)).*(e*e') > -icdf('norm',alpha,0,1));
    end
    disp(['Cell ' num2str(k) ' is completed']);
end