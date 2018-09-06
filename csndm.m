function ndm = csndm(data,alpha,boxsize,normalize)
%Construction of network degree matrix
%The function performs the transformation from gene expression matrix to
%network degree matrix (ndm).
%data: Gene expression matrix, rows = genes, columns = cells
%alpha: Significant level (eg. 0.001, 0.01, 0.05 ...), Default = 0.01
%boxsize: Size of neighborhood, Default = 0.1
%normalize: 1  result is normalized (Default)
%           0  result is not normalized

if nargin < 4
    normalize = 1;
end
if nargin < 3 || isempty(boxsize)
    boxsize = 0.1;
end
if nargin <2 || isempty(alpha)
    alpha = 0.01;
end

%Define the neighborhood of each plot
[n1,n2] = size(data);
upper = zeros(n1,n2);
lower = zeros(n1,n2);
for i = 1 : n1
    [s1,s2] = sort(data(i,:));
    n0 = n2-sum(sign(s1));
    h = round(boxsize/2*sum(sign(s1)));
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
            lower(i,s2(k:k+s)) = data(i,s2(max(n0*(n0>h)+1,k-h)));
        end
        k = k+s+1;
    end
end

%If gene expression matrix is sparse, use the sparse matrix will accelerate
%the calculation and reduce memory footprint 
%data = sparse(data); upper = sparse(upper); lower = sparse(lower);

%Construction of network degree matrix
ndm = zeros(n1,n2);
B = zeros(n1,n2);
for k = 1 : n2
    for j = 1 : n2
        B(:,j) = (data(:,j) <= upper(:,k) & data(:,j) >= lower(:,k));
    end
    a = sum(B,2);
    e = sign(data(:,k));
    d = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1));
    d = (d.*(1-eye(n1)).*(e*e') > -icdf('norm',alpha,0,1));
    ndm(:,k) = sum(d,2);
    disp(['Cell ' num2str(k) ' is completed']);
end

%Normalization of network degree matrix
if normalize
    a = sum(ndm);
    a0 = mean(sum(sign(data)))^2/2000;
    for i = 1 : n1
        ndm(i,:) = ndm(i,:)./a*a0;
    end
end