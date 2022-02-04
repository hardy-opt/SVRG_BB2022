function data = W8A(seed)

    M = load('w8a.mat'); 
    %data matrix M loads four files (x_train, y_train, x_test, y_test)
    
    [n,d] = size(M.x_train);
    
    
    D = M.x_train;
    
    %Data normalization (0 mean, unit varianve)
   s = std(D);
   s(s==0)=1;
   m=mean(D);
   D = (D-m)./s;
    
    D = [D  ones(n,1)];

    rng(seed);
    perm = randperm(n);
    A =  D(perm,:);%M.x_train(perm,:);
    B = M.y_train(perm);

    data.x_train = A';
    data.y_train = B';

%     f = floor(n/5);
% 
%     idx = f*(seed-1)+1:f*seed;
% 
%     idxv = ones(n,1);
%     idxv(idx) = 0;
%     data.x_val = A(idxv==0,:)';
%     data.y_val = B(idxv==0)';
%     data.x_train = A(idxv==1,:)';
%     data.y_train = B(idxv==1,:)';
    
    fprintf('This is W8A train data with n=%d, d=%d\n',size(data.x_train'));

    rng(seed);
    [n,~] = size(M.x_test);
    perm=randperm(n);
    data.x_test =  [M.x_test(perm,:) ones(n,1)]';
    data.y_test = M.y_test';


    fprintf('This is W8A test data with n=%d, d=%d\n',size(data.x_test'));

    %Initial point with different random seed
    rng(seed);
    data.w_init = randn(d+1,1);

end