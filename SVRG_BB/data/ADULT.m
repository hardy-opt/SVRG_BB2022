function data = ADULT(seed)

    M = load('adult.mat'); 
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
%     

     data.x_train = A';
     data.y_train = B';
    fprintf('This is Adult train data with n=%d, d=%d\n',size(data.x_train'));
    
    
    
   P = M.x_test;
   s = std(P);
   s(s==0)=1;
   m=mean(P);
   P = (P-m)./s;
   [e,~] = size(M.x_test);
   P = [P  ones(e,1)];
    
   rng(seed);
   per=randperm(e);
   data.x_test =  P(per,:)';
   data.y_test = M.y_test(per)';


    fprintf('This is Adult test data with n=%d, d=%d\n',size(data.x_test'));

    %Initial point with different random seed
    rng(seed);
    data.w_init = randn(d+1,1);

end