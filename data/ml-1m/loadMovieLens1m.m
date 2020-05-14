function [M_train_all, M_test_all] = loadDataML1m()
try
    load ml1m.mat
catch
    downloadMovieLens1m
end
n=length(ratings(:,1));

X=ratings(:,1:2);
y=ratings(:,3);

n_movies=max(X(:,1));
n_users=max(X(:,2));
M_train_all={};
M_test_all={};

fold_struct=cvpartition(n, 'KFold', 5);

for i = 1:5
    train=find(fold_struct.training(i));
    test=find(fold_struct.training(i)==0);
    
    
    M_train=sparse(X(train,1),X(train,2), y(train),n_movies,n_users);
    M_test=sparse(X(test,1),X(test,2), y(test),n_movies,n_users);
    
    M_train_all{i}=M_train;
    M_test_all{i}=M_test;
    
end

