function [A,b,M_train, M_test] = constructMovieLensProblem100k(fold)
URL='http://files.grouplens.org/datasets/movielens/ml-100k.zip';
try
    trainingFile=fopen("u"+fold+".base", 'r');
    testingFile=fopen("u"+fold+".test",'r');
    formatSpec="%i\t%i\t%i\t%i";
    A_test = fscanf(testingFile,formatSpec,[4, inf]);
    A_train = fscanf(trainingFile,formatSpec,[4, inf]);
catch 
    fprintf("Downloading MovieLens100k Dataset...\n");
    unzip(URL)
    fprintf("Done!\n")
    trainingFile=fopen("u"+fold+".base", 'r');
    testingFile=fopen("u"+fold+".test",'r');
    formatSpec="%i\t%i\t%i\t%i";
    A_test = fscanf(testingFile,formatSpec,[4, inf]);
    A_train = fscanf(trainingFile,formatSpec,[4, inf]);
end
I=A_train(1,:);
J=A_train(2,:);
m=max([A_train(1,:),A_test(1,:)]);
n=max([A_train(2,:),A_test(2,:)]);

M_train=sparse(A_train(1,:), A_train(2,:), A_train(3,:),m,n);
M_test=sparse(A_test(1,:), A_test(2,:), A_test(3,:),m,n);

A=sparse(1:length(I), sub2ind([m,n], I,J), ones(length(I), 1), length(I), n*m);
b=A_train(3,:)';


end

