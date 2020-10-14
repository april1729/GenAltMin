function [M_train_all,M_test_all]=loadJester()
try
    raw_data=xlsread('jester-data-1.xls');
catch 
    fprintf("Downloading Jester Dataset...\n")
    URL='https://goldberg.berkeley.edu/jester-data/jester-data-1.zip';
    websave('data\jester\jester-data-1.zip', URL)
    unzip('data\jester\jester-data-1.zip', 'data\jester\jester-data-1');
    addpath('data\jester\jester-data-1')


    raw_data=xlsread('jester-data-1.xls');
    fprintf("Done \n")
end


idx=find(raw_data<=10);

fold_struct=cvpartition(length(idx), 'KFold', 5);

for i = 1:5

	M_test=zeros(size(raw_data));
	M_train=zeros(size(raw_data));

	train=(find(fold_struct.training(i)));
	test=(find(fold_struct.training(i)==0));

	M_train(idx(train))=raw_data(idx(train));
	M_test(idx(test))=raw_data(idx(test));

	M_train_all{i}=M_train;
	M_test_all{i}=M_test;

end
