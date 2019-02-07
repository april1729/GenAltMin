function [rankTable, errorTable, runTimeTable]=displayResults(indepVar, rank, error,runtime, labels)


rankTable = array2table([indepVar', rank], 'VariableNames', labels);
errorTable = array2table([indepVar', error], 'VariableNames', labels);

plotLabels=strrep(labels,"_"," ");

disp("Rank:");
disp(rankTable);

disp("Error:");
disp(errorTable);

if ~isempty(runtime)
    runTimeTable = array2table([indepVar', runtime], 'VariableNames', labels);
    disp("Run Time:");
    disp(runTimeTable);
    
figure('DefaultAxesFontSize',14)
subplot(3,1,1)
    plot(indepVar, rank)
    legend(plotLabels{2:length(labels)}, 'Location','EastOutside')
    ylabel("Rank \n(0.001 tolerance)")
    
    subplot(3,1,2)
    plot(indepVar, error)
    ylabel("Relative \n Error")
    legend(plotLabels{2:length(labels)}, 'Location','EastOutside')
    xlabel(plotLabels{1})
    
    subplot(3,1,3)
    plot(indepVar, runtime)
    ylabel("Run Time \n (Seconds)")
    legend(plotLabels{2:length(labels)}, 'Location','EastOutside')
    xlabel(plotLabels{1})
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [.2, 0.2, .8, 0.8])
    
else
    runTimeTable=[]
    
figure('DefaultAxesFontSize',18) 
subplot(2,1,1)
    plot(indepVar, rank)
    legend(plotLabels{2:length(labels)}, 'Location','EastOutside')
    ylabel("Rank \n (0.001 tolerance)")
    
    subplot(2,1,2)
    plot(indepVar, error)
    ylabel("Relative \n Error")
    legend(plotLabels{2:length(labels)}, 'Location','EastOutside')
    xlabel(plotLabels{1})
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [.2, 0.2, .8, 0.8])
end
