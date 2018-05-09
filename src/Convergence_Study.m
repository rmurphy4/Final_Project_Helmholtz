%Grid Convergence study
%At every node in the "N" matrix, a value of U_new is calculated and input
%into the NormU matrix. 

N=[10 20 40 60 80 100 120 140 160 200 320];

NormU=[171.86735 298.9179 508.1781 712.7300 908.8700 1.0907e+03 1.2562e+03 1.4063e+03 1.5431e+03 1.7850e+03 2.3549e+03];

    
n=length(N);
for i=2:n-1
    NormDiff(i)=NormU(i+1)-NormU(i);
end

for i=1:10
    Error(i)=abs(NormDiff(i)-NormDiff(8))/(NormDiff(8));
end

plot(N(2:10),Error(2:10),'*')
title('Grid Convergence Study')
xlabel('Nodes');
ylabel('Error')