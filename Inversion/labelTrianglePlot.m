function labelTrianglePlot(Label)
m = size(Label,1);
for kk=1:m
    subplot(m,m,(m-1)*m+kk)
    xlabel(Label{kk})

    subplot(m,m,(kk-1)*m+1)
    ylabel(Label{kk})
end