clearvars
close all
clc
Color = brewermap(8,'Dark2');

%% depth of the layers (a,b,c,d)
z = [0 1.2 10 20 39];
w = diff(z)/sum(diff(z));

%% mean and std of deepest layer (d)
m = 2580;
s = 209;

%% number of samples
nos = 1e6;
%% pre-allocate array for bulk density vs depth
R = zeros(5,nos);
%% generate the samples
for kk=1:nos
    if ~mod(kk,100)
        fprintf('Sample %g/%g\r',kk,nos)
    end
    go = 1;
    while go == 1
        a = 1600+360*randn;         % sample for 1st layer
        b = 2300+130*randn;         % sample for 2nd layer
        d = 2680+(2900-2680)*rand;  % sample for 4th layer
        c = ((m+s*randn)-w(1)*a-w(2)*b-w(4)*d)/w(3);
        if b>a && c>b && d>c 
            go=0;
        end
    end
    R(:,kk) = [a b c d d];
end
mR = mean(R,2);
sR = std(R');

fprintf('Mean bulk density: %g\n',mR(3))
fprintf('Std. of bulk density: %g\n',sR(3))
fprintf('2 x std. of bulk density: %g\n',2*sR(3))

PlotInds = randi(nos,500,1);
figure, hold on
stairs(z,R(:,PlotInds),'Color',[Color(1,:) 0.1],'LineWidth',1)
stairs(z,mR,'Color',Color(2,:),'LineWidth',2)
ylim([2000 4000])
view([90 90])
set(gcf,'Color','w')
set(gca,'FontSize',16)
xlabel('Depth (km)')
ylabel('Bulk density (kg/m^3)')
box off


figure
histogram(R(3,:),'Normalization','pdf')
set(gcf,'Color','w')
set(gca,'FontSize',16)
ylabel('pdf')
xlabel('Bulk density (kg/m^3)')
box off


PlotPDF(R)