%% ------------------------------------------------------------------------
figure
histogram(RMSE,50,'Normalization','pdf','FaceColor',Colors(8,:))
xlabel('RMSE')
ylabel('pdf')
set(gcf,'Color','w')
set(gca,'FontSize',30)
box off
xlim([0 3])

figure
TrianglePlot(Xrs,0.8)
Label = cell(3,1);
Label{1} = '\alpha';
Label{2} = '\phi';
Label{3} = '\gamma_w';
Label{4} = '\kappa_m [GPa]';
Label{5} = '\mu_m [GPa]';
Label{6} = '\rho_{m} [kg/m^3]';
labelTrianglePlot(Label)
set(gcf,'Position',[2 1 1317 976])
%% ------------------------------------------------------------------------


%% Plots dependent on which data we use
%% ------------------------------------------------------------------------
if sum(H)==3 %% using all three data points
    figure, hold on
    histogram(Drs(1,:),50,'Normalization','pdf','FaceColor',Colors(8,:))
    errorbar(d(1),0,2*s(1),'horizontal','.','MarkerSize',30,...
        'CapSize',18,'LineWidth',4,'Color',Colors(2,:))
    set(gca,'FontSize',30)
    box off
    xlabel('V_p [km/s]'),ylabel('pdf')
    xlim([3 5])
    set(gcf,'Color','w')

    figure, hold on 
    histogram(Drs(2,:),50,'Normalization','pdf','FaceColor',Colors(8,:))
    errorbar(d(2),0,2*s(2),'horizontal','.','MarkerSize',30,...
        'CapSize',18,'LineWidth',4,'Color',Colors(2,:))
    set(gca,'FontSize',30)
    box off
    xlim([1.5 3.5])
    xlabel('V_s [km/s]'),ylabel('pdf')
    set(gcf,'Color','w')

    figure, hold on
    histogram(Drs(3,:),50,'Normalization','pdf','FaceColor',Colors(8,:))
    errorbar(d(3),0,2*s(3),'horizontal','.','MarkerSize',30,...
        'CapSize',18,'LineWidth',4,'Color',Colors(2,:))
    set(gca,'FontSize',30)
    box off
    xlabel('\rho_b [kg/m^3]'),ylabel('pdf')
    set(gcf,'Color','w')
    % xticks([1800 2200 2600])
end

if sum(H)==2
    if H(1)==1 && H(2) == 1 %% vp and vs
        figure
        subplot(121), hold on
        histogram(Drs(1,:),50,'Normalization','pdf','FaceColor',Colors(8,:))
        errorbar(d(1),0,2*s(1),'horizontal','.','MarkerSize',30,...
            'CapSize',18,'LineWidth',4,'Color',Colors(2,:))
        set(gca,'FontSize',16)
        box off
        xlabel('vp'),ylabel('pdf')
        xlim([2 6])

        subplot(122), hold on
        histogram(Drs(2,:),50,'Normalization','pdf','FaceColor',Colors(8,:))
        errorbar(d(2),0,2*s(2),'horizontal','.','MarkerSize',30,...
            'CapSize',18,'LineWidth',4,'Color',Colors(2,:))
        set(gca,'FontSize',16)
        box off
        xlim([1 3.5])
        xlabel('vs'),ylabel('pdf')

    elseif H(1)==1 && H(3) == 1 %% vp and rho
            figure
            subplot(121), hold on
            histogram(Drs(1,:),50,'Normalization','pdf','FaceColor',Colors(8,:))
            errorbar(d(1),0,2*s(1),'horizontal','.','MarkerSize',30,...
                'CapSize',18,'LineWidth',4,'Color',Colors(2,:))
            set(gca,'FontSize',16)
            box off
            xlabel('vp'),ylabel('pdf')
            xlim([2 6])

            subplot(122), hold on
            histogram(Drs(2,:),50,'Normalization','pdf','FaceColor',Colors(8,:))
            errorbar(d(2),0,2*s(2),'horizontal','.','MarkerSize',30,...
                'CapSize',18,'LineWidth',4,'Color',Colors(2,:))
            set(gca,'FontSize',16)
            box off
            xlabel('Bulk density'),ylabel('pdf')
            set(gcf,'Color','w')

    elseif H(2)==1 && H(3) == 1 %% vs and rho
            figure
            subplot(121), hold on
            histogram(Drs(1,:),50,'Normalization','pdf','FaceColor',Colors(8,:))
            errorbar(d(1),0,2*s(1),'horizontal','.','MarkerSize',30,...
                'CapSize',18,'LineWidth',4,'Color',Colors(2,:))
            set(gca,'FontSize',16)
            box off
            xlabel('vp'),ylabel('pdf')
            xlim([2 6])

            subplot(122), hold on
            histogram(Drs(2,:),50,'Normalization','pdf','FaceColor',Colors(8,:))
            errorbar(d(2),0,2*s(2),'horizontal','.','MarkerSize',30,...
                'CapSize',18,'LineWidth',4,'Color',Colors(2,:))
            set(gca,'FontSize',16)
            box off
            xlabel('Bulk density'),ylabel('pdf')
            set(gcf,'Color','w')
    end
end


%% --------------------------------------------