clearvars
close all
clc

load('Results.mat')
PlotScript

load('Results_4250.mat')
figure(2)
% subplot(6,6,1)
% hold on
% histogram(Xrs(1,:),20,'Normalization','pdf',...
%     'DisplayStyle','stairs','EdgeColor',Colors(8,:),'LineWidth',2)

subplot(6,6,8)
hold on
histogram(Xrs(2,:),20,'Normalization','pdf',...
    'DisplayStyle','stairs','EdgeColor',Colors(8,:),'LineWidth',2)

subplot(6,6,15)
hold on
histogram(Xrs(3,:),20,'Normalization','pdf',...
    'DisplayStyle','stairs','EdgeColor',Colors(8,:),'LineWidth',2)

% subplot(6,6,22)
% hold on
% histogram(Xrs(4,:),20,'Normalization','pdf',...
%     'DisplayStyle','stairs','EdgeColor',Colors(8,:),'LineWidth',2)
% 
% subplot(6,6,29)
% hold on
% histogram(Xrs(5,:),20,'Normalization','pdf',...
%     'DisplayStyle','stairs','EdgeColor',Colors(8,:),'LineWidth',2)
% 
subplot(6,6,36)
hold on
histogram(Xrs(6,:),20,'Normalization','pdf',...
    'DisplayStyle','stairs','EdgeColor',Colors(8,:),'LineWidth',2)