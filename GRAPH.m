% load('brainNetSet_PSRhub_show_0_1.mat')
clear;clc;
% x=abs(brainNetSet{37}(:,:,1));
load('g-12.mat')
x=ans; %abs(importdata('a.edge'));
%Create custom node labels

myLabel=importdata('AALatlas90.txt');
myLabel1 = cell(length(x));
 for i = 1:length(x)
myLabel1(i,1)=myLabel(i,1);
end
%load('node_90.mat')

%for i = 1:length(x)
%  myLabel1{i} = myLabel(1,i);
%end


%Create custom colormap
figure;
myColorMap = lines(length(x));

circularGraph(x,'Colormap',myColorMap,'Label',myLabel1);
%set(gcf,'CloseRequestFcn','closereq')