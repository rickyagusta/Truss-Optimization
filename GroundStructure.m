clear;clc;close all
iter=0
fitnessvalue=0
load data.mat
%Importing the node connection and member from the objective function
s = groundStructure{5}(:,1);
t = groundStructure{5}(:,2);
A(1:15)=14
area = A';
G = graph(s,t,area);
x_axis=groundStructure{4}(:,1);
y_axis=groundStructure{4}(:,2);

%Determining the thickness of the member in figure
LWidths = 10*G.Edges.Weight/10;
hold on

%Plotting the structure to figure
plotting = plot(G,'XData',x_axis,'YData', y_axis,'LineWidth', LWidths );

%Plotting the P (force) arrow
p1 = [x_axis(2) y_axis(2)];                         
p2 = [x_axis(2) y_axis(2)-25];                         
dp = p2-p1;                         
quiver(p1(1),p1(2),dp(1),dp(2),0, 'color',[0 0 0], 'LineWidth', 2, 'MaxHeadSize', 0.7)
text(x_axis(2), y_axis(2)-27,'P')

%Limiting the figure axis limit
xlim([-50 400])
ylim([-30 150])
load data.mat

%Looping to make the pins passive DOF
pdof=groundStructure{6,1};
for i=2:2:length(pdof)
    pdofJoint=pdof(i)/2;
    pdofLoc=groundStructure{4,1}(pdofJoint,:);
    x=pdofLoc(1);
    y=pdofLoc(2);
    fill([x x-5 x-5 x ], [y y-5 y+5 y], 'k')
end

title(sprintf('Ground structure of 15-bar truss optimization problem'))

%Set the image size
set(gcf, 'Position',  [100, 100, 900, 600])
%Saving the image to result folder turn on when needed
saveas(gcf,sprintf('result/%d.png', iter));
