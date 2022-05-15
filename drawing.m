%This function is made to plot and save the figure
function plotting = drawing(fitnessvalue, A, node_con,node_loc, iter)

close

%Importing the node connection and member from the objective function
s = node_con(:,1);
t = node_con(:,2);
area = A';
G = graph(s,t,area);
x_axis=node_loc(:,1);
y_axis=node_loc(:,2);

%Determining the thickness of the member in figure
LWidths = 10*G.Edges.Weight/10;
hold on

%Plotting the structure to figure
plotting = plot(G,'XData',x_axis,'YData', y_axis,'LineWidth', LWidths );

%Plotting the P (force) arrow
p1 = [x_axis(2) y_axis(2)];                         
p2 = [x_axis(2) y_axis(2)-30];                         
dp = p2-p1;                         
quiver(p1(1),p1(2),dp(1),dp(2),0, 'color',[0 0 0], 'LineWidth', 2, 'MaxHeadSize', 0.7)
text(x_axis(2), y_axis(2)-33,'P')

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

%Check whether there is fatal structure violation
if fitnessvalue < 10^7-1
    title({[sprintf('Iteration# %d',iter)],[sprintf('Total mass= %.1f lb',fitnessvalue)]})
else
    title({[sprintf('Iteration# %d',iter)],[sprintf('Structure is invalid')]})
end

%Set the image size
set(gcf, 'Position',  [100, 100, 900, 600])
%Saving the image to result folder turn on when needed
saveas(gcf,sprintf('result/%d.png',iter));
