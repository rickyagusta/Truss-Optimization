%If caseNum = 1 then it's 15-bar truss case
function [ub, lb]= boundary(caseNum)

%Try not to load data.mat multiple times during the code
%Load groundStructure matrix to know the ub and lb 
load data.mat

%Bar available on ground structure
barNumber = groundStructure{caseNum};

%Ub and lb for size and topology problems
%{
For the 15 bar truss, the first 15 cell (1-15) of ub and lb are made for the size
and topology, while the rest 8 ub and lb (16-23 are made for the shape
problem
%}

ub(1:barNumber)= size(groundStructure{2,caseNum},2);
lb(1:barNumber)= -size(groundStructure{2,caseNum},2);

%Number of nodes that possible to move
PnodeNumber = size(groundStructure{3,caseNum},1);

%Bb and lb for shape problem
ub(barNumber+1: barNumber+PnodeNumber)= groundStructure{3,caseNum}(:,1);
lb(barNumber+1: barNumber+PnodeNumber)= groundStructure{3,caseNum}(:,2);




