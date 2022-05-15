function [y, cons_total,int_stress, A, node_con, node_loc, C1, C2, C3, C4, C7, C8] = ObjectiveFunction(x, caseNum, iter, maxiter)

%{
This one is exception, need to be coded manually for each cases (see the
Tejani, 2017 ref. Because x5 must be the same with x6, also similar for
y3=y4
%}
x = [x(1:15) x(16) x(16) x(17) x(17) x(18:23)];

%Load groundStructure matrix to know the ground structure
load data.mat

%Steel quality and specification
E=groundStructure{9,caseNum}(1);
rho=groundStructure{9,caseNum}(2);
fy=groundStructure{9,caseNum}(3);
fu=groundStructure{9,caseNum}(4);
r=groundStructure{2,caseNum}(2,:);

%Ground structure, including node connections
node_loc = groundStructure{4,caseNum};
node_list = groundStructure{5,caseNum};

%To update shape based on metaheuristic algorithms
for updateShape = 1: size(groundStructure{11,caseNum},1)
    
    node_loc(groundStructure{11,caseNum}(updateShape,1),groundStructure{11,caseNum}(updateShape,2))= x(groundStructure{1,caseNum} + updateShape);
    
end

%Because discrete problem, then need rounding
%Rounding ONLY for size and topology problems (the first 15 cells)
xSize = round(x(1:groundStructure{1,caseNum}));
B = xSize > 0;
%Determining the topology
node_con=node_list(B,:);
%Changing from continuous variables to discrete
D=xSize(B);
A= groundStructure{2,caseNum}(1,D);

%Initialization
ksg=zeros(2*length(node_loc));
disp = zeros(2*length(node_loc),1);
int_stress=zeros(length(node_con),1);
int_force=zeros(length(node_con),1);
pu_tensile=zeros(length(node_con),1);
pu_comp=zeros(length(node_con),1);
cons_total=0;

%Input passive DOF
p_dof = groundStructure{6,caseNum};
node = unique(node_con);

%Constraint 1 (validity of structure)
for validity=1:length(groundStructure{8,caseNum})
    
    if sum(node==groundStructure{8,caseNum}(validity))==1
        C1=0;
    else
        C1=1;
        y=10^9;
        cons_total=y;
        return
    end
    
end

%Constraint 2 (grubler's criterion for 2D truss)
grubler = 2*size(node,1)-size(A,2)-length(p_dof);
if grubler <= 0
    C2 = 0;
else
    C2 = 1;
    y=10^8;
    cons_total=y;
    return
end

%Input forces
F=zeros(2*length(node_loc),1);
F(groundStructure{7,caseNum}(2),1)= groundStructure{7,caseNum}(1);

%Generate matrix consists of 2 vectors each nodes
id = 1;
for i = 1:size(node_loc,1)
    node_id(i,:) = [id, id+1];
    id = id+2;
end

%DIRECT STIFFNESS METHOD
for j = 1:size(node_con,1)
    %Numbering node connectivity
    con_id(j,:) = [ node_id(node_con(j,1),:)...
        node_id(node_con(j,2),:)];
    %length
    length_x(j,1) = (node_loc(node_con(j,2),1)-node_loc(node_con(j,1),1));
    length_y(j,1) = (node_loc(node_con(j,2),2)-node_loc(node_con(j,1),2));
    length_tot(j,1) = sqrt(length_x(j)^2+length_y(j)^2);
    
    %Bar stiffness
    k(j) = E*A(1,j)./ length_tot(j,1); %N/m
    
    %Length projection
    cx(j,1)=length_x(j,1)/length_tot(j,1);
    cy(j,1)=length_y(j,1)/length_tot(j,1);
    
    %Transforming to global
    T{j} = [ cx(j,1) cy(j,1) 0       0
        0       0       cx(j,1) cy(j,1)];
    
    kml{j} = k(j)* [1 -1; -1 1];
    kmg{j} = T{j}' *kml{j} *T{j};
    ksg(con_id(j,:),con_id(j,:)) = ksg(con_id(j,:),con_id(j,:)) + kmg{j};
    
    %Calculating the mass of members
    m(j) = rho*length_tot(j)*A(j);
end

%To locate a_dof in node
dof = unique(con_id)';
for t=1:length(p_dof)
    checker = dof ~= p_dof(t);
    dof= dof(checker);
end
a_dof = dof;

%Constraint 2 (positive def.)
stability = eigs(ksg(a_dof,a_dof),1,'SA');

if rank(ksg(a_dof,a_dof)) == size(ksg(a_dof,a_dof),1) && rank(ksg(a_dof,a_dof)) ~= 0 && stability > 0
    C2=0;
else
    C2=1;
    y=10^7;
    cons_total=y;
    return
end

%Displacement
disp(a_dof,1)=ksg(a_dof,a_dof)\F(a_dof,1);

%CONSTRAINT HANDLING

%Constraint 3 (stress)
C3=0;
% stressLimit = groundStructure{10,caseNum}(1);
% is = abs(int_stress);
% C3_check = is > stressLimit;
% C3=sum(abs((is(C3_check)-stressLimit)/stressLimit));

%Constraint 4 (disp) no displacement constraint in the 15-bar truss case
C4=0;
% dispLimit= groundStructure{10,caseNum}(2);
% da(1,:)=abs(disp);
% C4_check = da > dispLimit; %inch
% C4=sum(abs((da(C4_check)-dispLimit)/dispLimit));

%Internal force and stress
for q=1:size(node_con,1)
    
    int_stress(q,:)=E/length_tot(q,1) * [-1 1] * T{q} * disp(con_id(q,:),:);
    int_force(q,:)=k(q) *[-1 1] * T{q}  * disp(con_id(q,:),:);
    
    %Split internal force into tensile and compression
    if int_force(q,1) >= 0
        pu_tensile(q,1) = int_force(q,1);
        pu_comp(q,1) = 0;
    else
        pu_tensile(q,1) = 0;
        pu_comp(q,1) = abs(int_force(q,1));
    end
    
end

%New Constraint 7 (tensile strength) AISC 2010 D2
%D2(a) tensile yielding in the gross section
pn1 = fy * A; %Newton
phi_pn1 = 0.9*pn1; %Newton
%D2(b) tensile rupture in the net section
pn2 = fu * A; %Newton
phi_pn2 = 0.75*pn2; %Newton

phi_pntensile = min(phi_pn1,phi_pn2);
phi_pntensile=phi_pntensile';


C7_check = pu_tensile > phi_pntensile;
C7=sum(abs((pu_tensile(C7_check(:,1),1)-phi_pntensile(C7_check(:,1)))./phi_pntensile(C7_check(:,1))));


%New Constraint 8 (compression strength)AISC 2010 E3
k_comp = 1;

fe = pi.^2.*E ./(k_comp.*length_tot./(r(1,D))).^2; %N/m2
checkrumus = k_comp.*length_tot./(r(1,D)) <= 4.71 * (E / fy)^0.5;

for w=1:size(length_tot)
    if checkrumus(w) == 1
        fcr(w) = 0.658^(fy/fe(w))*fy;
    else
        fcr(w) = 0.877*fe(w);
    end
end

pn_comp = fcr.*A; %Newton
phi_pncomp = 0.9*pn_comp; %Newton
phi_pncomp =phi_pncomp';
C8_check = pu_comp > phi_pncomp;
C8=sum(abs((pu_comp(C8_check(:,1),1)-phi_pncomp(C8_check(:,1)))./phi_pncomp(C8_check(:,1))));

cons_total = C3 + C4 + C7 + C8 ; %only stress constraint

if cons_total == 0
    y= sum(m);
else
    %Penalty
    y= sum(m)*(1+(iter/maxiter)+(3*cons_total)^3);%tejani2018
end







