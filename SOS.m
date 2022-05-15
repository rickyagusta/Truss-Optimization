clear;clc;close all; warning off;

%Adjust FE number, popsize and loop
caseNum=1;
maxFE=40000;
FEperIter=4;
popsize=100;
maxIter=maxFE/popsize/FEperIter;

%Initialize the controlling parameters and boundaries
[ub,lb]=boundary(caseNum);
nVar=length(ub);


fobj=@ObjectiveFunction;

x=zeros(popsize,nVar);     % Ecosystem Matrix
fv =zeros(popsize,1);      % Fitness Matrix


% --- Ecosystem Initialization
for i=1:popsize
    % Initialize the organisms randomly in the ecosystem
    x(i,:)=rand(1,nVar).*(ub-lb)+lb;
    fv(i,:)=fobj(x(i,:),caseNum,0,maxIter);
end

% Update the best Organism
[bestFv,idx]=min(fv); bestX=x(idx,:);

% --- Main Looping
for iter=1:maxIter
    
    for i=1:popsize % Organisms' Looping
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mutualism Phase
        % Choose organism j randomly other than organism i
        j=i;
        while i==j
            j=randi(popsize);
        end
        
        % Determine Mutual Vector & Beneficial Factor
        mutVec=mean([x(i,:);x(j,:)]);
        BF1=round(1+rand); BF2=round(1+rand);
        
        % Calculate new solution after Mutualism Phase
        xNew1=x(i,:)+rand(1,nVar).*(bestX-BF1.*mutVec);
        xNew2=x(j,:)+rand(1,nVar).*(bestX-BF2.*mutVec);
        
        % replace the out-of-boundary location
        index1=find(xNew1(1,:)>ub); index2=find(xNew1(1,:)<lb);
        xNew1(1,index1)=ub(index1); xNew1(1,index2)=lb(index2);
        
        index1=find(xNew2(1,:)>ub); index2=find(xNew2(1,:)<lb);
        xNew2(1,index1)=ub(index1); xNew2(1,index2)=lb(index2);
        
        % Evaluate the fitness of the new solution
        fvNew1=fobj(xNew1,caseNum,iter,maxIter);
        fvNew2=fobj(xNew2, caseNum,iter,maxIter);
        
        % Accept the new solution if the fitness is better
        if fvNew1<fv(i)
            fv(i)=fvNew1;
            x(i,:)=xNew1;
        end
        if fvNew2<fv(j)
            fv(j)=fvNew2;
            x(j,:)=xNew2;
        end
        
        % End of Mutualism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Commensialism Phase
        
        % Choose organism j randomly other than organism i
        j=i;
        while i==j
            j=randi(popsize);
        end
        
        % Calculate new solution after Commensalism Phase
        xNew1=x(i,:)+(rand(1,nVar)*2-1).*(bestX-x(j,:));
        
        % Replace the out-of-boundary location
        index1=find(xNew1(1,:)>ub); index2=find(xNew1(1,:)<lb);
        xNew1(1,index1)=ub(index1); xNew1(1,index2)=lb(index2);
        
        % Evaluate the fitness of the new solution
        fvNew1=fobj(xNew1, caseNum,iter,maxIter);
        
        % Accept the new solution if the fitness is better
        if fvNew1<fv(i)
            fv(i)=fvNew1;
            x(i,:)=xNew1;
        end
        
        % End of Commensalism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parasitism Phase
        
        % Choose organism j randomly other than organism i
        j=i;
        while i==j
            j=randi(popsize);
        end
        
        % Determine Parasite Vector & Calculate the fitness
        parVec=x(i,:);
        seed=randperm(nVar);
        id=seed(1:ceil(rand*nVar));  % select random dimension
        parVec(:,id)=rand(1,length(id)).*(ub(id)-lb(id))+lb(id);
        fvPar=fobj(parVec, caseNum,iter,maxIter);
        
        % Kill organism j and replace it with the parasite
        % if the fitness is lower than the parasite
        if fvPar < fv(j)
            fv(j)=fvPar;
            x(j,:)=parVec;
        end
        % End of Parasitism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Update the best Organism
        
        [bestFv,idx]=min(fv); bestX=x(idx,:);
    end % End of Organisms' Looping
    
    %Saving the best values
    bestRecord(iter,:)=bestFv;
    bestRecordX(iter,:)=bestX;
    outmsg=['Iter #',num2str(iter),' gbest=',num2str(bestFv)];
    disp(outmsg)
    
    %Plot truss figure
    [~,~,~, A, node_con, node_loc] = fobj(bestX, caseNum, iter, maxIter);
    trussPlot=drawing(bestFv, A, node_con, node_loc, iter);
    
end % End of Main Looping
[y,cons_total]=ObjectiveFunction(bestX, caseNum, iter, maxIter);
% --- Display the result
disp(' ');
disp(['Best Fitness: ', num2str(bestFv)])
disp(['Best Organism: ', num2str(bestX)])
disp(' ');