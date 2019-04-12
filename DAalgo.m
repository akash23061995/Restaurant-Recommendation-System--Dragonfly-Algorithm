function [Best_pos]=DAalgo(datas20,tr20,tnr20,nar)
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers
format long;

SearchAgents_no=size(datas20,2);
Max_iteration=200;
lb=zeros(size(datas20,1),1);
lb=lb';
ub=zeros(size(datas20,1),1);
ub=ub';
for i=1:length(ub)
    lb(i)=0;
    ub(i)=2;
end
datas20=double(datas20);
dim=size(datas20,1);

display('DA is optimizing your problem');
if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end
w1=0.5;
w2=0.5;
%The initial radius of dragonflies neighbourhoods
r=(ub-lb)/10;
Delta_max=(ub-lb)/10;
Food_fitness=inf;
Food_pos=zeros(dim,1);
Enemy_fitness=-inf;
Enemy_pos=zeros(dim,1);
X=initialization(SearchAgents_no,dim,ub,lb);
Fitness=zeros(1,SearchAgents_no);
DeltaX=initialization(SearchAgents_no,dim,ub,lb);
for iter=1:Max_iteration
    
    r=(ub-lb)/4+((ub-lb)*(iter/Max_iteration)*2);
    
    w=0.9-iter*((0.9-0.4)/Max_iteration);
    
    for i=1:SearchAgents_no %Calculate all the objective values first
        %          Fitness(1,i)=CostFunction(X(:,i)');
        Fitness(1,i)=CostFunction(datas20(:,i)');
        All_fitness(1,i)=Fitness(1,i);
        if Fitness(1,i)<Food_fitness
            Food_fitness=Fitness(1,i);
            Food_pos=X(:,i);
        end
        
        if Fitness(1,i)>Enemy_fitness
            if all(X(:,i)<ub') && all( X(:,i)>lb')
                Enemy_fitness=Fitness(1,i);
                Enemy_pos=X(:,i);
            end
        end
    end
    
    for i=1:SearchAgents_no
        index=0;
        neighbours_no=0;
        
        clear Neighbours_V
        clear Neighbours_DeltaX
        %find the neighbouring solutions
        for j=1:SearchAgents_no
            Dist2Enemy=distance(X(:,i),X(:,j));
            if (all(Dist2Enemy<=r) && all(Dist2Enemy~=0))
                index=index+1;
                neighbours_no=neighbours_no+1;
                Neighbours_DeltaX(:,index)=DeltaX(:,j);
                Neighbours_X(:,index)=X(:,j);
            end
        end
        
        % Seperation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S=zeros(dim,1);
        if neighbours_no>1
            for k=1:neighbours_no
                S=S+(Neighbours_X(:,k)-X(:,i));
            end
            S=-S;
        else
            S=zeros(dim,1);
        end
        
        % Alignment%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if neighbours_no>1
            A=(sum(Neighbours_DeltaX')')/neighbours_no;
        else
            A=DeltaX(:,i);
        end
        
        % Cohesion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if neighbours_no>1
            C_temp=(sum(Neighbours_X')')/neighbours_no;
        else
            C_temp=X(:,i);
        end
        
        C=C_temp-X(:,i);
        
        % Attraction to food%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dist2Food=distance(X(:,i),Food_pos(:,1));
        if iter~=1
            if all(Dist2Food<=r)
                F=Food_pos-X(:,i);
            else
                F=0;
            end
        else
            if all(tr20<nar)
                F=tr20;
            else
                F=0;
            end
        end
        % Distraction from enemy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dist2Enemy=distance(X(:,i),Enemy_pos(:,1));
        if iter~=1
            if all(Dist2Enemy<=r)
                Enemy=Enemy_pos+X(:,i);
            else
                Enemy=zeros(dim,1);
            end
        else
            if all(tnr20>=1)
                Enemy=tnr20;
            else
                Enemy=zeros(dim,1);
            end
        end
        for tt=1:dim
            if X(tt,i)>ub(tt)
                X(tt,i)=lb(tt);
                DeltaX(tt,i)=rand;
            end
            if X(tt,i)<lb(tt)
                X(tt,i)=ub(tt);
                DeltaX(tt,i)=rand;
            end
        end
        
        if any(Dist2Food>r)
            if neighbours_no>1
                for j=1:dim
                    DeltaX(j,i)=w*DeltaX(j,i)+rand*A(j,1)+rand*C(j,1)+rand*S(j,1);
                    if DeltaX(j,i)>Delta_max(j)
                        DeltaX(j,i)=Delta_max(j);
                    end
                    if DeltaX(j,i)<-Delta_max(j)
                        DeltaX(j,i)=-Delta_max(j);
                    end
                    X(j,i)=X(j,i)+DeltaX(j,i);
                end
            else
                X(:,i)=X(:,i)+Levy(dim)'.*X(:,i);
                DeltaX(:,i)=0;
            end
        else
            for j=1:dim
                %                 DeltaX(j,i)=(a*A(j,1)+c*C(j,1)+s*S(j,1)+f*F(j,1)+e*Enemy(j,1)) + w*DeltaX(j,i);
                DeltaX(j,i)=w1*(A(j,1)+C(j,1)+S(j,1)+F(j,1)+Enemy(j,1)) + w2*DeltaX(j,i);
                if DeltaX(j,i)>Delta_max(j)
                    DeltaX(j,i)=Delta_max(j);
                end
                if DeltaX(j,i)<-Delta_max(j)
                    DeltaX(j,i)=-Delta_max(j);
                end
                X(j,i)=X(j,i)+DeltaX(j,i);
            end
        end
        
        Flag4ub=X(:,i)>ub';
        Flag4lb=X(:,i)<lb';
        X(:,i)=(X(:,i).*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;
        
    end
    Best_score=Food_fitness;
    Best_pos=Food_pos;
    
end
