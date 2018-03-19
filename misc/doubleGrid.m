function fgsout=doubleGrid(fgs)
% fgsout=doubleGrid(fgs);

%clear fgsout

newE=NaN*ones(fgs.ne*5,3);
newX=NaN*ones(fgs.nn*5,1);
newY=newX;
newZ=newX;

UsedEdges=sparse(fgs.ne*5,fgs.ne*5);

NewNodeCode=zeros(size(newX));

NewNodeCount=fgs.nn;

AllBndNodes=unique(fgs.bnd(:));

je1=[1 2]';
je2=[2 3]';
je3=[3 1]';

mm=floor(fgs.ne/100);

for j=1:fgs.ne
    
    if rem(j,mm)==1,fprintf('Pct Complete =   %4.0f\n',floor(100*j/fgs.ne)),end
    
    ThisE=fgs.e(j,:);
    ThisX=fgs.x(ThisE);
    ThisY=fgs.y(ThisE);
    ThisZ=fgs.z(ThisE);

    % note switched end point order to test for previous ...
    PreviouslyUsed1=full(UsedEdges(ThisE(2),ThisE(1)));
    PreviouslyUsed2=full(UsedEdges(ThisE(3),ThisE(2)));
    PreviouslyUsed3=full(UsedEdges(ThisE(1),ThisE(3)));
    
    %find([PreviouslyUsed1 PreviouslyUsed2 PreviouslyUsed3])
    
    l=length(find([PreviouslyUsed1 PreviouslyUsed2 PreviouslyUsed3]==0));
    
    if l==3   % unused element     
       
        n1=NewNodeCount+1;
        n2=NewNodeCount+2;
        n3=NewNodeCount+3;
        
        e1=[ThisE(1) n1 n3];
        e2=[ThisE(2) n2 n1];
        e3=[ThisE(3) n3 n2];
        e4=[n1 n2 n3];
        
        UsedEdges(ThisE(1),ThisE(2))=n1;
        UsedEdges(ThisE(2),ThisE(3))=n2;
        UsedEdges(ThisE(3),ThisE(1))=n3;
        
        xe=ThisX([je1 je2 je3]);
        ye=ThisY([je1 je2 je3]);
        ze=ThisZ([je1 je2 je3]);

        newE((j-1)*4+1:j*4,:)=[e1;e2;e3;e4];
        newX([n1 n2 n3])=mean(xe)';
        newY([n1 n2 n3])=mean(ye)';
        newZ([n1 n2 n3])=mean(ze)';
        
        NewNodeCount=NewNodeCount+3;
        
    elseif l==2  % one edge has been used
       
        % new node numbers
       
        n1=NewNodeCount+1;
        n2=NewNodeCount+2;
        
        if PreviouslyUsed1>0
     
            e1=[ThisE(1) PreviouslyUsed1 n2];
            e2=[ThisE(2) n1 PreviouslyUsed1];
            e3=[ThisE(3) n2 n1];
            e4=[PreviouslyUsed1 n1 n2];
            
            UsedEdges(ThisE(2),ThisE(3))=n1;
            UsedEdges(ThisE(3),ThisE(1))=n2;
        
            xe=ThisX([je2 je3]);
            ye=ThisY([je2 je3]);
            ze=ThisZ([je2 je3]);

        elseif PreviouslyUsed2>0
            
            e1=[ThisE(1) n2 n1];
            e2=[ThisE(2) n1 PreviouslyUsed2 ];
            e3=[ThisE(3) PreviouslyUsed2 n2];
            e4=[PreviouslyUsed2 n1 n2];
            
            UsedEdges(ThisE(1),ThisE(2))=n1;
            UsedEdges(ThisE(3),ThisE(1))=n2;
            
            xe=ThisX([je1 je3]);
            ye=ThisY([je1 je3]);
            ze=ThisZ([je1 je3]);
            
        else
            
            e1=[ThisE(1) n1 PreviouslyUsed3];
            e2=[ThisE(2) n2 n1 ];
            e3=[ThisE(3) PreviouslyUsed3 n2];
            e4=[PreviouslyUsed3 n1 n2];
                    
            UsedEdges(ThisE(1),ThisE(2))=n1;
            UsedEdges(ThisE(2),ThisE(3))=n2;
        
            xe=ThisX([je1 je2]);
            ye=ThisY([je1 je2]);
            ze=ThisZ([je1 je2]);

        end
        
        newE((j-1)*4+1:j*4,:)=[e1;e2;e3;e4];
        newX([n1 n2])=mean(xe)';
        newY([n1 n2])=mean(ye)';
        newZ([n1 n2])=mean(ze)';
        
        NewNodeCount=NewNodeCount+2;
        
    elseif l==1  % two edges have been used
        
        n1=NewNodeCount+1;

        if PreviouslyUsed1==0
        
            e1=[ThisE(1) n1 PreviouslyUsed3];
            e2=[ThisE(2) PreviouslyUsed2 n1];
            e3=[ThisE(3) PreviouslyUsed3 PreviouslyUsed2];
            e4=[n1 PreviouslyUsed2 PreviouslyUsed3];  
            UsedEdges(ThisE(1),ThisE(2))=n1;
            xe=ThisX(je1);
            ye=ThisY(je1);
            ze=ThisZ(je1);
        
        elseif PreviouslyUsed2==0
            
            e1=[ThisE(1) PreviouslyUsed1 PreviouslyUsed3];
            e2=[ThisE(2) n1 PreviouslyUsed1 ];
            e3=[ThisE(3) PreviouslyUsed3 n1];
            e4=[n1 PreviouslyUsed3 PreviouslyUsed1];  
            UsedEdges(ThisE(2),ThisE(3))=n1;
            xe=ThisX(je2);
            ye=ThisY(je2);
            ze=ThisZ(je2);
            
        else
            
            e1=[ThisE(1) PreviouslyUsed1 n1];
            e2=[ThisE(2) PreviouslyUsed1 PreviouslyUsed2 ];
            e3=[ThisE(3) n1 PreviouslyUsed2];
            e4=[n1 PreviouslyUsed2 PreviouslyUsed1];  
            UsedEdges(ThisE(3),ThisE(1))=n1;
            xe=ThisX(je3);
            ye=ThisY(je3);
            ze=ThisZ(je3);
            
        end
             
        newE((j-1)*4+1:j*4,:)=[e1;e2;e3;e4];
        newX(n1)=mean(xe)';
        newY(n1)=mean(ye)';
        newZ(n1)=mean(ze)'; 
        
        NewNodeCount=NewNodeCount+1;

    else  % all edges used; new elements entirely from previously split edges
        
        e1=[ThisE(1) PreviouslyUsed1 PreviouslyUsed3];
        e2=[ThisE(2) PreviouslyUsed2 PreviouslyUsed1];
        e3=[ThisE(3) PreviouslyUsed3 PreviouslyUsed2];
        e4=[PreviouslyUsed3 PreviouslyUsed1 PreviouslyUsed2];

        newE((j-1)*4+1:j*4,:)=[e1;e2;e3;e4];                
 
    end
    
end

newX=[fgs.x;newX];
newY=[fgs.y;newY];
newZ=[fgs.z;newZ];
newX(isnan(newX))=[];
newY(isnan(newY))=[];
newZ(isnan(newZ))=[];
newE(isnan(newE))=[];
newE=reshape(newE,length(newE)/3,3);


% create/fill output struct
fgsout.name=[fgs.name '_doubled'];
fgsout.e=newE;
fgsout.x=newX;
fgsout.y=newY;
fgsout.z=newZ;
fgsout.nn=length(fgsout.x);
fgsout.ne=length(fgsout.e);
fgsout.bnd=detbndy(fgsout.e);
if isfield(fgs,'A')
    fgsout=belint(fgsout);
end
if isfield(fgs,'ar')
    fgsout=el_areas(fgsout);
end
if isfield(fgs,'xecen')
    fgsout=attach_elem_centroids(fgsout);
end

%%
if isfield(fgs,'nopen')
    disp('Fixing open boundary ...')
    fgsout.nopen=fgs.nopen;
    newob=NaN*ones(fgs.nopennodes{1}*3);
    j=0;
    for i=1:fgs.nopennodes{1}-1
       % bnd segment pairs assumed in order of fgs.ob
       p1=fgs.ob{1}(i);
       p2=fgs.ob{1}(i+1);
       %disp([i p1 p2])
       if i==1 
           newob(j+1)=p1;
           newob(j+2)=full(UsedEdges(p1,p2));
           newob(j+3)=p2;
           j=j+3;
       else
           newob(j+1)=full(UsedEdges(p1,p2));
           newob(j+2)=p2;
           j=j+2;
       end
    end
    newob(isnan(newob))=[];
    fgsout.nopennodes={length(newob)};
    fgsout.ob{1}=newob;
    fgsout.elevation=length(newob);
end

%%
if isfield(fgs,'nland')
    disp('Fixing land boundaries ...')

    fgsout.nland=fgs.nland;
    nlandnodestotal=0;
    for l=1:fgs.nland
        j=0;
        newln=NaN*ones(fgs.nlandnodes(l)*3);
        for i=1:fgs.nlandnodes(l)-1
            % bnd segment pairs assumed in order of fgs.ln{l}
            p1=fgs.ln{l}(i);
            p2=fgs.ln{l}(i+1);
            %disp([i p1 p2])
            if i==1
                newln(j+1)=p1;
                newln(j+2)=full(UsedEdges(p1,p2));
                newln(j+3)=p2;
                j=j+3;
            else
                newln(j+1)=full(UsedEdges(p1,p2));
                newln(j+2)=p2;
                j=j+2;
            end
        end
        newln(isnan(newln))=[];
        fgsout.nlandnodes(l)=length(newln);
        fgsout.ln{l}=newln';
        nlandnodestotal=nlandnodestotal+length(newln);
    end
    fgsout.nlandnodestotal=nlandnodestotal;
    fgsout.ibtype=fgs.ibtype;
end

%%



%  
% test grid;
g.x=[1 0 2 1 3 0 4 3]';
g.y=[1 2 2 3 1 1 2 3]';
g.e=[1 3 2; 3 4 2; 1 5 3; 1 2 6; 3 5 7; 8 3 7; 8 4 3];
g.ne=length(g.e);
g.nn=length(g.x);
g.bnd=detbndy(g.e);
g.name='x';
g.z=zeros(size(g.x));


