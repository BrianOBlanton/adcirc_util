function f13Out=MapFort13(FgsIn,f13In,FgsOut,j)
% f13Out=MapFort13(FgsIn,f13In,FgsOut,j)

if ~is_valid_struct(FgsIn) || ~is_valid_struct(FgsOut)
    error('One or both of the input FGS structs are not valid')
end

if ~isfield(FgsIn,'strtree')
    error('Input FGS must have strtree field attached.  Call FgsIn.strtree=ComputeStrTree(FgsIn)')
end

% locate elements in FgsIn that contain nodes in FgsOut 
if ~exist('j','var') || isempty(j) 
    tic;
    fprintf('Finding elements in FgsIn containing nodes in FgsOut ... ');
    j=FindElementsInStrTree(FgsIn,FgsOut.x,FgsOut.y);
    v=toc;
    fprintf('Done. Took %f secs\n',v);
end

f13Out.header=FgsOut.name;
f13Out.nn=length(FgsOut.x);
f13Out.natts=f13In.natts;
f13Out.att_names=f13In.att_names;
f13Out.att_units=f13In.att_units;
f13Out.att_numdefvals=f13In.att_numdefvals;
f13Out.att_defvals=f13In.att_defvals;

tic
fprintf('Interpolating attributes ... ');
for i=1:length(f13In.att_names)
   
    f13Out.att_names{i}=strtrim(f13Out.att_names{i});
    
    f13Out.atts{i}=interp_scalar(FgsIn,f13In.atts{i},FgsOut.x,FgsOut.y,j);
    
    % attribute-specific processing
    if strcmp(f13Out.att_names{i},'sea_surface_height_above_geoid')
        f13Out.atts{i}=round(f13Out.atts{i}*10000)/10000; 
        f13Out.att_setnodes{i}=0; 
    end
         
    % mannings_n_at_sea_floor
    if strcmp(f13Out.att_names{i},'mannings_n_at_sea_floor')
        f13Out.atts{i}=round(f13Out.atts{i}*100)/100; 
        ikeep=abs(f13Out.atts{i}-f13Out.att_defvals{i})>0.001;
        f13Out.att_setnodes{i}=find(ikeep);

    end
    
    % primitive_weighting_in_continuity_equation
    if strcmp(f13Out.att_names{i},'primitive_weighting_in_continuity_equation')
        f13Out.atts{i}=round(f13Out.atts{i}*1000)/1000;
        ikeep=abs(f13Out.atts{i}-f13Out.att_defvals{i})>0.0001;
        f13Out.att_setnodes{i}=find(ikeep);
    end
    
    % surface_canopy_coefficient
    if strcmp(f13Out.att_names{i},'surface_canopy_coefficient')
        f13Out.atts{i}=round(f13Out.atts{i}); 
        ikeep=~(f13Out.atts{i}==f13Out.att_defvals{i});
        f13Out.att_setnodes{i}=find(ikeep);
    end
    
    % surface_submergence_state
    if strcmp(f13Out.att_names{i},'surface_submergence_state')
        f13Out.atts{i}=round(f13Out.atts{i}); 
    end
    
    % initial_river_elevation
    if strcmp(f13Out.att_names{i},'initial_river_elevation')
        idx=f13Out.atts{i}<-1000;
        f13Out.atts{i}(idx)=-99999; 
    end    
    
    % surface_directional_effective_roughness_length
    if strcmp(f13Out.att_names{i},'surface_directional_effective_roughness_length')
        i0=find(~all(f13Out.atts{i}'<1e-10));
        f13Out.att_setnodes{i}=i0(:);
%     else
%         f13Out.att_setnodes{i}=find(abs(f13Out.atts{i}-f13Out.att_defvals{i})>1e-5);
    end
    
    % elemental_slope_limiter
    if strcmp(f13Out.att_names{i},'elemental_slope_limiter')
        idx=f13Out.atts{i}>1000;
        f13Out.atts{i}(idx)=99999; 
        idx=find(f13Out.atts{i}<=1000);
        for n=1:length(idx)
            jn=idx(n);
            qn=FgsIn.e(j(jn),:);
            v=f13In.atts{i}(qn);
            dx=FgsIn.x(qn)-FgsOut.x(idx(n));
            dy=FgsIn.y(qn)-FgsOut.y(idx(n));
            [a,b]=min(abs(dx+1j*dy));
            f13Out.atts{i}(idx(n))=f13In.atts{i}(qn(b));
        end
        ikeep=abs(f13Out.atts{i}-f13Out.att_defvals{i})>1000;
        f13Out.att_setnodes{i}=find(ikeep);
    end    

    
end

v=toc;
fprintf('Done. Took %f secs\n',v);
