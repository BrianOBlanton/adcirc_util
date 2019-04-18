function f13Out=MapFort13(FgsIn,f13In,FgsOut,j)

if ~is_valid_struct(FgsIn) || ~is_valid_struct(FgsOut)
    error('One or both of the input FGS structs are not valid')
end

if ~isfield(FgsIn,'strtree')
    error('Input FGS must have strtree field attached.  Call FgsIn.strtree=ComputeStrTree(FgsIn)')
end

% locate elements in FgsIn that contain nodes in FgsOut 
if ~exist('j','var') || isempty(j) 
    fprintf('Finding elements in FgsIn containing nodes in FgsOut ... ');
    j=FindElementsInStrTree(FgsIn,FgsOut.x,FgsOut.y);
    fprintf('done.\n');
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
    
    % mannings_n_at_sea_floor
    
    if strcmp(f13Out.att_names{i},'mannings_n_at_sea_floor')
        f13Out.atts{i}=round(f13Out.atts{i}*100)/100; 
    end
    
    % primitive_weighting_in_continuity_equation
    if strcmp(f13Out.att_names{i},'primitive_weighting_in_continuity_equation')
        f13Out.atts{i}=round(f13Out.atts{i}*1000)/1000; 
    end
    
    % surface_canopy_coefficient
    if strcmp(f13Out.att_names{i},'surface_canopy_coefficient')
        f13Out.atts{i}=round(f13Out.atts{i}); 
    end
    
    % surface_submergence_state
    if strcmp(f13Out.att_names{i},'surface_submergence_state')
        f13Out.atts{i}=round(f13Out.atts{i}); 
    end
    
    % surface_directional_effective_roughness_length
    if strcmp(f13Out.att_names{i},'surface_directional_effective_roughness_length')
        i0=find(~all(f13Out.atts{i}'<1e-10));
        f13Out.att_setnodes{i}=i0(:);
        
    else
        f13Out.att_setnodes{i}=find(abs(f13Out.atts{i}-f13Out.att_defvals{i})>1e-5);
    end
end

v=toc;
fprintf('Done. Took %f secs\n',v);
