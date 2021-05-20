function err=write_adcirc_fort13(f13In,fname)
% WRITE_ADCIRC_FORT13 output an ADCIRC nodal attribute file (fort.13)
% err=write_adcirc_fort13(f13In);

err=0;

if ~exist('fname'),fname='fort.13';end

if exist(fname,'file')
    fprintf('output filename exists. Remove/move it first.\n');
    err=-1;
    return
end

fid=fopen(fname,'w');

% section 1
fprintf(fid,'%s\n',f13In.header);
fprintf(fid,'%d\n',f13In.nn);
fprintf(fid,'%d\n',f13In.natts);
for i=1:f13In.natts
    fprintf(fid,'%s\n',f13In.att_names{i});
    fprintf(fid,'%s\n',strtrim(f13In.att_units{i}));
    fprintf(fid,'%d\n',f13In.att_numdefvals(i));
    fmtstr=repmat('%f ',[1 f13In.att_numdefvals(i)]);
    fprintf(fid,[fmtstr '\n'],f13In.att_defvals{i});
end

%section 2
for i=1:f13In.natts
    fprintf('%s\n',f13In.att_names{i})
    fprintf(fid,'%s\n',f13In.att_names{i});
    if f13In.att_setnodes{i} == 0
        fprintf(fid,'0\n');
    else
        fprintf(fid,'%d\n',length(f13In.att_setnodes{i}));
        fmtstr=repmat(' %f',[1 f13In.att_numdefvals(i)]);
        iset=f13In.att_setnodes{i};
        out=f13In.atts{i}(iset,:);
        fprintf(fid,['%d' fmtstr '\n'],[iset out]');
    end
    
end

fclose(fid);
