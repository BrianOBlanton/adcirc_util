function out=ArrayListToCell(in)

out=cell(in.size,1);

for i=0:in.size-1
   out{i+1}=in.get(i);
end
