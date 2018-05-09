function s = chkin(s, field)
% Purpose:
% Adds or creates one or more empty fields specified by cell array "field"


for k=1:length(field)
   s.(field{k}) = [];
end

end  