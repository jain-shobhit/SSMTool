%% extract cell array of names from sym array
function xn=sco_names_from_sym(xs)
xn=cell(size(xs));
for i=1:numel(xs)
    xn{i}=char(xs(i));
end
end
