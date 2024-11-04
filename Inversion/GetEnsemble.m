function cX = GetEnsemble(Step,n,Ne,X)
cX = zeros(n,Ne);
for kk=1:Ne
    cX(:,kk) = X{kk}(:,Step);
end