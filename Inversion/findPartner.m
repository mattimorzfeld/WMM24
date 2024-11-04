function pind = findPartner(ind,Ne)
go = 1;
while go == 1
    pind = randi(Ne); % index of random partner
    if pind ~= ind
        go = 0;
    end
end