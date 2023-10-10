function new_index = jk_mod(index,total_n)
%
new_index = mod(index,total_n);

if new_index == 0
    new_index = total_n;
end

end

