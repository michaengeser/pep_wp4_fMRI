function resp = getDAQKey(daq_var, bits)

resp = read(daq_var,1, "OutputFormat", "Matrix");
while ~any(resp(bits))
    resp = read(daq_var,1, "OutputFormat", "Matrix");
end
resp = find(resp(bits));

return

persistent keybits;
persistent key_1_bits;
persistent key_2_bits;
persistent key_4_bits;
persistent key_3_bits;

if isempty(keybits)
    keybits = sum(2.^(0:6)) + sum(2.^(8:10));
end

if isempty(key_1_bits)
    key_1_bits = sum(2.^(0));
end
if isempty(key_2_bits)
    key_2_bits = sum(2.^(1:2));
end
if isempty(key_4_bits)
    key_4_bits = sum(2.^(3:6));
end
if isempty(key_3_bits)
    key_3_bits = sum(2.^(8:10));
end


disp(keybits);


%resp = bitand(res, keybits);
