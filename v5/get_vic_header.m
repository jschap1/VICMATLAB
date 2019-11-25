% Gets the header from a VIC 5 classic output file

function header_names = get_vic_header(fname, lineno)

eb_ID = fopen(fname);
for i=1:lineno % cycles through to the third line
    eb_header = fgetl(eb_ID);
end
header_names = strsplit(eb_header);

end