function generateLookupCPPandHeaderFiles(Type,x_query,fx_interp)
% Arrange it to comma-separated string
str1 = num2str(x_query,'%0.4f\n');
str1 = regexprep(str1,'\s+','f,');

str2 = num2str(fx_interp,'%0.4f\n');
str2 = regexprep(str2,'\s+','f,');

LookupTableDir = 'LookupTables';
if ~exist(LookupTableDir, 'dir')
    mkdir(LookupTableDir)
end

% Save as 'TypelookupTable.cpp' header file
fid = fopen([LookupTableDir,filesep,Type,'lookupTable.cpp'],'w');
fprintf(fid,'// JLK store lookup table as header file for absolute B1 map calculation. \n');
fprintf(fid,['#include "',Type,'lookupTable.h"\n']);
fprintf(fid,'#include <vector>\n');
fprintf(fid,['std::vector <float> fx_lookupTable',Type,' = {%sf};\n'],str1);
fprintf(fid,'\n');
fprintf(fid,['std::vector <float> ffx_lookupTable',Type,' = {%sf};\n'],str2);
fprintf(fid,'\n');
fprintf(fid,['extern float flookupTableMaximum',Type,' = ',num2str(max(fx_interp),'%0.4f\n'),'f;\n']);
fprintf(fid,['extern float flookupTableMinimum',Type,' = ',num2str(min(fx_interp),'%0.4f\n'),'f;\n']);
fprintf(fid,['extern int ilookupTableLength',Type,' = ',num2str(length(fx_interp)),';\n']);
fclose(fid);

% Save as 'TypelookupTable.h' cpp file
fid = fopen([LookupTableDir,filesep,Type,'lookupTable.h'],'w');
fprintf(fid,'// JLK store lookup table as header file for absolute B1 map calculation. \n');
fprintf(fid,'#pragma once\n');
fprintf(fid,'#include <vector> \n');
fprintf(fid,['extern std::vector<float> fx_lookupTable',Type,';\n']);
fprintf(fid,['extern std::vector<float> ffx_lookupTable',Type,';\n']);
fprintf(fid,['extern float flookupTableMaximum',Type,';\n']);
fprintf(fid,['extern float flookupTableMinimum',Type,';\n']);
fprintf(fid,['extern int ilookupTableLength',Type,';\n']);
fclose(fid);
end