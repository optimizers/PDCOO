pdcoo_home = pwd;

addpath(pdcoo_home)
addpath(fullfile(pdcoo_home, 'Formulations'))
addpath(fullfile(pdcoo_home, 'Solvers'))
addpath(fullfile(pdcoo_home, 'readmps'))
addpath(fullfile(pdcoo_home, 'Tools'))
if exist('Variants', 'dir') ~= 7
  mkdir('Variants');
end
addpath(fullfile(pdcoo_home, 'Variants'))
addpath(fullfile(pdcoo_home, 'Test'))
