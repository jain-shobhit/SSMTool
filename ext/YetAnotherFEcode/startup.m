clc

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

addpath(genpath(strcat(pwd,sslash,'src')));
addpath(genpath(strcat(pwd,sslash,'external')));
addpath(genpath(strcat(pwd,sslash,'examples',sslash,'Meshes')));

disp('              _____ _____     ')
disp('  _   _  __ _|  ___| ____|___ ')
disp(' | | | |/ _` | |_  |  _| / __|')
disp(' | |_| | (_| |  _| | |__| (__ ')
disp('  \__, |\__,_|_|   |_____\___|')
disp('  |___/       YetAnotherFEcode')
fprintf('\n\n') 