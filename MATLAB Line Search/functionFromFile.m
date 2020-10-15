% Reads variables, objective function, and intialization point from file
%
% Format expected:
%      variables as column vector \n
%      objective function with respect to variables \n
%      initial point as column vector \n
% Example:
%      [u; v]
%      (2 - u*v)^2 + (v-u^2)^2
%      [1;1]

function [vars,obj,x0] = functionFromFile(filename)

fid = fopen(filename);

% Read all lines & collect in cell array
txt = textscan(fid,'%s','delimiter','\n');

s = size(txt{1});
if(s(1,1) ~= 3)
    error('Invalid input. Expected variables, objective, and initial x delimited by \n')
end
    
vars = str2sym(txt{1}{1});
obj = str2sym(txt{1}{2});
x0  = str2num(txt{1}{3});

if (size(x0) ~= size(vars))
    error('Invalid initial point. Must specify scalar value for each symbolic variable')  
end

end