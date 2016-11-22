function [S] = safe_rmfield(S,fields)

fields = fields(isfield(S,fields));

S=rmfield(S,fields);