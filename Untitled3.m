subs = [99999,99999,99999]; %<-- Subscripts of the nonzeros.
vals = [0]; %<-- The values of the nonzeros.
X = sptensor(subs,vals) %<-- Create a sparse tensor with 3 nonzeros.