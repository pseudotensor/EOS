% sets-up grid related to eos_extract.m
function [min max step funout] = setupgrid(numindex)
  
% base grid upon integers first, so really have correct number of grid points
  min=0;
  max=1.0;
  step=(max-min)/(1.0d0*(numindex-1));
  
  grid=0:1:numindex-1;
  % now varies from 0..1
  grid=min + grid.*step;
  % offset so no issues with interpolations
  grid(1)=grid(1)+1E-15;
  grid(numindex)=grid(numindex)-1E-15;
  funout = grid;
  
end



