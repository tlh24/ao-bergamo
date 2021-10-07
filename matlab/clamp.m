function [k] = clamp(x, lo, hi)
k = max(min(x, hi),lo); 