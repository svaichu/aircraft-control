function [Gsym, v] = tf2sym(G, varName)
%TF2SYM Symbolic transfer-function matrix from an LTI model.
%   Gs = tf2sym(G) converts the SISO or MIMO model G (tf, zpk or ss) into a
%   symbolic matrix in the variable z, sized ny-by-nu like G itself. A SISO
%   model gives back a 1-by-1 sym, so the result feeds straight into
%   zeroRectangleTf and the rest of the symbolic tools.
%
%   Gs = tf2sym(G, 'x') uses the variable x instead ('s' for a continuous
%   model, say).
%
%   [Gs, v] = tf2sym(...) also returns the symbolic variable used.
%
%   Coefficients are rationalised on the way in (sym's default for doubles),
%   so 0.1 becomes 1/10 and gcd/lcm/simplify stay exact rather than drifting
%   on floating-point noise. Each entry keeps its own num/den exactly as the
%   model states it, scaled to a monic denominator - nothing is cancelled, so
%   a pole/zero pair that overlaps inside an entry survives the conversion.
%
%   Transport delays are carried over: z^-d for a discrete model, exp(-v*d)
%   for a continuous one, with d the total I/O delay of that entry.
%
%   See also SYM2TF, ZERORECTANGLETF, POLY2SYM, TFDATA.

    if nargin < 2 || isempty(varName)
        varName = 'z';
    end
    v = sym(varName);

    if ~isa(G, 'tf')
        G = tf(G);          % zpk, ss, frd-free models all convert cleanly
    end

    [num, den] = tfdata(G);  % cell arrays, ny-by-nu, even in the SISO case
    delays = totalDelays(G);

    [ny, nu] = size(num);
    Gsym = sym(zeros(ny, nu));
    for i = 1:ny
        for j = 1:nu
            Gsym(i, j) = entry2sym(num{i, j}, den{i, j}, v) ...
                       * delayFactor(delays(i, j), v, G.Ts);
        end
    end
end

function e = entry2sym(numCoef, denCoef, v)
%ENTRY2SYM One num/den coefficient pair as a symbolic rational function.
%   The pair is scaled to a monic denominator and otherwise left alone: no
%   simplify, so a pole/zero pair that cancels inside an entry stays visible.
%   Smith-McMillan needs it that way - cancelling here would drop poles from
%   the pole polynomial.
    numCoef = numCoef(:).';
    denCoef = denCoef(:).';

    lead = find(denCoef ~= 0, 1);
    if isempty(lead)
        error('tf2sym:ZeroDenominator', 'Transfer function entry has a zero denominator.');
    end
    numCoef = sym(numCoef) / sym(denCoef(lead));
    denCoef = sym(denCoef) / sym(denCoef(lead));

    e = poly2sym(numCoef, v) / poly2sym(denCoef, v);
end

function f = delayFactor(d, v, Ts)
%DELAYFACTOR Delay d rendered in the model's own variable.
%   A discrete model counts d in samples, so the factor is z^-d; a
%   continuous one counts it in seconds, giving exp(-v*d).
    if d == 0
        f = sym(1);
    elseif Ts ~= 0
        f = v^(-sym(round(d)));
    else
        f = exp(-v * sym(d));
    end
end

function d = totalDelays(G)
%TOTALDELAYS Total I/O delay per entry, zeros if the model carries none.
    try
        d = totaldelay(G);
    catch
        d = zeros(size(G));
    end
end
