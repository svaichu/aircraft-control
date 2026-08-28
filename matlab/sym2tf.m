function [G, v] = sym2tf(Gsym, varargin)
%SYM2TF Transfer-function model from a symbolic transfer-function matrix.
%   G = sym2tf(Gs) converts the symbolic scalar or matrix Gs into a tf model
%   of the same size, entry by entry. It is the inverse of TF2SYM: symbolic
%   work goes back to a model you can bode/step/margin.
%
%   The variable is picked up from Gs itself. A model in z comes back
%   discrete with an unspecified sample time (Ts = -1); anything else comes
%   back continuous.
%
%   G = sym2tf(Gs, Ts) sets the sample time explicitly: a positive Ts for a
%   discrete model, 0 for a continuous one.
%
%   G = sym2tf(Gs, ..., 'x') reads Gs as a function of x, which settles the
%   ambiguity when the expression carries more than one symbol.
%
%   [G, v] = sym2tf(...) also returns the variable that was used.
%
%   Each entry is split with numden, so a common denominator is fine and
%   z^-k delay factors simply land in the denominator. Every entry must be
%   rational in one variable with numeric coefficients - a leftover
%   parameter or an exp(-s*T) delay is an error, since neither has a
%   coefficient-vector form. Use SUBS for the first and PADE for the second.
%
%   See also TF2SYM, ZERORECTANGLETF.

    Gsym = sym(Gsym);
    [Ts, v] = parseArgs(Gsym, varargin{:});

    [ny, nu] = size(Gsym);
    num = cell(ny, nu);
    den = cell(ny, nu);
    for i = 1:ny
        for j = 1:nu
            [num{i, j}, den{i, j}] = entryCoeffs(Gsym(i, j), v, i, j);
        end
    end

    G = tf(num, den, Ts);
end

function [numCoef, denCoef] = entryCoeffs(e, v, i, j)
%ENTRYCOEFFS Numerator and denominator coefficient vectors of one entry.
    [n, d] = numden(simplifyFraction(e));

    numCoef = polyCoeffs(n, v, i, j);
    denCoef = polyCoeffs(d, v, i, j);
    lead = find(denCoef ~= 0, 1);
    if isempty(lead)
        error('sym2tf:ZeroDenominator', 'Entry (%d,%d) has a zero denominator.', i, j);
    end

    % numden clears fractions, which leaves entries scaled by whatever the
    % common denominator was; rescale so the denominator comes back monic and
    % a tf2sym/sym2tf roundtrip returns the coefficients it started with.
    numCoef = numCoef / denCoef(lead);
    denCoef = denCoef / denCoef(lead);
end

function c = polyCoeffs(p, v, i, j)
%POLYCOEFFS Descending-power coefficients of p, checked for convertibility.
    leftover = setdiff(symvar(p), v);
    if ~isempty(leftover)
        error('sym2tf:NonNumericCoefficients', ...
              ['Entry (%d,%d) still depends on %s. A tf holds numeric ' ...
               'coefficients only - substitute a value first.'], ...
              i, j, strjoin(cellstr(string(leftover)), ', '));
    end

    try
        c = sym2poly(expand(p));
    catch
        error('sym2tf:NotRational', ...
              ['Entry (%d,%d) is not rational in %s (a delay term such as ' ...
               'exp(-%s*T) is not). Approximate it with pade first.'], ...
              i, j, char(v), char(v));
    end
end

function [Ts, v] = parseArgs(Gsym, varargin)
%PARSEARGS Trailing char argument is the variable, a numeric one is Ts.
    varName = '';
    if ~isempty(varargin) && (ischar(varargin{end}) || isstring(varargin{end}))
        varName = char(varargin{end});
        varargin(end) = [];
    end

    if isempty(varName)
        v = mainVariable(Gsym);
    else
        v = sym(varName);
    end

    if ~isempty(varargin)
        Ts = varargin{1};
        if ~isnumeric(Ts) || ~isscalar(Ts)
            error('sym2tf:BadSampleTime', 'Ts must be a numeric scalar.');
        end
    elseif strcmp(char(v), 'z')
        Ts = -1;   % discrete, sample time unspecified
    else
        Ts = 0;
    end
end

function v = mainVariable(Gsym)
%MAINVARIABLE The variable Gsym is expressed in (defaults to z).
    v = symvar(Gsym);
    if isempty(v)
        v = sym('z');
    elseif numel(v) > 1
        error('sym2tf:AmbiguousVariable', ...
              'Expression holds several symbols (%s); name the variable explicitly.', ...
              strjoin(cellstr(string(v)), ', '));
    else
        v = v(1);
    end
end
