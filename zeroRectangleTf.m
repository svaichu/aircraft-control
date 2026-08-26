function [zerosOut, zeroPoly, polePoly] = zeroRectangleTf(G)
%ZERORECTANGLETF Transmission zeros of a (possibly rectangular) transfer matrix.
%   z = zeroRectangleTf(G) returns the transmission zeros of the symbolic
%   transfer-function matrix G, repeated according to their multiplicity.
%
%   [z, zPoly, pPoly] = zeroRectangleTf(G) also returns the zero polynomial
%   and the pole polynomial. Both are defined only up to a non-zero constant
%   factor and are returned in whatever factored form gcd/lcm produce.
%
%   Method (Smith-McMillan, cf. Skogestad & Postlethwaite Sec. 4.5):
%     r       = normal rank of G
%     p(s)    = least common denominator of ALL non-zero minors of G of
%               every order 1..r                              -> pole polynomial
%     z(s)    = greatest common divisor of the numerators of the order-r
%               minors, after each has been rewritten over p(s)
%
%   Note the pole polynomial needs the minors of every order, not just the
%   order-r ones; using only the order-r minors loses the poles that cancel
%   inside the large minors and corrupts z(s) as well.

    G = sym(G);
    v = mainVariable(G);

    r = rank(G);
    if r == 0
        zeroPoly = sym(1);
        polePoly = sym(1);
        zerosOut = sym.empty(0, 1);
        return;
    end

    allMinors = sym.empty(0, 1);
    for k = 1:r
        allMinors = [allMinors; getMinors(G, k)]; %#ok<AGROW>
    end
    topMinors = getMinors(G, r);

    % Both polynomials are only defined up to a non-zero constant, so they
    % are left exactly as gcd/lcm produce them - factored, not expanded.
    polePoly = lcdDenom(allMinors);
    zeroPoly = gcdNumAdjusted(topMinors, polePoly);

    zerosOut = polyRoots(zeroPoly, v);
end

function v = mainVariable(G)
%MAINVARIABLE The variable G is expressed in (defaults to s).
    v = symvar(G);
    if isempty(v)
        v = sym('s');
    else
        v = v(1);
    end
end

function minors = getMinors(M, order)
    [nRows, nCols] = size(M);
    minors = sym.empty(0, 1);
    if order < 1 || order > min(nRows, nCols)
        return;
    end

    rowComb = nchoosek(1:nRows, order);
    colComb = nchoosek(1:nCols, order);
    for i = 1:size(rowComb, 1)
        for j = 1:size(colComb, 1)
            minors(end+1, 1) = simplify(det(M(rowComb(i, :), colComb(j, :)))); %#ok<AGROW>
        end
    end
    minors = minors(~isAlways(minors == 0, 'Unknown', 'false'));
end

function lcd = lcdDenom(minors)
%LCDDENOM Least common denominator across symbolic expressions.
    lcd = sym(1);
    if isempty(minors)
        return;
    end

    [~, den] = numden(minors);
    lcd = den(1);
    for k = 2:numel(den)
        lcd = lcm(lcd, den(k));
    end
    lcd = simplify(lcd);
end

function gcdPoly = gcdNumAdjusted(minors, polePoly)
%GCDNUMADJUSTED Rewrite each minor over polePoly, then take the GCD of the
%   resulting numerators.
    gcdPoly = sym(1);
    if isempty(minors)
        return;
    end

    adjusted = sym.empty(0, 1);
    for k = 1:numel(minors)
        [numPart, denPart] = numden(minors(k));
        % polePoly is a common denominator, so polePoly/denPart is a polynomial.
        adjusted(end+1, 1) = expand(simplify(numPart * polePoly / denPart)); %#ok<AGROW>
    end

    gcdPoly = adjusted(1);
    for k = 2:numel(adjusted)
        gcdPoly = gcd(gcdPoly, adjusted(k));
    end
    gcdPoly = simplify(gcdPoly);
end

function rootsOut = polyRoots(p, v)
%POLYROOTS Roots of p, repeated by multiplicity (exact where possible).
    rootsOut = sym.empty(0, 1);
    if isempty(symvar(p))
        return;  % constant: no zeros
    end

    % solve may already repeat roots; collapse first, then count multiplicity
    % from the derivatives so the count is right either way.
    distinct = uniqueRoots(reshape(solve(p == 0, v, 'MaxDegree', 4), [], 1));
    for k = 1:numel(distinct)
        m = rootMultiplicity(p, v, distinct(k));
        rootsOut = [rootsOut; repmat(distinct(k), m, 1)]; %#ok<AGROW>
    end
    rootsOut = simplify(rootsOut);

    if numel(rootsOut) ~= polynomialDegree(p, v)
        warning('zeroRectangleTf:IncompleteRoots', ...
                'Found %d zeros for a degree-%d zero polynomial.', ...
                numel(rootsOut), polynomialDegree(p, v));
    end
end

function out = uniqueRoots(r)
%UNIQUEROOTS Drop duplicates, comparing symbolically rather than textually.
    out = sym.empty(0, 1);
    for k = 1:numel(r)
        isNew = true;
        for j = 1:numel(out)
            if isAlways(r(k) == out(j), 'Unknown', 'false')
                isNew = false;
                break;
            end
        end
        if isNew
            out(end+1, 1) = r(k); %#ok<AGROW>
        end
    end
end

function m = rootMultiplicity(p, v, r)
%ROOTMULTIPLICITY Smallest k with the k-th derivative of p non-zero at r.
    m = 0;
    d = p;
    while m <= polynomialDegree(p, v)
        if ~isAlways(subs(d, v, r) == 0, 'Unknown', 'true')
            break;
        end
        d = diff(d, v);
        m = m + 1;
    end
    m = max(m, 1);
end

function n = polynomialDegree(p, v)
    n = double(feval(symengine, 'degree', p, v));
end
