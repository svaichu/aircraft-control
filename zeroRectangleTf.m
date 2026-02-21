function zerosOut = zeroRectangleTf(G)
%ZERORECTANGLETf Compute invariant zeros from a transfer matrix.
%   zerosOut = zeroRectangleTf(G) returns the zeros derived from the
%   transfer-function matrix G (symbolic). The routine mirrors the Python
%   zeroRectangleTf() implementation.

    G = sym(G);
    r = rank(G);
    minors = getMinors(G, r);
    polePoly = lcdDenom(minors);
    zeroPoly = gcdNumAdjusted(minors, polePoly);

    z = sym('z');
    zerosOut = solve(zeroPoly, z);
    zerosOut = reshape(zerosOut, [], 1);
end

function minors = getMinors(M, order)
%GETMINORS Compute all rank-order minors (determinants) of M.
    [nRows, nCols] = size(M);
    if order > max(nRows, nCols)
        error('getMinors:OrderTooLarge', ...
              'Order %d exceeds matrix dimensions (%d x %d).', order, nRows, nCols);
    end

    if order == nRows
        comb = nchoosek(1:nCols, order);
        minors = sym.empty(0, 1);
        for k = 1:size(comb, 1)
            cols = comb(k, :);
            minorMatrix = M(:, cols);
            minors(end+1, 1) = det(minorMatrix); %#ok<AGROW>
        end
        return;
    elseif order == nCols
        comb = nchoosek(1:nRows, order);
        minors = sym.empty(0, 1);
        for k = 1:size(comb, 1)
            rows = comb(k, :);
            minorMatrix = M(rows, :);
            minors(end+1, 1) = det(minorMatrix); %#ok<AGROW>
        end
        return;
    else
        error('getMinors:UnsupportedCase', ...
              'Order must match the row count or column count.');
    end
end

function lcd = lcdDenom(minors)
%LCDDENOM Least common denominator across symbolic expressions.
    if isempty(minors)
        lcd = sym(1);
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
%GCDNUMADJUSTED Normalize minors to a common denominator then take GCD.
    if isempty(minors)
        gcdPoly = sym(1);
        return;
    end

    adjusted = sym.empty(0, 1);
    for k = 1:numel(minors)
        [numPart, denPart] = numden(minors(k));
        if ~isequal(denPart, polePoly)
            quo = simplify(denPart / polePoly);
        else
            quo = sym(1);
        end
        [quoNum, quoDen] = numden(quo);
        if quoNum ~= sym(1)
            error('gcdNumAdjusted:UnsupportedAdjustment', ...
                  'Cannot adjust minor %d to the common denominator.', k);
        end
        adjusted(end+1, 1) = numPart * quoDen; %#ok<AGROW>
    end

    gcdPoly = adjusted(1);
    for k = 2:numel(adjusted)
        gcdPoly = gcd(gcdPoly, adjusted(k));
    end
    gcdPoly = simplify(gcdPoly);
end
