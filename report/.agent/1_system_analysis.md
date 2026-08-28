<!-- @.agent/AGENT.md is the entry point for the agent -->
\section{System Analysis}

Task: Find zeros of a given MIMO system for two cases:
\begin{itemize}
    \item Number of inputs is equal to the number of outputs.
    \item Number of inputs is not equal to the number of outputs.
\end{itemize}
Either case, MATLAB's control system toolbox or python's control library are not used.

\subsection{State space}

When the system is given in the state space form, Rosenbrock-matrix is used.

\begin{equation}
    P(s) = \begin{bmatrix}
        sI - A & -B \\ C & D
    \end{bmatrix}
    \label{eq:rosenbrock}
\end{equation}

Determinant of the above matrix is computed, which upon solving for $s$ gives the system zeros.

\begin{align*}
    \text{det}(P(s)) &= 0 \\
    (sI - A)(D) - (-B)(C) &= 0 \\
\end{align*}

Matlab implementation of the above procedure is given in the code snippet below. The function \texttt{getZeros} takes the state space matrices $A$, $B$, $C$, and $D$ as input and returns the system zeros.

\begin{lstlisting}
function zero = getZeros(A, B, C, D)
    syms s; % Define symbolic variable
    n = size(A,1);
    M = [s*eye(n)-A, -B; C, D];
    zero_eq = det(M); % Compute determinant symbolically
    zero = vpasolve(zero_eq, s); % Use vpasolve for numerical solutions
    zero = double(zero); % Convert symbolic results to double
end
\end{lstlisting}

The above function can be compared with MATLAB's built-in \texttt{tzero} function.
\begin{lstlisting}
p = 4; % input size
m = p; % output size
s = tf('s');
sys = rss(6, m, p);
A = sys.A; B = sys.B; C = sys.C; D = sys.D;

z1 = getZeros(A, B, C, D);
z2 = tzero(sys);
z1 = sort(z1(:));
z2 = sort(z2(:));
\end{lstlisting}

% assert(numel(z1) == numel(z2) && all(abs(z1 - z2) < 1e-6), ...
%     'getZeros and tzero do not match within tolerance');
% disp('getZeros and tzero match within tolerance');



\section{Transfer Function}

We maynot be always be able to find determinant of the Rosenbrock matrix, especially when the number of inputs is not equal to the number of outputs. In such cases, we can use the transfer function representation of the system to find the zeros. Smith-McMillan method is used to determine transmission zeros of the given MIMO transfer function. Since the task requires not to use MATLAB's control system toolbox, I have implemented using symbolic math toolbox and converted the transfer function to symbolic form by replacing the complex variable $s$ with a symbolic variable.


\begin{lstlisting}
function [zerosOut, zeroPoly, polePoly] = zeroRectangleTf(G)

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
\end{lstlisting}


And evaluating for the given example transfer function, we get the following results:

\begin{lstlisting}
s = tf('s');
G = [ (s-1)*(s+2), 0, (s-1)^2;
    -(s+1)*(s+2), (s-1)*(s+1), (s-1)*(s+1) ]/((s+1)*(s+2)*(s-1));

G_sym = tf2sym(G);
zeroRectangleTf(G_sym)
\end{lstlisting}
\begin{lstlisting}[style=matlabout]
ans=1
\end{lstlisting}


\section{System Perturbation}

\subsection{Internal stability}
Internal stabilty of a closed loop system is given if, 
\begin{equation}
    \begin{bmatrix}
        I & K(s) \\ -G(s) & I
    \end{bmatrix}^{-1} =
    \begin{bmatrix}
        (I + K(s)G(s))^{-1} & -K(s)(I + G(s)K(s))^{-1} \\ 
        G(s)(I + K(s)G(s))^{-1} & (I + G(s)K(s))^{-1}
    \end{bmatrix}
\end{equation}  

all four matrices in the LHS are stable ($RH_\infty$). Now in the given example, only 2 input and 2 output system is considered.

\begin{lstlisting}
s = tf('s');
W = [7 8; 6 7];
N = [1/(s+1) 0; 0 2/(s+2)];
G = W*N*inv(W);
K = -eye(2);
\end{lstlisting}

So, here all the four matrices need not be checkedz since there are no hidden modes or unstable pole-zero cancellations. Only $$(I + K(s)G(s))^{-1}$$ needs to be checked for stability. 