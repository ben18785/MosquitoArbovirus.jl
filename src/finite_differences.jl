module finite_differences
export second_derivative

"""
Construct the finite difference form of the second spatial derivative with
user-dictated boundary conditions at ``x=1``.

We now describe how the boundary conditions enter into the matrix. At ``x=0``,
we assume the flux of infection into the layer equals that out of it. So we
have the flux conservation equation:

``incoming flux = - D \\frac{\\partial m}{\\partial x}|_{x=0}``

Letting the LHS equal ``\\lambda(t)`` and using the central difference formula
for the RHS spatial derivative, we obtain:

``\\lambda(t) = -D \\frac{u_1 - u_{-1}}{\\Delta x},``

which can be rearranged to:

``u_{-1} = u_1 + \\Delta x \\frac{\\lambda(t)}{D(t)}.``

We then obtain an expression for the ODE governing this boundary:

``\\frac{d u_0}{d t} = D \\frac{u_{1} - 2 u_0 + u_{-1}}{\\Delta x^2} = D \\frac{2 u_1 - 2 m_0}{\\Delta x^2} + \\frac{\\lambda(t)}{\\Delta x}``

The latter term above does not enter the Laplacian matrix and is handled in the ODE RHS.

Considering now the ``x=1`` boundary, two possibilities are available to the user:
either "neumann" or "dirichlet". For the neumann case, this is "no flux":

``\\frac{\\partial u(x,t)}{\\partial t}|_{x=1} = 0``

then using the central finite difference approximation, we obtain:

``\\frac{u_{T+1} - u_{T-1}}{\\Delta x} = 0,``

which implies:

``u_{T+1}=u_{T-1}.``

So the method of lines (MOL) time derivative for ``u_T`` is given by:

``\\frac{d u_T}{d t} = D \\frac{u_{T+1} - 2 u_T + u_{T-1}}{\\Delta x^2} = D \\frac{2 u_{T-1} - 2 u_T}{\\Delta x^2}.``

For the Dirichlet case, we assume that the quantity is pinned at a given value,
``dirichlet_val`` which we'll shorten to ``dv`` in what follows.

Thus ``u_T=dv`` so ``\\frac{d u_T}{d t}=0``. The time-derivative for ``u_{T-1}`` is:

``\\frac{d u_{T-1}}{d t} = D \\frac{u_{T} - 2 u_{T-1} + u_{T-2}}{\\Delta x^2}= D \\frac{- 2 u_{T-1} + u_{T-2}}{\\Delta x^2} + D \\frac{dv}{\\Delta_x^2}``

In the above, the last term on the RHS is not included in the second derivative matrix.

"""
function second_derivative(delta_x, x_max, bc_rhs="neumann")
    x = 0:delta_x:x_max
    N = length(x)
    A = zeros(N, N)
    for i in 1:N, j in 1:N
        if abs(i - j) == 1
            A[i, j] = 1
        end
        if i == j
            A[i, j] = -2
        end
    end

    # lhs flux conversation boundary
    A[1, 2] = 2

    if bc_rhs == "neumann"
        A[N, N - 1] = 2
    end

    if bc_rhs == "dirichlet"
        # dirichlet
        A[N - 1, N] = 0.0

        # ensuring that time-derivate at x=1 is zero
        A[N, N - 1] = 0
        A[N, N] = 0
    end

    A = A / (delta_x^2);
    A
end

end
