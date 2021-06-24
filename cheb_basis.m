function val = cheb_basis(y,n)
    %Returns value of Chebyshev polynomial of order n evaluated at y
    val = cos(n*acos(y));
end