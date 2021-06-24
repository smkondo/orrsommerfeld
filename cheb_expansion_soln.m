function val = cheb_expansion_soln(y,a)
    %function that evaluates the Chebyshev expansion for given coeffs
    % y = points to evaluate Chebyshev expansion at
    % a = vector of coeffs for the Chebyshev expansion
    val = zeros(length(y),1);
    for i=1:length(y)
        for j=1:length(a)
            val(i) = val(i) + a(j)*cheb_basis(y(i),j-1);
        end
    end
end