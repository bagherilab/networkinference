% -------------------------------------------------------------------------
% GET_CHAOS calculates the chaos metric form the given matrix normalized to
% the value of the chaos metric calculated from checkerboard matrix of the
% same size.
% -------------------------------------------------------------------------

function c = get_chaos(A)

n = length(A);
i = 2:n;
j = 1:n-1;
c = zeros(n);

c(i,:) = c(i,:) + abs(A(i,:) - A(j,:)); % top
c(j,:) = c(j,:) + abs(A(j,:) - A(i,:)); % bottom
c(:,i) = c(:,i) + abs(A(:,i) - A(:,j)); % left
c(:,j) = c(:,j) + abs(A(:,j) - A(:,i)); % right

c(i,i) = c(i,i) + abs(A(i,i) - A(j,j)); % top left
c(j,j) = c(j,j) + abs(A(j,j) - A(i,i)); % bottom right
c(i,j) = c(i,j) + abs(A(i,j) - A(j,i)); % top right
c(j,i) = c(j,i) + abs(A(j,i) - A(i,j)); % bottom left

c = sum(sum(c))/(4*(n - 2)^2 + 3*(n - 2)*4 + 4*2);

end