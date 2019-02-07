function ptr=reconstructPoints(B, ptGiven)
  %This function reconstructs a set of points, X, from the product B=X*X'
  %Given where the first 3 points should be, it will also rotate the
  %reconstructed points so that they are in the same orientation as the
  %original set of points.
  %
  %Input:
  % B= product X*X' where X is the set of points to be reconstructed 
  %
  % ptGiven = 3x3 matrix containing the location of the first 3 points from
  % the original data
  %
  %Output:
  %
  % ptr= reconstructed points
[v,s]=eigs(B,3);
n=length(B);
ptr=v*sqrt(s);
if norm(imag(ptr), 'fro')>sqrt(eps)
    warning('Point reconstruction contains imaginary points.  Function will return real part')
end
ptr=real(ptr);

%rotationMatrix=inv(ptr(1:3,:))*ptGiven;
%ptr=ptr*rotationMatrix;

end