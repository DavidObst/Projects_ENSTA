function g = vraisemblance(q_star,q_k,lambda)
	%D = distance(q_star,q_k)
	R = sum(sqrt(q_star .* q_k));
	D = sqrt(1-R);
	g = exp(-lambda*D^2);
end