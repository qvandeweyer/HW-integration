
# module HW_int


	# question 1 b) 
	# here are the packages I used

	using FastGaussQuadrature
	using Roots
	using Sobol
	using PyPlot
	using Distributions

	# here are some functions I defined for useage 
	# in several sub questions

	# demand function

	# gauss-legendre adjustment factors for map change

	# eqm condition for question 2
	# this is the equilibrium condition: total demand = supply, 


	# weighted sum for integration from the slides.



	# function question_1b(n)

	# end


	# function question_1c(n)

	# end


	# function question_1d(n)

	# end


	# function question_2a(n)

	# end


	# function question_2b(n)

	# end	


	# function to run all questions
	# function runall(n=10)
	# 	println("running all questions of HW-integration:")
	# 	println("results of question 1:")
	# 	question_1b(n)	# make sure your function prints some kind of result!
	# 	question_1c(n)
	# 	question_1d(n)
	# 	println("")
	# 	println("results of question 2:")
	# 	q2 = question_2a(n)
	# 	println(q2)
	# 	q2b = question_2b(n)
	# 	println(q2b)
	# 	println("end of HW-integration")
	# end
	###############################################################
	###############################################################
	
	println("Collaborators: @jlkeene, @aelgouacem, @jilnur, @qvandeweyer")

	# Question 1

	# demand function
	function demandFun(P)

		2*exp(-0.5*log(P))

	end

	# price vector
	P = 0.5:10;
	# demand
	Q = demandFun(P);
	Q1 = ones(size(P))*demandFun(1); 
	Q4 = ones(size(P))*demandFun(4); 
	plot(P, [Q Q1 Q4])


	# number of nodes for approximation
	np = [10 100 1000]
	# Question 1.a : Gauss Legendre Integration
	function GLIntegration(NP::Int, LB::Real, UB::Real)
		N, W = gausslegendre(NP)
		newN = (UB-LB)/2*N+(UB+LB)/2
		CSApprox = dot(W, 1.5*demandFun(newN))

		return CSApprox
	end

	CSApproxGL = GLIntegration(np[1], 1,4);
	println("GL Integration Solution with 10 nodes = ", CSApproxGL)
	CSApproxGL = GLIntegration(np[2], 1,4);
	println("GL Integration Solution with 100 nodes = ", CSApproxGL)
	CSApproxGL = GLIntegration(np[3], 1,4);
	println("GL Integration Solution with 1000 nodes = ", CSApproxGL)

	# srand(12345);
	 
	# Question 1.b: Monte Carlo Integration
	function MCIntegration(NP::Int, LB::Real, UB::Real)
		srand(12345)
		N = rand(NP,1);
		# rescale shift by 1 and multiply by 3
		newN = (UB-LB)*N+1;
		CSApprox = sum(demandFun(newN)*(UB-LB))/NP

		return CSApprox
	end

	CSApproxMC = MCIntegration(np[1], 1,4);
	println("MC Integration Solution with 10 nodes = ", CSApproxMC)
	CSApproxMC = MCIntegration(np[2], 1,4);
	println("MC Integration Solution with 100 nodes = ", CSApproxMC)
	CSApproxMC = MCIntegration(np[3], 1,4);
	println("MC Integration Solution with 1000 nodes = ", CSApproxMC)
	CSApproxMC = MCIntegration(10000, 1, 4);
	println("MC Integration Solution with 10000 nodes = ", CSApproxMC)


	# Question 1.c: Quasi Monte Carlo Integration
	function QMCIntegration(NP::Int, LB::Real, UB::Real)
		s = SobolSeq(1,LB,UB);
		N = hcat([next(s) for i = 1:NP]...)'
		
		CSApprox = sum(demandFun(N)*(UB-LB))/NP

		return CSApprox
	end

	CSApproxQMC = QMCIntegration(np[1],1,4);
	println("QMC Integration Solution with 10 nodes = ", CSApproxQMC)
	CSApproxQMC = QMCIntegration(np[2],1,4);
	println("QMC Integration Solution with 100 nodes = ", CSApproxQMC)
	CSApproxQMC = QMCIntegration(np[3],1,4);
	println("QMC Integration Solution with 1000 nodes = ", CSApproxQMC)



	###############################################################
	###############################################################
	# Question 2.a: Gauss Hermite Expected Value and Variance

	logTheta1, weights1 = gausshermite(10);
	logTheta2, weights2 = gausshermite(10);

	stdLogTheta1 	= sqrt(0.02);
	stdLogTheta2 	= sqrt(0.01);
	covLogThetas 	= 0.01;
	meanLogThetas 	= [0 0];
	SigmaLogThetas  = [0.02 covLogThetas; 0.01 covLogThetas];
	Omega 			= chol(SigmaLogThetas,Val{:U})

	tildeTheta1 = zeros(size(logTheta1))
	tildeTheta2 = zeros(size(logTheta2))

	for i = 1:length(logTheta1)
			a = Omega*[logTheta1[i] logTheta2[i]]';
			tildeTheta1[i] = a[1];
			tildeTheta2[i] = a[2];
	end


	prc = zeros(length(tildeTheta1),length(tildeTheta2))
	for i = 1: length(tildeTheta1)
		for j = 1: length(tildeTheta1)

			prcFun(p) = log(exp(tildeTheta1[i])/p + exp(tildeTheta2[j])/p^(0.5))-log(2);
			prc[i, j] = fzero(prcFun, [0,10])  

		end
	end

	kronWeights = kron(weights1, weights2');

	expectedPrice = exp(pi^(-1.5)*sum(prc.*kronWeights));
	variancePrice = sum((prc[:] - expectedPrice).^2)/length(prc[:]);
	println("Expected Value of the Price = ", expectedPrice)
	println("Variance of the Price = ", variancePrice)

	# Question 2.b: Monte Carlo Integration 

	






# end

