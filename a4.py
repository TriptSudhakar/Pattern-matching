import random
import math

#To generate random prime less than N
def randPrime(N):
	primes = []
	for q in range(2,N+1):
		if(isPrime(q)):
			primes.append(q)
	return primes[random.randint(0,len(primes)-1)]

# To check if a number is prime
def isPrime(q):
	if(q > 1):
		for i in range(2, int(math.sqrt(q)) + 1):
			if (q % i == 0):
				return False
		return True
	else:
		return False
		
#pattern matching
def randPatternMatch(eps,p,x):
	N = findN(eps,len(p))
	q = randPrime(N)
	return modPatternMatch(q,p,x)

#pattern matching with wildcard
def randPatternMatchWildcard(eps,p,x):
	N = findN(eps,len(p))
	q = randPrime(N)
	return modPatternMatchWildcard(q,p,x)

# return appropriate N that satisfies the error bounds
def findN(eps,m):
	threshold = (2*m*math.log2(26)/eps)**2
	return max(8,math.ceil(threshold))

	# let A,B be substrings such that A!=B and h(A) = f(A) (mod q) = f(B) (mod q) = h(B)
	# let f(A) - f(B) be d
	# then h(A) = h(B) if q is one of the prime factors of d
	# we know that prime factors of a number d are at most log2(d)
	# let N be the upper bound for the random prime chosen
	# then probability of choosing such q is log2(d)/pi(N) where pi(N) is the number of primes less than or equal to N
	# we know that pi(N) >= N/2log2(N)
	# d is strictly less than 26^m
	# so the probability of choosing q such that h(A) = h(B) is strictly less than log2(26^m)/(N/2log2(N))
	# this probability should be less than less than or equal to eps
	# so m*log2(26)/(N/2log2(N))<=eps
	# => N/log2(N) >= 2*(m/eps)*log2(26)
	# for N >= 8 , sqrt(N)>log2(N)
	# which implies N/log2(N)> sqrt(N) for N>=8
	# so if we take sqrt(N) >= 2*(m/eps)*log2(26), then our condition for N is stil satisified
	# => N = (2*(m/eps)*log2(26))^2

# Return sorted list of starting indices where p matches x
def modPatternMatch(q,p,x):
	output = [] # initialization of output list
	m = len(p) # length of the pattern
	n = len(x) # length of the string given for pattern matching 
	hash_value = 0 # hash_value stores the (value of the pattern interpreted as a 26-ary number) modulo q where q is the prime number used for hashing
	mod_exp = 1 # mod_exp stores the value of 26^m modulo q

	for i in range(m):
		hash_value = (hash_value*26 + ord(p[i]) - ord('A')) % q
		mod_exp = (mod_exp*26)%q

	hash_iter = 0 # hash_iter stores the (value of the hash function) modulo q corresponding to any substring of length m of the string x 
	for i in range(m):
		hash_iter = (hash_iter*26 + ord(x[i]) - ord('A')) % q

	for i in range(n-m+1):
		if hash_iter==hash_value: 
			output.append(i) # if value of hash function modulo q is same for the pattern and the substring of length m starting at index i

		if i+m==n: break

		hash_iter = (26*hash_iter + ord(x[i+m]) - ord('A') - mod_exp*(ord(x[i]) - ord('A'))) % q 
		# hash_iter is now assigned the value of hash function modulo q for next substring of length m starting at index (i+1)
	return output

	# Time complexity:
	# Computation of hash_value takes O(mlog2(q)) time
	# Computation of hash_iter for all possible susbstrings of length m takes O(nlog2(q)) time 
	# So time complexity is given by O((m+m)log2(q))


	# Space complexity:
	# Storing output list of size k takes O(k) space
	# Storing indices i,j takes O(log2(n)) space
	# Storing hash_value and hash_iter takes O(log2(q)) space
	# So space complexity is given by O(k+log2(q)+log2(n))

# Return sorted list of starting indices where p matches x
def modPatternMatchWildcard(q,p,x):
	output = [] # initialization of output list
	m = len(p) # length of the pattern
	n = len(x) # length of the string given for pattern matching 
	mod_exp = 1 # mod_exp stores the value of 26^m modulo q
	hash_value = 0 # hash_value stores the (value of the pattern interpreted as a 26-ary number) modulo q where q is the prime number used for hashing
	wildcard = 0 # index of the wildcard character in the pattern 
	mod_wildcard = 1 # mod_exp stores the value of 26^(m-wildcard-1) modulo q

	for i in range(m):
		mod_exp = (mod_exp*26) % q
		if ord(p[i])>ord('Z') or ord (p[i])<ord('A'):
			wildcard = i # store the index of the wildcard character 
			# the value of  wildcard character is treated as 0
			hash_value = (hash_value*26) % q
		else:
			hash_value = (hash_value*26 + ord(p[i]) - ord('A')) % q

	# To match the pattern with the wildcard character, while calculating the hash function modulo q of the substring of x of length m,
	# we ignore the value of the character (take it 0) at the index correspinding to index of the wildcard in the pattern
	# To do this we take 2 hash iterables which store the value of the hash function modulo q in 2 parts such that 
	# value of the modified hash function of the substring of x is given by (hash_iter1 + hash_iter2) modulo q
	hash_iter1 = hash_iter2 = 0

	for i in range(wildcard):
		hash_iter1 = (hash_iter1*26 + ord(x[i]) - ord('A')) % q
	hash_iter1 = (hash_iter1*26)%q

	for i in range(wildcard+1,m):
		hash_iter2 = (hash_iter2*26 + ord(x[i]) - ord('A')) % q
		hash_iter1 = (hash_iter1*26) % q
		mod_wildcard = (mod_wildcard*26) % q

	i = 0
	for j in range(wildcard+1,n-m+wildcard+2):
		if (hash_iter1+hash_iter2)%q==hash_value:
			output.append(i) # if value of modified hash function modulo q is same for the pattern and the substring of length m starting at index i

		if j+m==n+wildcard+1: break

		# Updating the hash iterables
		hash_iter1 = (26*hash_iter1 + ((26*mod_wildcard)%q)*(ord(x[i+wildcard]) - ord('A')) - (mod_exp*(ord(x[i]) - ord('A')))%q) % q
		hash_iter2 = (26*hash_iter2 + ord(x[j+m-wildcard-1]) - ord('A') - mod_wildcard*(ord(x[j]) - ord('A'))) % q
		i+=1
	
	return output

	# Time complexity:
	# Computation of hash_value takes O(mlog2(q)) time
	# Computation of hash_iter1 and hash_iter1 for all possible susbstrings of length m takes O(nlog2(q)) time 
	# So time complexity is given by O((m+m)log2(q))

	# Space complexity:
	# Storing output list of size k takes O(k) space
	# Storing indices i,j takes O(log2(n)) space
	# Storing hash_value, hash_iter1 and hash_iter2 takes O(log2(q)) space
	# So space complexity is given by O(k+log2(q)+log2(n))

# print(randPatternMatchWildcard(0.05, 'JA?S', 'AMADBOXERSHOTQUICKGLOVEDJABSTOTHEJAWSOFHISDIZZYOPPONENTATTHEJAMSROCKSHUFFLE'))