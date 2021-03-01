# Set parameters
@with_kw mutable struct SirePars
	# No. of sires, dams per sire, and offspring per pair, and groups
	nsire::Int = 10
	dpsire::Int = 5
	oppair::Int = 2
	ngroups::Int = 2

	# Derived population numbers
	ndam::Int = nsire * dpsire
	nparent::Int = nsire + nsire * dpsire
	nprogeny::Int = nsire * dpsire * oppair
	nind::Int = nprogeny

	# TODO: I think it should be this, but we'll keep it as is for the moment
	# nind = nparent + nprogeny

	# TODO: what if this isn't an integer?
	group_size::Int = nprogeny ÷ ngroups

	# Traits
	traitnames = "susceptibility", "infectivity", "recoverability"
	# No. of traits to simulate.
	ntraits = length(traitnames)

	# Trait means
	tmeans = fill(0, ntraits)
	# Trait additive genetic variance
	tVA = [0.5, 0.1, 0.1]
	# Trait non-additive (environmental) variance
	# with no association to the pedigree
	tVE = [0.3, 0.3, 0.3]
	# Trait non-additive correlation matrix
	trhoE = I(ntraits)
	# Trait additive genetic correlation matrix.
	# In this case correlations are zero.
	trhoG = I(ntraits)

	# Infection and recoverability, targeting R0 = 2.5?
	PairwiseBeta = 0.1
	RRgamma = 0.2

	# MCMC settings
	nsample = 5000
	burnin = 500

	# Additional options
	trim_parents = false
	DEBUG = false
end #params


function summarise_pars(sp::SirePars)
	println("Generating data with:")
	println(
		sp.nsire, " sires, ",
		sp.ndam, " dams, ",
		sp.nprogeny, " progeny, (",
		sp.nind, " individuals), ",
		sp.ngroups, " groups")
	println("R0 = ", sp.PairwiseBeta / sp.RRgamma, "?")

	nothing
end

