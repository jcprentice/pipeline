using Parameters
using DataFrames

function ped_grm(sp::SirePars, pedigree = nothing)
	@unpack ntotal, nsire, ndam, ngroup, nparent, nprogeny, dpsire, oppair = sp

	if pedigree == nothing
		ped = DataFrame(id = 1:ntotal
						sire = 0, dam = 0, group = 0, dso = :none)

		ped.dso  .= [fill(:sire, nsire); fill(:dam, ndam); fill(:offspring, nprogeny)]

		prog_ids = (1:nprogeny) .+ nparent

		ped.sire[prog_ids]  .= repeat(1:nsire,           inner = dpsire * oppair)
		ped.dam[prog_ids]   .= repeat((1:ndam) .+ nsire, inner = oppair)
		ped.group[prog_ids] .= repeat(1:ngroup,          outer = nprogeny ÷ ngroup)
	else
		ped = pedigree
	end #if

	ped2 = copy(ped)


	# construct the GRM
	GRM = zeros(ntotal, ntotal)

	for i in 1:ntotal, j in i + 1:ntotal
		if ped.sire[j] == i || ped.dam[j] == i
			# i is j's parent
			GRM[i, j] = 1/2
		elseif ped.dso[i] != :offspring || ped.dso[j] != :offspring
			# not both offspring, at least one is a parent
		else
			# check if parents match
			same_sire = ped.sire[i] == ped.sire[j]
			same_dam = ped.dam[i] == ped.dam[j]

			if same_sire && same_dam
				# both parents match, full sibs
				GRM = 1/2
			elseif same_sire || same_dame
				# one parent matches, half sibs
				GRM = 1/4
			end #if
		end #if
	end #for

	# combine with transpose and identity
	GRM .= GRM + I + GRM'


	ped, GRM
end #function
