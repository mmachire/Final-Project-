#function for StrucVar

function Struc_var(genesofint,vcf,num_of_var,how_far)
	chr = convert(Array,genesofint.a)
	posit = convert(Array,genesofint.b)
	reader = open(GeneticVariation.VCF.Reader,vcf)
	chromosome = Array{String}(undef,num_of_var)
	position = Array{Int64}(undef,num_of_var)
	a = length(chr)
	i = 0
	k = 0
	for record in reader 
		for inf = 1:a
			global k = 0
			for pos = posit[inf]:posit[inf]+100
				if pos == GeneticVariation.VCF.pos(record)
					if chr[inf] == 	GeneticVariation.VCF.chrom(record)
						i = i+1
						chromosome[i] = chr[inf]
						position[i] = pos
					end 
				end 
				if i == num_of_var
					break
				end 
			end 
			if i == num_of_var                                                                                    break                                                                                end  
		end 
		if i == num_of_var                                                                                    break                                                                                end  
	end 	
	df = DataFrame(a = chromosome, b = position)
	return df
end 

