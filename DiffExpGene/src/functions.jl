#functions for gene count 

function Gene_Count_Matrix(gff3, bam...)

nsam = length(bam)
num_feat = 0
features = open(collect, GenomicFeatures.GFF3.Reader,gff3)
filter!(x -> GenomicFeatures.GFF3.featuretype(x) == "mRNA", features)
for feature in features
        global num_feat = num_feat+1
end

genecount = zeros(num_feat,nsam)
for i = 1:nsam
	j = 0
	reader = open(XAM.BAM.Reader,bam(i))
	for feature in features
		global j = j+1
		for record in BioAlignments.eachoverlap(reader,feature)
			global genecount(j,i) = genecount(j,i)+1
		end 
	end 
end 

#normalize
genecount_n = log.(transpose(genecount))
return genecount_n
end 

function Genes_of_Interest(gff3, genecount,n)
	Result = NMF.nnmf(genecount,n)
	H = Result.H
	#get 10 genes of interest 
	indicies = zeros(10,n)
	for i = 1:n
		X = [n,:]
		ind = sortperm(X)
		y = length(ind)
		topten = ind[y-9:y]
		global indicies[:,i] = topten
	end 
	indi = convert(Array, indicies)
	num_feat = 0
	features = open(collect, GenomicFeatures.GFF3.Reader,gff3)
	filter!(x -> GenomicFeatures.GFF3.featuretype(x) == "mRNA", features)
	for feature in features 
               global num_feat = num_feat+1  
        end 	
	
	A = Array{String}(undef,num_feat)
	B = Array{Int}(undef,num_feat)
	j = 0
	for feature in features
        	global j = j+1
        	global B[j] = GenomicFeatures.GFF3.seqstart(feature)
        	global  A[j] = GenomicFeatures.GFF3.seqid(feature)

	end

	ng = length(indi)
	C = Array{String}(undef,ng)
	D = Array{Int}(undef,ng)

	k = 0
	for i in 1:ng
		global k = k+1
		b = convert(Int64,indi[i])
		C[k] = A[b]
        	D[k] = B[b]
       	end
	df = DataFrame(a = C, b = D)
	return df
end 






	
	




