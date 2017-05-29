library(methods)
library(optparse)
library(parallel)
library(Rcpp)
library(RcppArmadillo)


sourceCpp('~/break_composition/break_composition.cpp', cacheDir='~/break_composition')


get_input = function(){
    
    option_list = list(        
        make_option('--data', help='Input matrix'),
	make_option('--dim', help='Dimension to use', default='rows'),
	make_option('--samples', help='Samples to use', default=''),
	make_option('--vars', help='Variables to use', default=''),
	make_option('--minu', help='Minimum number of unique samples per variable', default=10, type='integer'),
	make_option('--cores', help='(P) Number of cores', default=1, type='integer'),
	make_option('--max_size', help='(P) Maximum samples per core', default=1e5, type='integer'),
	make_option('--contig', help='(P) Divide samples into contiguous blocks?', default=FALSE, action='store_true'),
	make_option('--print', help='Print commands? (for cluster)', default=FALSE, action='store_true'),
	make_option('--out', help='Output prefix', default=''),
	make_option('--merge', help='Merge regex', default=''),
	make_option('--norm', help='Normalize samples?', default=FALSE, action='store_true')
    )
    parse_args(OptionParser(option_list=option_list))
}

print_command = function(data, dim, samples, vars, minu, cores, max_size, contig, out){
    
    for(i in 1:length(samples)){
	indices = paste(samples[[i]], collapse=',')
	outi = paste0(out, '.', i)
	cmd = paste('Rscript ~/break_composition/break_composition.r --data', data, '--dim', dim, '--samples', indices, '--minu', minu, '--cores 1 --max_size', max_size, '--contig', contig, '--out', outi)
	if(vars != ''){
	    cmd = paste(cmd, '--vars', vars)
	}
	cat(paste0(cmd, '\n'))
    }
}


break_composition = function(x, dim='rows', samples='', vars='', minu=10, p.cores=1, max_size=1e10, p.contig=TRUE, print=FALSE, out='', norm=FALSE){

    x = subset_data(x, dim=dim, samples=samples, vars=vars, norm=norm)
    
    indices = get_sample_indices(nrow(x), max_size=max_size, p.cores=p.cores, p.contig=p.contig)

    if(print == TRUE){
        print_command(args$data, dim, indices, vars, minu, p.cores, max_size, p.contig, args$out)
	return(NULL)
    }
		  
    if(p.cores == 1){
	
	res = lapply(indices, function(i){
	    
	    # Remove non-varying columns
	    j = apply(x[i,], 2, function(a){length(unique(a))}) >= minu
	    
	    # Objective function
            f = function(s){totCorC(s*x[i,j])}
	    
	    # Estimate size factors for subset
            optim(rep(1, length(i)), f, method='L-BFGS-B', lower=.1, upper=100)$par
	})
	
    } else {

	res = mclapply(indices, function(i){
	    
	    # Remove non-varying columns
	    j = apply(x[i,], 2, function(a){length(unique(a))}) >= minu
	    
	    # Objective function
	    f = function(s){totCorC(s*x[i,j])}
	    
	    # Estimate size factors for subset
	    optim(rep(1, length(i)), f, method='L-BFGS-B', lower=.1, upper=100)$par
	}, mc.cores=p.cores)	
    }
    sizes = data.frame(size = unlist(res), row.names=rownames(x)[unlist(indices)])
    sizes = sizes[rownames(x),,drop=F]
    
    if(out != ''){
        write.table(sizes, file=paste0(out, '.sizes.txt'), sep='\t', quote=F, col.names=F)
    }
    
    return(sizes)
}


subset_data = function(x, dim='rows', samples='', vars='', minu=10, norm=FALSE){
    
    # Normalize samples
    if(norm == TRUE){
        if(dim == 'rows'){
            x = t(scale(t(x), center=F, scale=rowSums(x)))
	} else {
	    x = scale(x, center=F, scale=colSums(x))
	}
    }
    
    # Convert to matrix
    x = as.matrix(x)
    
    # Transpose
    if(dim == 'cols'){
        x = t(x)
    }
    
    # Subset rows
    if(samples != ''){
        i = as.numeric(strsplit(as.character(samples), ',')[[1]])
	if(length(i) == 1){
	    i = sample(1:nrow(x), i)
	}
	x = x[i,]
    }
    
    # Subset columns
    if(vars != ''){
        j = as.numeric(strsplit(as.character(vars), ',')[[1]])
	if(length(j) == 1){
	    j = order(colMeans(x), decreasing=T)[1:j]
	}
	x = x[,j]
    }
    
    # Remove non-varying columns
    j = apply(x, 2, function(a){length(unique(a))}) >= minu
    x = x[,j]

    return(x)
}


get_sample_indices = function(n, p.cores=1, max_size=1e10, p.contig=FALSE){

    # Calculate number of parallel tasks
    num_iter = max(p.cores, ceiling(n/max_size))
    
    # Split sample indices into equally sized groups
    if(p.contig == TRUE){
        i = 1:n
    } else {
        i = sample(1:n, n)
    }
    split(i, ceiling(seq_along(1:n)/(n/num_iter)))
}


merge_files = function(regex, out){

    # Get mean sizes
    files = list.files(pattern=regex)
    sizes = lapply(files, function(a){read.table(a, sep='\t', stringsAsFactors=F)})
    sizes = do.call(rbind, sizes)
    sizes = tapply(sizes[,2], list(sizes[,1]), mean)
    sizes = data.frame(sizes)

    # Write output
    write.table(sizes, file=paste0(out, '.sizes.merged.txt'), sep='\t', quote=F)
}    


if(!interactive()){
    
    # Get command line arguments
    args = get_input()

    if(args$merge != ''){
        
	merge_files(regex=args$merge, out=args$out)
        
    } else {
    
        # Read data
        x = read.table(args$data, sep='\t', header=T, row.names=1)
    
        # Break composition
        sizes = break_composition(x, dim=args$dim, samples=args$samples, vars=args$vars, minu=args$minu, p.cores=args$cores, max_size=args$max_size, p.contig=args$contig, print=args$print, out=args$out, norm=args$norm)
    
    }
} else {

    args = list(data='test', dim='rows', samples='', vars='', minu=10, p.cores=1, max_size=Inf, p.contig=FALSE, print=FALSE, out='test', norm=FALSE)
}
