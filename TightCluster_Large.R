###############################################################################
###############################################################################
#	FUNCTIONS : PartitionClust(I), f(II), g(III), 	 delrow(IV)				#
#			merger(V), 		 FinalClust(VI), 	 tightClust(VII)			#
#	Ext FUNS  : tight.clust(IE) 	(tightClust package),					#
#			spantree(IIE) 	(vegan 	package)								#
#	FUNS MAP  : VII --> I, VI												#
#			I   --> IE														#
#			VI  --> V, II, III												#
#																			#
#	USE CASE  : x = cbind(rmvnorm(1000,c(0,0)), rmvnorm(1000,c(10,10)))		#
#			tightClust(x,2)													#
###############################################################################
######################MAIN FUNC DATAILS########################################
#tightClust (VII)															#
#Input 	: x = data matrix   , target  = wanted no of clusters.				#
#		  K = number of reps, maxPart = max allowed size of parts			#
#		  alpha, beta inputs for tight.clust								#
#																			#
#target 	: no of clusters finally appearing in the variable "Cluster"	#
#																			#
#Output	: Cluster = a list of clusters, size = dim of the clusters,			#
#		  CsLabel = class labels											#
###############################################################################


############UPDATE ############################################################
#Fix error 
#	“Error in mst(x[a, ]) : INTEGER() can only be applied to a 'integer', not a 'NULL'”  #
#	Changes made in functions PartitionClust (lines 270 - 273)
#							  f(lines 291 - 297), merger (lines 427)
#								g ( line 356 ), tightClust(lines 450)
##############################################################################
#library("tightClust")
#library("vegan")
#library("foreach")
#library('rpud')
library("nnclust")

tight.clust <-
function (x, target, k.min, alpha = 0.1, beta = 0.6, top.can = 7,
seq.num = 2, resamp.num = 10, samp.p = 0.7, nstart = 1, remain.p = 0.1,
k.stop = 5, standardize.gene = TRUE, random.seed = NULL)
{
    find.candidates <- function(x, k, alpha = 0.1, top.can = 7,
    resamp.num = 10, samp.p = 0.7, nstart = 1) {
        kmeans.classify <- function(x, k, p = 0.7, ...) {
            size <- dim(x)[1]
            partial.size <- round(p * size)
            partial <- sample(rownames(x), partial.size)
            centers <- kmeans(x[partial, ], k, ...)$centers
            return(suppressWarnings(kmeans(x, centers, iter.max = 1,
            algorithm = "Lloyd"))$cluster)
        }
        find.candidates.one <- function(x) {
            tmp <- apply(x == 1, 1, sum)
            return(which(x[, which(tmp == max(tmp))[1]] == 1))
        }
        extend.candidate <- function(D, can, alpha = 0.1) {
            can.ex <- which(apply(as.matrix(D[, can] >= 1 - alpha),
            1, all))
            D.temp <- D[can.ex, can.ex]
            if (!is.matrix(D.temp)) {
                D.temp <- as.matrix(D.temp)
                colnames(D.temp) <- names(can.ex)
            }
            D.bad <- apply(as.matrix(D.temp < 1 - alpha), 1,
            sum)
            while (sum(D.bad) > 0) {
                index <- which(D.bad == max(D.bad))[1]
                D.temp <- D.temp[-index, -index]
                D.bad <- apply(as.matrix(D.temp < 1 - alpha),
                1, sum)
            }
            return(can.ex[colnames(D.temp)])
        }
        N <- dim(x)[1]
        Dbar <- matrix(0, N, N)
        for (i in 1:resamp.num) {
            cl <- kmeans.classify(x, k, p = samp.p, nstart = nstart,
            iter.max = 100)
            D <- outer(cl, cl, function(a, b) a == b)
            Dbar = Dbar + D
        }
        Dbar = Dbar/resamp.num
        colnames(Dbar) <- 1:N
        rownames(Dbar) <- 1:N
        i = 1
        D.temp <- Dbar
        res <- list()
        while (i <= top.can * 2 && dim(D.temp)[1] > 0) {
            candidate.one <- find.candidates.one(D.temp)
            candidate <- extend.candidate(D.temp, candidate.one,
            alpha = alpha)
            D.temp <- D.temp[-candidate, -candidate]
            res[[i]] <- names(candidate)
            mode(res[[i]]) <- "numeric"
            i = i + 1
        }
        res <- res[order(unlist(lapply(res, length)), decreasing = TRUE)][1:top.can]
        return(res)
    }
    if (!is.null(random.seed))
    set.seed(random.seed)
    original.data <- x
    if (standardize.gene)
    x <- t(scale(t(x)))
    k.max <- k.min + 10
    id <- rownames(x)
    N <- dim(x)[1]
    #write(paste("Number of points:", N, "\tDimension:", dim(x)[2],
    #    "\n"), "")
    rownames(x) <- 1:N
    index.m <- as.matrix(expand.grid(lapply(1:seq.num, function(x) 1:top.can)))
    remain <- N
    nfound <- 0
    found <- TRUE
    k0 <- k.min
    k <- k0
    candidates <- list()
    tclust <- list()
    while (nfound < target && remain/N >= remain.p && (found ||
    k <= k.max)) {
        if (found) {
            #write(paste("Looking for tight cluster", nfound +
            #    1, "..."), "")
            k <- k0
            for (i in 1:seq.num) {
                #   write(paste("k =", k + i - 1), "")
                candidates[[i]] <- find.candidates(x, k + i -
                1, alpha = alpha, top.can = top.can, resamp.num = resamp.num,
                samp.p = samp.p, nstart = nstart)
            }
        }
        else {
            candidates <- candidates[-1]
            candidates[[seq.num]] <- find.candidates(x, k + seq.num -
            1, alpha = alpha, top.can = top.can, resamp.num = resamp.num,
            samp.p = samp.p, nstart = nstart)
        }
        calc.beta <- function(y) {
            temp <- lapply(1:seq.num, function(z) candidates[[z]][[y[z]]])
            i.temp <- temp[[1]]
            u.temp <- temp[[i]]
            for (j in 2:seq.num) {
                i.temp <- intersect(i.temp, temp[[j]])
                u.temp <- union(u.temp, temp[[j]])
            }
            return(length(i.temp)/length(u.temp))
        }
        beta.temp <- unlist(apply(index.m, 1, calc.beta))
        if (any(beta.temp >= beta)) {
            found = TRUE
            nfound = nfound + 1
            #write(paste(nfound, "tight cluster(s) found!"), "")
            if (k0 > k.stop)
            k0 = k0 - 1
            found.temp <- candidates[[seq.num]][[index.m[which(beta.temp >=
            beta)[1], seq.num]]]
            tclust[[nfound]] <- rownames(x)[found.temp]
            mode(tclust[[nfound]]) <- "numeric"
            x <- x[-found.temp, ]
            remain <- remain - length(tclust[[nfound]])
            #write(paste("Cluster size:", length(tclust[[nfound]]),
            #    "\tRemaining number of points:", remain, "\n"),
            #    "")
        }
        else {
            found = FALSE
            k = k + 1
        }
    }
    clust.id <- rep(-1, N)
    size <- unlist(lapply(tclust, length))
    for (i in 1:length(tclust)) clust.id[tclust[[i]]] <- i
    res <- list(data = original.data, cluster = clust.id, size = size)
    class(res) <- "tight.clust"
    return(res)
}

mst<-function (X, rebuild = sqrt(nrow(X))/4)
{
    redo <- as.integer((1:rebuild) * nrow(X)/(rebuild + 1))
    rval <- .Call("call_primq", X, redo)
    rval[[1]] <- rval[[1]][-1] + 1
    rval[[2]] <- rval[[2]][-1] + 1
    rval[[3]] <- rval[[3]][-1]
    names(rval) <- c("from", "to", "dist")
    rval$n <- NROW(X) - 1
    rval$p <- NCOL(X)
    class(rval) <- "mst"
    rval
}

###############################(I)#############################################
#function	: partitions the data in L parts and gets tight most cluster 	#
#	 	  from each part. rep it K times and put them all in bucket b		#
#Input	: x   = data matrix(last column contains the index					#
#		  b   = bucket of partitions										#
#		  L   = no of partitions											#
#		  K   = no of times the process is repeated							#
#		  col = dim of the data												#
#Output	: b   = updated bucket of clusters									#
###############################################################################

PartitionClust.new <- function(idx, b, K, L, col, alpha, beta){
	#row  = nrow(x)   #no of rows in x
	#col = ncol(x)   #no of columns in x
	len  = length(b) #no of elements in the basket
	
	foreach(i=1:K) %dopar%
	{
		# cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
		cat(">>Running for repetition number ", i,"of ",K,".\n")
		# cat(">>Partitioning data of size ", length(idx)," partitioning to ", L," parts.\n")
		a   = as.integer(length(idx)/L)
		res = length(idx)  - a*L
		#cat(">>Partition size =", a,".\n")
		p   = sample(idx) #sampling before partition

		if(res > 0){left = p[((a*L)+1):length(idx)]}

		foreach(j=1:L) %dopar%
		{
            	#cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
			#cat("Considering partition number ",j,"of Rep No.",i,"\n")
			m = p[(a*(j-1)+1):(a*j)]
			if(j < (res+1)){ m = c(m,left[j])}
			#cat(">>Pulling the tight most cluster of the part no. ", j,", no. data points = ", length(m),".\n")
			cluster = tight.clust(x[m,],1,5,alpha,beta,top.can=2)$cluster
			b[[len + L*(i-1) + j]] = m[cluster==1]
			cat(">>Part No.",j,"of ",L," part size ",length(m)," cluster size ",sum(cluster==1),"\n")
		}

	}
	return(b) #returning the basket of clusters
}


PartitionClust <- function(idx, b, K, L, col, alpha, beta){
	#row  = nrow(x)   #no of rows in x
	#col = ncol(x)   #no of columns in x
	len  = length(b) #no of elements in the basket
	
	for(i in 1:K)
	{
		# cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
		cat(">>Running for repetition number ", i,"of ",K,".\n")
		# cat(">>Partitioning data of size ", length(idx)," partitioning to ", L," parts.\n")
		a   = as.integer(length(idx)/L)
		res = length(idx)  - a*L
		#cat(">>Partition size =", a,".\n")
		p   = sample(idx) #sampling before partition

		if(res > 0){left = p[((a*L)+1):length(idx)]}

		for(j in 1:L)
		{
            	#cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
			#cat("Considering partition number ",j,"of Rep No.",i,"\n")
			m = p[(a*(j-1)+1):(a*j)]
			if(j < (res+1)){ m = c(m,left[j])}
			#cat(">>Pulling the tight most cluster of the part no. ", j,", no. data points = ", length(m),".\n")
			cluster = tight.clust(x[m,],1,5,alpha,beta,top.can=2)$cluster
			b[[len + L*(i-1) + j]] = m[cluster==1]
			cat(">>Part No.",j,"of ",L," part size ",length(m)," cluster size ",sum(cluster==1),"\n")
		}

	}
	if(len >0)
		b = c(b[1:len],b[(len+1):length(b)][which(0<sapply(b[(len+1):length(b)],FUN=length))])
	else
		b = b[which(sapply(b,FUN=length)>0)]
	return(b) #returning the basket of clusters
}

###############################(II)############################################
#f		: function to calculate distance between two clusters				#
#Inputs	: two cluster matrices a,b,dimension of gene data					#
#output	: whether a,b belong to same cluster or not							#
###############################################################################

f = function(a,b){
	#require("vegan")
	#require('nnclust')
	#calculation of d: the distance of the two matrices a and b	
	
	#D1 = dist(x[a,])
	#D2 = dist(x[b,])
	#D  = dist(x[unique(c(a,b)), ])
	if(length(a)== 0 | length(b) == 0) return(1)
	debuga <<- a
	debugb <<- b
	debugab <<- unique(c(a,b))
	MD = ifelse(length(a)==1,0,max(mst(as.matrix(x[a,]))[['dist']]))				#max(spantree(D)$dist)
	MD1= ifelse(length(b)==1,0,max(mst(as.matrix(x[b,]))[['dist']]))				#max(spantree(D1)$dist)
	MD2= ifelse(length(unique(c(a,b)))==1,0,max(mst(as.matrix(x[unique(c(a,b)),]))[['dist']]))	#max(spantree(D2)$dist)	

#	MD= mean(spantree(D)$dist)
#	MD1= mean(spantree(D1)$dist)
#	MD2= mean(spantree(D2)$dist)	
	
	MD12= (MD1+MD2)/2
	if(MD > MD12)
		{return (0) }
	else 
		{return (1)}
}

# f.new = function(a,b){
	# require("vegan")
	# #calculation of d: the distance of the two matrices a and b	
	
	# D1 = dist.matrix(x[a,])
	# D2 = dist.matrix(x[b,])
	# D  = dist.matrix(x[unique(c(a,b)), ])
	# MD = max(spantree(D)$dist)
	# MD1= max(spantree(D1)$dist)
	# MD2= max(spantree(D2)$dist)	

# #	MD= mean(spantree(D)$dist)
# #	MD1= mean(spantree(D1)$dist)
# #	MD2= mean(spantree(D2)$dist)	
	
	# MD12= (MD1+MD2)/2
	# if(MD > MD12)
		# {return (0) }
	# else 
		# {return (1)}
# }

# f.new1 = function(a,b){
	# require("vegan")
	# #calculation of d: the distance of the two matrices a and b	
	
	# D1 = rpuDist(x[a,])
	# D2 = rpuDist(x[b,])
	# D  = rpuDist(x[unique(c(a,b)), ])
	# MD = max(spantree(D)$dist)
	# MD1= max(spantree(D1)$dist)
	# MD2= max(spantree(D2)$dist)	

# #	MD= mean(spantree(D)$dist)
# #	MD1= mean(spantree(D1)$dist)
# #	MD2= mean(spantree(D2)$dist)	
	
	# MD12= (MD1+MD2)/2
	# if(MD > MD12)
		# {return (0) }
	# else 
		# {return (1)}
# }
#g=function to determine tightness of a cluster

g=function(a){
	b = ifelse(length(a)>1,mean(mst(as.matrix(a))[['dist']]),0)	#mean(spantree(dist(a))$dist)
	return(b)
}

# g.new=function(a){
	# b = mean(spantree(dist.matrix(a))$dist)
	# return(b)
# }

# g.new1=function(a){
	# b = mean(spantree(rpuDist(a))$dist)
	# return(b)
# }
###############################(IV)############################################
#delrow	: function to delete rows from the basket after tightmost			#
#		  cluster is found,delete "u" from "v"								#
#Input	: two matrices u and v												#	
#Output	: matrix v without the rows of u									#
###############################################################################

delrow=function(u,v,dimension){
	u1=unique(matrix(u,ncol=dimension))  	#u1:unique rows of u
	v1=unique(matrix(v,ncol=dimension))  	#v1:unique rows of v
	m=nrow(u1)	   					#m:total rows of u1
	n=nrow(v1)	   					#n:total rows of v1
	y=duplicated(rbind(u1,v1))
	z=y[(m+1):(m+n)]
	return(v[z==FALSE,])      			#return v without rows of u
}

###############################(V)#############################################
#merger     : function to merge elements of the bucket b based on			#
#		  dexision function FUN												#
#Input	: b   = the bucket of tight-most clusters 							#
#			  (contains the index in last column)							#
#		  col = dim of data													#
#		  FUN = returns 0 or 1 based on whether not or whether to 			#
#			  join two clusters. 											#
#			  input: a and b two matrices with (col+1) columns				#
#					   last column contains indexes							#
#				   col = no of cols - 1										#
#				   t   = the number of tight cluster 						#
#					   in order to be extracted (unused)					#
#			 output: 0 or 1													#
#Output	: Redefined bucket													#
###############################################################################



merger <- function(b, col, FUN)
{
	len = length(b)
	if(len <= 1) return(b)
	i = 1
	repeat{
		j = i+1
		repeat{
			if(j > len){break}
			if(FUN(b[[i]], b[[j]])==1){
				b[[i]] = c(b[[i]],b[[j]])
				b[[i]] = unique(b[[i]]) 
				b      = b[-j]
				j      = j-1
				len    = len-1
			}
			j = j+1
			if(j>len){break}
		}
		i = i+1
		if(i > len){break}
	}
	b = b[which(sapply(b,FUN=length)>0)]
	return(b)
}

# merger.new <- function(b, col, FUN)
# {
	# len = length(b)
	# if(len <= 1) return(b)
	
	# g <- function(a,b){
		# dec = FUN(a[[1]],b)
		# if(dec)	
			# return(list(unique(c(a[[1]],b)), c(a[[2]],dec)))
		# else return(list(a[[1]], c(a[[2]],dec)))
	# }
	
	# i = 1
	# while( i < length(b)){
		# m <- Reduce( g, b[-(1:i)], list(b[i][[1]], rep(0,i)))
		# b[i][[1]] = m[[1]]
		# b <- b[!m[[2]]]
		# i = i+1
	# }
	# return(b)
# }
###############################(VI)############################################
#FinalClust : function to chose the tightmost cluster and 					#
#		  delete the rows of tightmost cluster from basket					#
#Input	: b   = the bucket of tight-most clusters 							#
#			  (contains the index in last column)							#
#		  col = dim of data													#
#		  t   = the number of tight cluster in order to be extracted 		#
#Output	: a list with two elements											#
#		  res = the cluster extracted										#	
#		  b   = new bucket with all extracted elementes removed				#
###############################################################################

FinalClust <- function(b, col, t, len, L){
	#cat(".................................................................................\n")
	cat(">>>Identifying distinct clusters on bucket of length ", length(b),".\n")
	
	K  = (length(b) - len)/L  ###
	b1 = lapply(1:K, FUN = function(i) merger(b[(len + (i-1)*L+1):(len + i*L)], col, f) )
	b = b[-((1+len):length(b))]
	for(i in 1:K) b = c(b, b1[[i]])
	
	b = merger(b, col, f)		#merging based on cluster merger decision function 'f'

	#cat("bucket has",length(b),"element(s) after decision wise combining\n")

	FUN = function(a, b)
			(length(unique(c(a,b)))) <(length(a)*1.05)

	b = merger(b, col, FUN)		#merging based on number of common elements in two clusters
	
	# pick tightmost cluster and remove from basket, b
    	cat(">>>Picking the tight-most cluster from ",length(b)," distinct clusters.\n")
	Cl  = which.min(sapply(b, FUN = function(idx) g(x[idx,])))
	res = b[[Cl]]
	b   = b[-Cl]

#	cat(">>>A tight cluster of size ", length(res)," has been picked up from the bucket.\n")
#	cat(">>>Removing the points of the cluster from the elements of the bucket.\n")
#	cat(">>>Current bucket size is ", length(b),".\n")

	#remove elements of tightmost cluster
	if(length(b) > 0){
		for(i in 1:length(b)){
            #	cat(">>>Deleting cluster from the ", i,"th element of the bucket.\n")
			b[[i]] = setdiff(b[[i]], res)			
		}
	}

	return(list(b=b,res=res)) 
}

###############################(VII)###########################################
#main function																#
#Input 	: x = data matrix   , target  = wanted no of clusters.				#
#		  K = number of reps, maxPart = max allowed size of parts			#
#		  alpha, beta inputs for tight.clust								#
#																			#
#target 	: no of clusters finally appearing in the variable "Cluster"	#
#																			#
#Output	: Cluster = a list of clusters, size = dim of the clusters,			#
#		  CsLabel = class labels											#
###############################################################################

tightClust <- function(x, target, K, alpha = 0.1, beta = 0.6, maxPart = 2000){
	x 			<<- x
	b 			= list()    		#create a basket for clusters
	Cluster 	= list() 			#create list for final clusters
	size		= list()
	row 		= nrow(x)   		#no of rows of x
	col 		= ncol(x)   		#no of columns of x
	#x		= cbind(x, 1:row)		#adding the index colum
	#data_m 	= x         		#store x as a new variable
	index 	= rep(0,row)		#matrix(0,row,1)
	idx		= 1:row

	#every iteration gives a tight-most cluster	
	
	for (t in 1:target)
	{
		L = ceiling(length(idx)/maxPart)	#Allow atmax 2000 points in each part

		#cat("---------------------------------------------------------------------------\n")
		#cat(">Looking for ", t," tight cluster with ",dim(x)[1]," data points,\n existing bucket size = ", length(b),", K = ",K,", L = ", L,"\n")
		len = length(b) ###
		b = PartitionClust(idx, b, K, L, col, alpha, beta)
		
		if(length(idx)<30){break}		#break if very few data points are left
		
		Rt  = FinalClust(b, col, t, len, L)  	#picking up tighmost cluster and the modified basket.
		b   = Rt$b
		b   = b[which(sapply(b,FUN=length)>0)]
		idx = setdiff(idx, Rt$res) 		#x[! x[,(col+1)] %in% Rt$res[,(col+1)], ]

		index[Rt$res] = t
		Cluster[[t]]  = Rt$res			#get a cluster in the list "Cluster"
		size[[t]]     = length(Cluster[[t]])

        	cat(">>>>Found tight cluster ", t," with ",  length(Rt$res)," data points.\n")
		#cat(">Remaining Bucket Size = ", length(Rt$b),"\n")
	}

	return(list(Cluster = Cluster, size = size, CsLabel = index))	
}
