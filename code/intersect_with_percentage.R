# Intersect with a percentage of common genes in Peaks File
# Author: Nitesh Turaga


# 1. Make a list of lists
# How: lists =  c(list1,list2,list3...)

intersect.threshold = function(lists,threshold) {
		    allids = unique(unlist(lists))
		    idcount = sapply(lists, function(x,y) {y %in% x}, y = allids)
		    idcount_threshold = apply(idcount,1,function(x,thresh){ 
		    		       	sum(x)/length(allids) >thresh},
				       	thresh = threshold)
		    
		    return(idcount_threshold)
		    
		    }