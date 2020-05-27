# R functions used to calculate MAPS, CIs and P-values

calculate_MAPS <- function(data)
{
	expected_singleton=sum(data$Adjusted_contribution)
	observed_singleton=sum(data$gnomAD_AC[data$gnomAD_AC == 1])
	total_variants=nrow(data)
	MAPS=(observed_singleton-expected_singleton)/total_variants
	return(MAPS)
}

permute_MAPS <- function(data)
{
	temp <- vector()
	for (i in 1:nsim)
	{
		sampled_data=sample_n(data, nrow(data), replace = TRUE)
		MAPS_estimate=calculate_MAPS(sampled_data)
		temp[i] <- MAPS_estimate
	}
	sorted_temp<-sort(temp, decreasing=FALSE)
	lower_CI=sorted_temp[nsim/100*5]
	upper_CI=sorted_temp[nsim/100*95]
	return(list(lower_CI,upper_CI))
}

add_estimates <- function(category, class, data)
{
	MPS=calculate_MAPS(data)
	CI_list=permute_MAPS(data)
	lowerCI=CI_list[[1]]
	upperCI=CI_list[[2]]
	to_add=data.frame(variantSet=category, Region=class, MAPS=MPS, lowerCI=lowerCI, upperCI=upperCI, n_variants=nrow(data))
	new_data=rbind(MAPS_data, to_add)
	return(new_data)
}

MAPS_data<-data.frame(variantSet=character(), Region=character(), MAPS=numeric(0), lowerCI=numeric(0), upperCI=numeric(0), n_variants=numeric(0))

# category and class are just descriptors used for plotting
# data is a list of variants observed in gnomAD, including the gnomAD allele count (gnomAD_AC) and the proportion of variants of that type (i.e. given mutational context and methylation status) that are expected to be singletons (Adjusted_contribution) - the coefficients for this are determined using 'MAPSmodel.R'

calculate_P <- function(data1, data2)
{
	n=0
	for (i in 1:nsim)
	{
		sampled_data1=sample_n(data1, nrow(data1), replace = TRUE)
		sampled_data2=sample_n(data2, nrow(data2), replace = TRUE)
		MAPS1=calculate_MAPS(sampled_data1)
		MAPS2=calculate_MAPS(sampled_data2)
		MAPS_diff=MAPS2-MAPS1
		if (MAPS_diff < 0)
		{
			n=n+1
		}
	}
	Pval=n/nsim
 	return(Pval)
}