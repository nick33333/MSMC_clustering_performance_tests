README. Data and code pertaining to Germain, Feng et al. (2023) - Species-specific traits mediate avian demographic responses under past climate change


Species Ne data (folder)
	The following .txt files represent PSMC-based estimates of effective population size (Ne) for 325 bird species
	over the past ~1 million years.

	Each file contains two columns depeciting real time in years (e.g 199.77 represents 199.77 years ago, column 1) 
	and effective population size (column 2). Ne estimates are estimated as x10^4 (e.g. 0.041 represents an effective
	population size of 410).


"Summary_Stats_Ne_Curves_All_Species.csv" contains the names of each species and their corresponding code in the B10k Genomes Project (https://b10k.genomics.cn/) 
as well as basic summary stats on Ne values during the selected time period of 30,000 - 1 million years ago:

	Time_Ne_max = Time point when the species reached their relative maximum effective population size
	Time_Ne_min = Time point when the species reached their relative minimum effective population size
	Range_Ne  = Difference between maximum effective population size and minimum effective population size (Ne_max - Ne_min) 
	Var_Ne = Variance in effective population size
	Mean_Ne = Mean effective population size
	SD_Ne = Standard Deviation in effective population size
	count_Ne = Number of point estimates of Ne used in the previous calculations


////////// Individual data files ////////////

"30_1000.P121.normalized.txt" - Summarized Ne curves (as above) for each species (delineated by B10k ID) from 30kya to 1000kya, stored as a continues string. 
		Used for cluster analysis (see code)

"263.spe.name.txt" - codex file linking B10k ID (column 1) to species name (column 2) and order (column 3)
		Used for cluster analysis (see code)


"B10K_Trait_Data_031319.csv" - Species specific trait data compiled as part of the B10k project. Contains the following columns
		species - species name
		family -  Family
		order -	  Order
		unsexed mass -  mean unsexed mass or average of female and male masses (g)
		b_brain_size -  average of the two sexes or mean unsexed brain volume (mL)
		clutch size - mean number of eggs per clutch
		egg mass - mean fresh egg mass (g)
		maximum longevity - maximum recorded longevity for the species (years)
		mortality both - both sexes, annual mortality rate
		inc duration - incubation duration of the clutch (days)
		min elevation - minimum elevation in which species occurs (m)
		max elevation - maximum elevation in which species occurs (m)
		range size km2 - range size of species based on distributional data for the breeding range (km2)
		migration - 1= sedentary, 2= partially migratory, 3=migratory
		bill total culmen - length from the tip of the bill to the base of the skull (mm)
		bill nares - length from the anterior edge of the nostrils to the tip of the bill (mm)
		bill width - width of the bill at the anterior edge of the nostrils (mm)
		bill depth - depth of the bill at the anterior edge of the nostrils (mm)
		tarsus length - distance from the notch at the knee to the third crease at the ankle (mm)
		kipps distance - distance from the tip of the first secondary to the longest primary (mm)
		wing chord - length from the bend of the wing to the longest primary of the unflattened wing (mm)
		hand wing index - Kipps sidance divided by wing chord length, times 100
		threat status - category description of threat statsus from birdlife international (see code for updates)
		
		
"clusters_different_k.csv" - Identity of Cluster Group each species falls under, using 5 separate values of cluster number (k).  Columns indicate 
		ID - species ID in B10k genome project
		K=3 - Method separating species among 3 cluster groups
		K=4 - Method separating species among 4 cluster groups
		K=5 - Method separating species among 5 cluster groups
		K=6 - Method separating species among 6 cluster groups
		K=7 - Method separating species among 7 cluster groups




"Contemporary data for mix model.csv" - Summary trait data and centroid of breeding range for all bird species
		IUCN species name - species name
		centroid long - centroid of breeding range longitude
		latitude - centroid of breeding range latitude
		abs latitude - absolute value of centroid of breeding range latitude
		region - description of whether or not species central breeding point is in the tropical regions or not
		clutch size - mean number of eggs per clutch
		egg mass (g) -  mean fresh egg mass (g)
		incubation d - incubation duration of the clutch (days)
		body mass - mean unsexed mass or average of female and male masses (g)
		beak length (culmen) - length from the tip of the bill to the base of the skull (mm)
		HWI - Hand Wing Index - Kipps sidance divided by wing chord length, times 100


"Log_Ne_Full_Time_Period_equidistant_time_points.csv" - summaries of full demographic history (log Ne) for each species, separated into 121 equidistance time points
 		ID - species ID in the B10k genome project
		columns V1-V121 - Estimated Log Ne value for each species at that time point


"Normalized_Ne_Full_Time_Period_equidistant_time_points.csv"  - summaries of full demographic history (normalized Ne) for each species, separated into 121 equidistance time points
 		ID - species ID in the B10k genome project
		columns V1-V121 - Estimated Normalized Ne value for each species at that time point



"Pearson_correlation_Summary_Stats_No_Filter_May_2022.csv" - summaries of pearson correlation analysis for each species during both climate warming (IN1) and climate cooling (DC1) periods
		ID - species ID in the B10k genome project
		DC1_coef - correlation coefficient of species demographic response during climate cooling
		DC1_p_value - p value for correlation of species demographic response during climate cooling
		DC1_sample_size	- number of point estimates of Ne used to evaluate relationship between Ne and climate cooling
		DC1_response - overall response (positive/negative/no) of Ne to changing climate during climate cooling
		DC1_Ne_direction - direction (increase/decrease/not signficant) of Ne response to climate cooling		
		IN1_coef - correlation coefficient of species demographic response during climate warming
		IN1_p_value - p value for correlation of species demographic response during climate warming
		IN1_sample_size	- number of point estimates of Ne used to evaluate relationship between Ne and climate warming
		IN1_response - overall response (positive/negative/no) of Ne to changing climate during climate warming
		IN1_Ne_direction - direction (increase/decrease/not signficant) of Ne response to climate warming



"Realms.csv" - Zoogeographic realm for each species, based on Holt, B. G. et al. An Update of Wallace’s Zoogeographic Regions of the World. Science 339, 74–78 (2013)
		B10k ID - species ID in the B10k genome project
		Species Name Standard - Standard scientific name for each species
		Realm_fixed - primary realm, based on breeding location and/or location of sample used in whole genome sequencing
		realms_all - all realms in which the species is found
		





