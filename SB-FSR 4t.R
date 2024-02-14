### A) Proteomics data filtering and protein FSR calculation ###

### Load libraries ###
library(tidyverse)
library(enviPat)
data(isotopes)

### Set working directory ###

# In R Studio, this line of code sets the working directory to the same path as this script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

### Custom functions ###

# This function reads each amino acid in the peptide sequence x and sums the theoretical deuterium labelling
# probabilites for each amino acid in the sequence. The function requires that the table Amino_acids.txt is
# available in the working directory. y defines which column in the Amino_acids.txt table to use for the labelling
# probabilities, hence allowing different values to be used for different species.
deuterium_labels<-function(x,y){
  if(exists("amino_acids")==FALSE){amino_acids<<-read_tsv("Amino_acids.txt")}
  mass_count<-0
  for(i in seq(1,str_length(x),1)){
    mass_count<-mass_count+amino_acids[amino_acids$Code==substring(x,i,i),y]
  }
  print(as.numeric(mass_count))
}

### Import data ###

# Import experimental design information 
sampling_times<-read.csv("Study Design.csv")

# Import FTIR gradients (scale = hours) and intercepts
gradients<-read.csv("gradients.csv")

# List of subject IDs for collection of Skyline data
subject_ID<-unique(sampling_times$Subject)
 skyline_output = NULL
 skip_to_next = FALSE

# # Code to import Skyline data based on list of subjects above
for (i in subject_ID) {
   tryCatch({import <- read.csv(paste(i,".csv",sep="")) 
       import$Subject<-i
       skyline_output<-bind_rows(skyline_output,import)}, error = function(e) { print(paste("Error importing",i)); skip_to_next <<- TRUE})
   if(skip_to_next==TRUE) {next}
 }

rm(list=c("import","i","skip_to_next"))

### First step of QC criteria ###

# Import QC criteria from text file
QC_criteria<-read_tsv("QC_filters.txt")

M1_M2_min<-QC_criteria$Value[QC_criteria$Parameter=="M1_M2_min"]
minM0<-QC_criteria$Value[QC_criteria$Parameter=="minM0"]
minM4<-QC_criteria$Value[QC_criteria$Parameter=="minM4"]

# Find all subject/peptide/charge combinations where M+1 or M+2 or M+3 intensity for the first time point is less than M1_M2_min in QC_criteria file
t0_M1_M2_low_intensity<-skyline_output %>%
  filter(Fragment.Ion %in% c("precursor [M+1]","precursor [M+2]","precursor [M+3]") & 
           (X0.Area<M1_M2_min)) %>%
  select(Peptide.Modified.Sequence.Three.Letter.Codes,Product.Charge,Subject)

# Find all subject/peptide/charge combinations where M+1 OR M+2 OR M+3 intensity for time points after t0 are less than M1_M2_min 
# in QC_criteria file AND where the number of timepoints after t0 above this threshold are fewer than 3
tn_M1_M2_low_intensity<-skyline_output %>%
  mutate(x1.valid=case_when(X1.Area>M1_M2_min~1,TRUE~0),
         x2.valid=case_when(X2.Area>M1_M2_min~1,TRUE~0),
         x3.valid=case_when(X3.Area>M1_M2_min~1,TRUE~0)) %>%
  mutate(total.valid=x1.valid+x2.valid+x3.valid) %>%
  filter(Fragment.Ion %in% c("precursor [M+1]","precursor [M+2]","precursor [M+3]") & 
           total.valid < 3 &
           (X1.Area<M1_M2_min |
              X2.Area<M1_M2_min |
              X3.Area<M1_M2_min)) %>%
  select(Peptide.Modified.Sequence.Three.Letter.Codes,Product.Charge,Subject)


# Export the results of the above two queries to a CSV file

M1_M2_low_intensity<-distinct((rbind(t0_M1_M2_low_intensity,tn_M1_M2_low_intensity))) %>%
  left_join(skyline_output)

write_csv(M1_M2_low_intensity %>%
            mutate(Sequence_Charge=paste(str_replace_all(Peptide.Modified.Sequence.Three.Letter.Codes, "\\[1Ac]", "\\[+42.0105]"),strrep("+",Product.Charge),sep="")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[Oxi]","\\[+15.994915]")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[CAM]","\\[+57.021464]")),"M1_M2_M3_low_intensity.csv")

# Now remove these results from the skyline_output file. Note this code remove ALL time points (see note above
# about additional code needed both here and later in script to remove single timepoints).

skyline_output<-skyline_output %>%
  anti_join(M1_M2_low_intensity)

rm(M1_M2_low_intensity,t0_M1_M2_low_intensity,tn_M1_M2_low_intensity)

# Find all subject/peptide/charge combinations where t0 precursor intensity is less than minM0 in QC_criteria file

t0_low_minM0<-skyline_output %>%
  filter(Fragment.Ion=="precursor" & X0.Area<minM0) %>%
  select(Peptide.Modified.Sequence.Three.Letter.Codes, Subject, Product.Charge) %>%
  left_join(skyline_output)

# Export the results of the above to a CSV file

write_csv(t0_low_minM0 %>%
            mutate(Sequence_Charge=paste(str_replace_all(Peptide.Modified.Sequence.Three.Letter.Codes, "\\[1Ac]", "\\[+42.0105]"),strrep("+",Product.Charge),sep="")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[Oxi]","\\[+15.994915]")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[CAM]","\\[+57.021464]")),"low_t0_precursor.csv")

# Now remove these results from the skyline_output file

skyline_output<-skyline_output %>%
  anti_join(t0_low_minM0)

rm(t0_low_minM0)

# Find all subject/peptide/charge combinations where precursor intensity for time points after t0 are less than minM0 
# in QC_criteria file AND where the number of timepoints after t0 above this threshold are fewer than 3
tn_low_minM0<-skyline_output %>%
  mutate(x1.valid=case_when(X1.Area>minM0~1,TRUE~0),
         x2.valid=case_when(X2.Area>minM0~1,TRUE~0),
         x3.valid=case_when(X3.Area>minM0~1,TRUE~0)) %>%
  mutate(total.valid=x1.valid+x2.valid+x3.valid) %>%
  filter(Fragment.Ion == "precursor" & 
           total.valid < 3 &
           (X1.Area<minM0 |
              X2.Area<minM0 |
              X3.Area<minM0 )) %>%
  select(Peptide.Modified.Sequence.Three.Letter.Codes,Product.Charge,Subject)

# Export the results to a CSV file

tn_low_minM0<-tn_low_minM0 %>%
  left_join(skyline_output)

write_csv(tn_low_minM0 %>%
            mutate(Sequence_Charge=paste(str_replace_all(Peptide.Modified.Sequence.Three.Letter.Codes, "\\[1Ac]", "\\[+42.0105]"),strrep("+",Product.Charge),sep="")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[Oxi]","\\[+15.994915]")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[CAM]","\\[+57.021464]")),"low_tn_precursor.csv")

# Now remove these results from the skyline_output file.

skyline_output<-skyline_output %>%
  anti_join(tn_low_minM0)

rm(tn_low_minM0)

# Select all precursor M+4 intensities less than minM4 in QC_criteria file

low_minM4<-skyline_output %>%
  filter(Fragment.Ion =="precursor [M+4]" &
           (X0.Area<minM4 |
             X1.Area<minM4 |
              X2.Area<minM4 |
              X3.Area<minM4 ))

# Export the results to a CSV file

write_csv(low_minM4 %>%
            mutate(Sequence_Charge=paste(str_replace_all(Peptide.Modified.Sequence.Three.Letter.Codes, "\\[1Ac]", "\\[+42.0105]"),strrep("+",Product.Charge),sep="")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[Oxi]","\\[+15.994915]")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[CAM]","\\[+57.021464]")),"low_minM4.csv")

# Now change all M+4 intensitives to zero for the whole time series where one or more precursor M+4 intensities in
# the time series are less than minM4 in QC_criteria file

low_minM4$X0.Area<-0
low_minM4$X1.Area<-0
low_minM4$X2.Area<-0
low_minM4$X3.Area<-0

skyline_output<-skyline_output %>%
  anti_join(low_minM4 %>% select(Peptide.Modified.Sequence.Three.Letter.Codes,Product.Charge,Subject,Fragment.Ion))

skyline_output<-rbind(skyline_output,low_minM4)

rm(low_minM4)

# Now apply the retention time filters
within_run_RT<-QC_criteria$Value[QC_criteria$Parameter=="within_run_RT"]
between_run_RT<-QC_criteria$Value[QC_criteria$Parameter=="between_run_RT"]

RT_filter<-skyline_output %>%
  select(21,6,9,12,14,16,18,20) %>%
  gather(day_long,RT,5:8) %>% 
  mutate(Peptide.Sequence=Peptide.Modified.Sequence.Three.Letter.Codes) %>%
  arrange(Subject,Peptide.Sequence,Product.Charge) %>%
  group_by(Subject,Peptide.Sequence,Product.Charge,day_long) %>%
  summarise(within_RT_CV=sd(RT)/mean(RT),mean_RT=mean(RT)) %>%
  group_by(Subject,Peptide.Sequence,Product.Charge) %>%
  summarise(mean_within_run_RT_CV=mean(within_RT_CV),between_run_meanRT_CV=sd(mean_RT)/mean(mean_RT)) %>%
  filter(mean_within_run_RT_CV>within_run_RT | between_run_meanRT_CV>between_run_RT)

# Export the results to a CSV file

write_csv(RT_filter %>%
            mutate(Sequence_Charge=paste(str_replace_all(Peptide.Sequence, "\\[1Ac]", "\\[+42.0105]"),strrep("+",Product.Charge),sep="")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[Oxi]","\\[+15.994915]")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[CAM]","\\[+57.021464]")),"RT_filter.csv")

# Now remove these results from the skyline_output file. Removes whole time series.

skyline_output<-skyline_output %>%
  anti_join(RT_filter)

rm(RT_filter)

### Perform peptide calculations ###

# Extract unique peptide sequences, calculate theoretical deuterium labelling and then join back on to Skyline output
# At the moment, this just uses albumin values for chickens, however can customise here to use other values
peptide_calculations_list<-sapply(unique(skyline_output$Peptide.Sequence),deuterium_labels,y="Exch_deut_chicken_albumin")
peptide_calculations<-as_tibble(cbind(Sequence=as.character(unique(skyline_output$Peptide.Sequence)),Deut_labs=peptide_calculations_list))
skyline_output<-skyline_output %>%
  full_join(peptide_calculations,by=c("Peptide.Sequence"="Sequence")) %>%
  mutate(Deut_labs=as.numeric(Deut_labs))

### Now calculate the theoretical abundances using the EnviPat package

# Import amino acid formulae. Note these are condensed amino acids with C2H3NO modification to all cysteines.
aa_formulae<-read_tsv("aa_formulae.txt")

# Generage an enviPat compatible formula for all peptides

aa_sequence<-peptide_calculations$Sequence

# Turn sequence string into list of characters

formula_list=NULL

for(j in aa_sequence){

  split_sequence<-NULL

  for(i in seq(1,str_length(j))){
    split_sequence<-c(split_sequence,substr(j,i,i))
  }

  # Join sequence on to amino acid table
  formula_list<-formula_list %>%
    bind_rows(as_tibble(split_sequence) %>%
                rename(Single_letter=value) %>%
                mutate(peptide_sequence=j) %>%
                left_join(aa_formulae))

}

rm(split_sequence,i,j,aa_sequence,amino_acids)

# Summarise atoms, add H30
formula_list<-formula_list %>%
  group_by(peptide_sequence) %>%
  summarise(C=sum(C),H=sum(H)+3,N=sum(N),O=sum(O)+1,S=sum(S))

# Now join back on to skyline_output to find all the peptides with 1Ac, Oxi and CAM modifications
# Remember that CAM modification (C2H3N0 has already been added to the cysteine formula)
# Also count the 1Ac and Oxi modifications in each peptide, add as new columns and then add
# C2H2O for each 1Ac modification and O for each Oxi modification and create enviPat compatible formula

formula_list<-skyline_output %>%
  select(Peptide.Modified.Sequence.Three.Letter.Codes,peptide_sequence=Peptide.Sequence) %>%
  distinct() %>%
  left_join(formula_list) %>%
  rowwise() %>%
  mutate(Oxi=str_count(Peptide.Modified.Sequence.Three.Letter.Codes,"\\[Oxi]"),
         Ac1=str_count(Peptide.Modified.Sequence.Three.Letter.Codes,"\\[1Ac]"))%>%
  mutate(O=O+Oxi+Ac1,C=C+2*Ac1,H=H+2*Ac1) %>%
  mutate(pattern=paste("C",C,"H",H,"N",N,"O",O,"S",S,sep=""))

# Find all formulae with no sulphurs in them and remove the S0 at the end of the formula
S_is_zero<-formula_list %>%
  filter(S==0) %>%
  mutate(pattern=gsub('.{2}$', '', pattern))

# Replace formulae ending in S0 in the original table with the formula without the S0
formula_list<-formula_list %>%
  filter(S>0) %>%
  bind_rows(S_is_zero)

rm(S_is_zero)

# Perform enviPat calculations
pattern<-isopattern(
  isotopes,
  unique(formula_list$pattern),
  threshold=0.1,
  plotit=FALSE,
  charge=FALSE,
  emass=0.00054858,
  algo=1
)

formulae<-names(pattern)
enviPat_output<-NULL
counter<-1

# Extract enviPat abundance calculations with m/z at integer level
for(i in pattern){
  enviPat_output<-enviPat_output %>%
    bind_rows(as_tibble(i) %>%
                mutate(pattern=formulae[counter])%>%
                select(pattern,mz='m/z',abundance) %>%
                mutate(mz=round(mz,0)) %>%
                group_by(pattern,mz) %>%
                summarise(abundance=sum(abundance)))
  counter=counter+1
}

  enviPat_output<-enviPat_output %>%
  left_join(enviPat_output %>%
              group_by(pattern) %>%
              summarise(minimum_mz=min(mz))) %>%
  filter(mz<=minimum_mz+5)

enviPat_output<-enviPat_output %>%
  left_join(enviPat_output %>%
              group_by(pattern) %>%
              summarise(sum_abundance=sum(abundance))) %>%
  mutate(relative_abundance=abundance/sum_abundance) %>%
  left_join(formula_list %>%
              select(Peptide.sequence=Peptide.Modified.Sequence.Three.Letter.Codes,pattern))

rm(i,formulae,counter,pattern,peptide_calculations_list,peptide_calculations,formula_list,aa_formulae)

# This bit of code removes all the Fragment.Ion fields labelled precursor -64 and reduces the Transition.Count by 1
precursor64<-skyline_output %>%
  right_join(skyline_output %>%
  filter(Fragment.Ion=="precursor -64") %>%
  select(Subject,Peptide.Modified.Sequence.Three.Letter.Codes,Product.Charge)) %>%
  mutate(Transition.Count=Transition.Count-1) %>%
  filter(Fragment.Ion!="precursor -64")
         
skyline_output <- skyline_output %>%
  anti_join(skyline_output %>%
              filter(Fragment.Ion=="precursor -64") %>%
              select(Subject,Peptide.Modified.Sequence.Three.Letter.Codes,Product.Charge))   %>%
  bind_rows(precursor64)

rm(precursor64)

# Now remove all rows of data where deut_labs is less than deut_labs ins the QC_criteria file.
deut_labs_QC<-QC_criteria$Value[QC_criteria$Parameter=="deut_labs"]

deut_labs_filter<-skyline_output %>%
  filter(Deut_labs<deut_labs_QC)

skyline_output<-skyline_output %>%
  anti_join(deut_labs_filter)

# Export the results to a CSV file

write_csv(deut_labs_filter%>%
            mutate(Sequence_Charge=paste(str_replace_all(Peptide.Modified.Sequence.Three.Letter.Codes, "\\[1Ac]", "\\[+42.0105]"),strrep("+",Product.Charge),sep="")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[Oxi]","\\[+15.994915]")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[CAM]","\\[+57.021464]")),"deut_labs_filter.csv")

rm(deut_labs_filter)

### Prepare Skyline output for FSR calculations ###
# Take wide Skyline data and convert to long format. Beware, this is dangerous code as it
# uses column numbers, so potential for errors if Skyline output changes. Currently
# configured for 4 timepoints over one phase and M+3
phase1<-skyline_output %>%
  select(Subject,Protein.Gene,Protein.Accession,Protein.Name,Protein.Description,Peptide.Modified.Sequence.Three.Letter.Codes,Deut_labs,Fragment.Ion,Product.Charge,Transition.Count,X0.Area,X1.Area,X2.Area,X3.Area) %>%
  gather(day_long,Area,11:14) %>%
  spread('Fragment.Ion','Area')

# This code inserts a zero into all the empty precursor[M+4] fields, but not used at present, as limited to [M+3]
# phase1$`precursor [M+4]`[is.na(phase1$`precursor [M+4]`)]<-0

# This code calculates the mole_perc
phase1<-phase1 %>%
  mutate(precursor=as.numeric(precursor),`precursor [M+1]`=as.numeric(`precursor [M+1]`),`precursor [M+2]`=as.numeric(`precursor [M+2]`),`precursor [M+3]`=as.numeric(`precursor [M+3]`)) %>%
  mutate(day_seq=as.ordered(as.numeric(str_extract(day_long,regex("\\d+")))),
    sumM=rowSums(.[11:14],na.rm=TRUE),
    mole_perc=100*(.[[12]]/rowSums(.[11:14],na.rm=TRUE)+2*.[[13]]/rowSums(.[11:14],na.rm=TRUE)+3*.[[14]]/rowSums(.[11:14],na.rm=TRUE))) %>%
  rename(Peptide.Sequence=Peptide.Modified.Sequence.Three.Letter.Codes)  %>%
  arrange(Subject,Peptide.Sequence,Product.Charge,day_seq)

# Calculate the mean Transition.Count for each Subject/Peptide/Charge combination and join back onto phase1
phase1<-phase1 %>%
  left_join(phase1 %>%
  group_by(Subject,Peptide.Sequence,Product.Charge) %>%
  summarise(Mean_Trans=mean(Transition.Count)))

# Check to see whether there are any sample/peptide/charge sets within a time series that have different transition counts
phase1 %>%
  filter(Mean_Trans != Transition.Count)

### Also worth at this stage just checking the phase1 tibble to ensure that all peptides have a transition
### count of at least 4 and there are values in the M+1, M+2 and M+3 cells.

# This code calculates the mpe
mole_perc_baseline_1<-phase1 %>%
  filter(day_long=="X0.Area") %>%
  select(Subject,Peptide.Sequence,Product.Charge,baseline_mole_perc=mole_perc)

phase1<-phase1 %>%
  full_join(mole_perc_baseline_1) %>%
  mutate(mpe=mole_perc-baseline_mole_perc) %>%
  left_join(sampling_times%>%
              select(Subject,day_seq=Sequence,calc_interval=Interval..h.)%>%
              mutate(day_seq=as.ordered(day_seq)))

rm(mole_perc_baseline_1)

# Now calculate theoretical baseline up to M+3 and MPE
theoretical_baseline<-enviPat_output %>%
  filter(mz<=minimum_mz+3) %>%
  arrange(Peptide.sequence,mz)

theoretical_baseline$mass<-rep(c('M0','M+1','M+2','M+3'),length(unique(enviPat_output$Peptide.sequence)))

phase1<-phase1 %>%
  left_join(theoretical_baseline %>%
  pivot_wider(id_cols=Peptide.sequence, names_from=mass, values_from=abundance) %>%
  mutate(sumM=rowSums(.[2:5],na.rm=TRUE),
         theory_baseline=100*(.[[3]]/rowSums(.[2:5],na.rm=TRUE)+2*.[[4]]/rowSums(.[2:5],na.rm=TRUE)+3*.[[5]]/rowSums(.[2:5],na.rm=TRUE))) %>%
  select(Peptide.Sequence=Peptide.sequence,theory_baseline)) %>%
  mutate(baseline_actual_over_theory=baseline_mole_perc/theory_baseline,mpe_theory=mole_perc-theory_baseline)

# Set the theoretical MPE on day 0 to 0

phase1$mpe_theory[phase1$day_seq==0]<-0

### Perform FSR calculation

# Performed by finding the minimum RMS (root mean of the squares) of the difference
# between the observed vs calculated mpe for any value of k1. Starts looking between 0 and 0.5 (i.e. 50%/h)

FSR<-phase1 %>%
#  filter(Peptide.Sequence %in% c("DLQHRLDEAEQLALK","DPLNETVIGLYQK")) %>%
  group_by(Subject,Protein.Gene, Protein.Accession,Protein.Description,Peptide.Sequence,Product.Charge,baseline_mole_perc,theory_baseline,baseline_actual_over_theory) %>%
  select(Subject,Protein.Gene, Peptide.Sequence,Product.Charge,Deut_labs,baseline_mole_perc,theory_baseline,baseline_actual_over_theory) %>%
  #group_by(Subject,Protein.Gene, Protein.Accession,Protein.Description,Peptide.Sequence,Product.Charge) %>%
  #select(Subject,Protein.Gene, Peptide.Sequence,Product.Charge,Deut_labs) %>%
    summarise(Deut_labs=mean(Deut_labs)) %>%
  left_join(gradients) %>%
  rowwise() %>%
  mutate(fsr_temp = paste(optim(c(0,0.5),function(k1,c=0,d=0,r=1,k2=Gradient,n=Deut_labs,p=exp(Intercept),
                                                ph=Product.Charge,subj=Subject,pep_seq=Peptide.Sequence){
    print(pep_seq)
    as.numeric(phase1 %>%
                 filter(Product.Charge==ph,Subject==subj,Peptide.Sequence==pep_seq) %>%
                 select(Subject,Product.Charge,Peptide.Sequence,mpe,calc_interval) %>%
                 mutate(Fit=r*n*p*exp(k2*(calc_interval))*(1-exp(-k1*(calc_interval-d)))) %>%
                 mutate(Diff=(Fit-mpe)^2) %>%
                 summarise(RMS=sqrt(mean(Diff))) %>%
                 select(RMS))
  }),collapse=",")) %>%
  mutate(fsr_min_hr=100*as.numeric(substr(fsr_temp,3,regexpr(",",fsr_temp)-1)),
         fsr_max_hr=100*as.numeric(substr(fsr_temp,regexpr(",",fsr_temp)+2,regexpr(")",fsr_temp)-1)),
         Fit=as.numeric(substr(fsr_temp,regexpr(")",fsr_temp)+2,regexpr(")",fsr_temp)+7))) %>%
  select(-fsr_temp) %>%
  mutate(fsr_hr=(fsr_max_hr+fsr_min_hr)/2,fsr_day=24*(fsr_max_hr+fsr_min_hr)/2) 
  # %>%
  # mutate(fsr_temp_theory = paste(optim(c(0,0.5),function(k1,c=0,d=0,r=1,k2=Gradient,n=Deut_labs,p=exp(Intercept),
  #                                                 ph=Product.Charge,subj=Subject,pep_seq=Peptide.Sequence){
  #   print(pep_seq)
  #   as.numeric(phase1 %>%
  #                filter(Product.Charge==ph,Subject==subj,Peptide.Sequence==pep_seq) %>%
  #                select(Subject,Product.Charge,Peptide.Sequence,mpe_theory,calc_interval) %>%
  #                mutate(Fit_theory=r*n*p*exp(k2*(calc_interval))*(1-exp(-k1*(calc_interval-d)))) %>%
  #                mutate(Diff_theory=(Fit_theory-mpe_theory)^2) %>%
  #                summarise(RMS_theory=sqrt(mean(Diff_theory))) %>%
  #                select(RMS_theory))
  # }),collapse=",")) %>%
  # mutate(fsr_min_hr_theory=100*as.numeric(substr(fsr_temp_theory,3,regexpr(",",fsr_temp_theory)-1)),
  #        fsr_max_hr_theory=100*as.numeric(substr(fsr_temp_theory,regexpr(",",fsr_temp_theory)+2,regexpr(")",fsr_temp_theory)-1)),
  #        Fit_theory=as.numeric(substr(fsr_temp_theory,regexpr(")",fsr_temp_theory)+2,regexpr(")",fsr_temp_theory)+7))) %>%
  # select(-fsr_temp_theory) %>%
  # mutate(fsr_hr_theory=(fsr_max_hr_theory+fsr_min_hr_theory)/2,fsr_day_theory=24*(fsr_max_hr_theory+fsr_min_hr_theory)/2)

# Remove high FSR calculations or those with a big difference between min and max

QC_High_FSR<-QC_criteria$Value[QC_criteria$Parameter=="FSR_max"]
QC_FSR_range<-QC_criteria$Value[QC_criteria$Parameter=="FSR_diff"]

FSR_high<-FSR %>%
  filter(fsr_hr>QC_High_FSR | (fsr_min_hr-fsr_max_hr)^2>QC_FSR_range^2)

FSR<-FSR %>%
  anti_join(FSR_high)

# Export these high FSRs and those with a big difference between min and max

write_csv(FSR_high %>%
            mutate(Sequence_Charge=paste(str_replace_all(Peptide.Sequence, "\\[1Ac]", "\\[+42.0105]"),strrep("+",Product.Charge),sep="")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[Oxi]","\\[+15.994915]")) %>%
            mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[CAM]","\\[+57.021464]")),"FSR_high.csv")

# Now recalculate these very high FSRs using a tighter starting range of 0 to 10% /hr

FSR_high<-FSR_high %>%
  select(Subject, Peptide.Sequence, Product.Charge)

phase_recalc<-phase1 %>%
  inner_join(FSR_high)

FSR_recalc<-phase_recalc %>%
  group_by(Subject,Protein.Gene, Protein.Accession,Protein.Description,Peptide.Sequence,Product.Charge,baseline_mole_perc,theory_baseline,baseline_actual_over_theory) %>%
  select(Subject,Protein.Gene, Peptide.Sequence,Product.Charge,Deut_labs,baseline_mole_perc,theory_baseline,baseline_actual_over_theory) %>%
  #group_by(Subject,Protein.Gene, Protein.Accession,Protein.Description,Peptide.Sequence,Product.Charge) %>%
  #select(Subject,Protein.Gene, Peptide.Sequence,Product.Charge,Deut_labs) %>%
  summarise(Deut_labs=mean(Deut_labs)) %>%
  left_join(gradients) %>%
  rowwise() %>%
  mutate(fsr_temp = paste(optim(c(0,0.1),function(k1,c=0,d=0,r=1,k2=Gradient,n=Deut_labs,p=exp(Intercept),
                                                  ph=Product.Charge,subj=Subject,pep_seq=Peptide.Sequence){
    print(pep_seq)
    as.numeric(phase_recalc %>%
                 filter(Product.Charge==ph,Subject==subj,Peptide.Sequence==pep_seq) %>%
                 select(Subject,Product.Charge,Peptide.Sequence,mpe,calc_interval) %>%
                 mutate(Fit=r*n*p*exp(k2*(calc_interval))*(1-exp(-k1*(calc_interval-d)))) %>%
                 mutate(Diff=(Fit-mpe)^2) %>%
                 summarise(RMS=sqrt(mean(Diff))) %>%
                 select(RMS))
  }),collapse=",")) %>%
  mutate(fsr_min_hr=100*as.numeric(substr(fsr_temp,3,regexpr(",",fsr_temp)-1)),
         fsr_max_hr=100*as.numeric(substr(fsr_temp,regexpr(",",fsr_temp)+2,regexpr(")",fsr_temp)-1)),
         Fit=as.numeric(substr(fsr_temp,regexpr(")",fsr_temp)+2,regexpr(")",fsr_temp)+7))) %>%
  select(-fsr_temp) %>%
  mutate(fsr_hr=(fsr_max_hr+fsr_min_hr)/2,fsr_day=24*(fsr_max_hr+fsr_min_hr)/2)
  # %>%
  # mutate(fsr_temp_theory = paste(optim(c(0,0.1),function(k1,c=0,d=0,r=1,k2=Gradient,n=Deut_labs,p=exp(Intercept),
  #                                                        ph=Product.Charge,subj=Subject,pep_seq=Peptide.Sequence){
  #   print(pep_seq)
  #   as.numeric(phase1 %>%
  #                filter(Product.Charge==ph,Subject==subj,Peptide.Sequence==pep_seq) %>%
  #                select(Subject,Product.Charge,Peptide.Sequence,mpe_theory,calc_interval) %>%
  #                mutate(Fit_theory=r*n*p*exp(k2*(calc_interval))*(1-exp(-k1*(calc_interval-d)))) %>%
  #                mutate(Diff_theory=(Fit_theory-mpe_theory)^2) %>%
  #                summarise(RMS_theory=sqrt(mean(Diff_theory))) %>%
  #                select(RMS_theory))
  # }),collapse=",")) %>%
  # mutate(fsr_min_hr_theory=100*as.numeric(substr(fsr_temp_theory,3,regexpr(",",fsr_temp_theory)-1)),
  #        fsr_max_hr_theory=100*as.numeric(substr(fsr_temp_theory,regexpr(",",fsr_temp_theory)+2,regexpr(")",fsr_temp_theory)-1)),
  #        Fit_theory=as.numeric(substr(fsr_temp_theory,regexpr(")",fsr_temp_theory)+2,regexpr(")",fsr_temp_theory)+7))) %>%
  # select(-fsr_temp_theory) %>%
  # mutate(fsr_hr_theory=(fsr_max_hr_theory+fsr_min_hr_theory)/2,fsr_day_theory=24*(fsr_max_hr_theory+fsr_min_hr_theory)/2)

# Join these recalculated FSRs back on to the main FSR table

FSR<-rbind(FSR,FSR_recalc)

### Change FIT parameter from RMSE to RRMS

FSR<-FSR %>%
  left_join(phase1 %>%
              select(Subject,Peptide.Sequence, Product.Charge, mpe) %>%
              group_by(Subject,Peptide.Sequence, Product.Charge) %>%
              summarise(ave_mpe=mean(mpe))) %>%
  mutate(Fit=Fit/ave_mpe*100) %>%
  select(-ave_mpe)

### Generate QC parameters

colnames(phase1)[12:15]<-c("precursorM1","precursorM2","precursorM3","precursorM4")

FSR<-inner_join(FSR, phase1 %>%
                  select(Subject,Peptide.Sequence,Product.Charge,precursor,precursorM1,precursorM2,precursorM3,precursorM4,sumM) %>%
                  group_by(Subject,Peptide.Sequence,Product.Charge) %>%
                  summarise(minM0=min(precursor),minM1=min(precursorM1),minM2=min(precursorM2),minM3=min(precursorM3),minM4=min(precursorM4),minSumM=min(sumM)))

phase1_QC<-skyline_output %>%
  select(21,6,9,12,14,16,18,20) %>%
  gather(day_long,RT,5:8) %>% 
  mutate(Peptide.Sequence=Peptide.Modified.Sequence.Three.Letter.Codes) %>%
  arrange(Subject,Peptide.Sequence,Product.Charge) %>%
  group_by(Subject,Peptide.Sequence,Product.Charge,day_long) %>%
  summarise(within_RT_CV=sd(RT)/mean(RT),mean_RT=mean(RT)) %>%
  group_by(Subject,Peptide.Sequence,Product.Charge) %>%
  summarise(mean_within_run_RT_CV=mean(within_RT_CV),between_run_meanRT_CV=sd(mean_RT)/mean(mean_RT))

FSR<-FSR %>%
  inner_join(phase1_QC)

FSR<-FSR %>%
  mutate(Sequence_Charge=paste(str_replace_all(Peptide.Sequence, "\\[1Ac]", "\\[+42.0105]"),strrep("+",Product.Charge),sep="")) %>%
  mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[Oxi]","\\[+15.994915]")) %>%
  mutate(Sequence_Charge=str_replace_all(Sequence_Charge,"\\[CAM]","\\[+57.021464]"))

# Remove the low FSR calculations and any recalculated FSRs that are still high or have a large min/max range

QC_low_FSR<-QC_criteria$Value[QC_criteria$Parameter=="FSR_low"]

FSR_low_high<-FSR %>%
  filter(fsr_hr<QC_low_FSR | fsr_hr>QC_High_FSR | (fsr_min_hr-fsr_max_hr)^2>QC_FSR_range^2)

write_csv(FSR_low_high,"FSR_low_high.csv")

FSR<-FSR %>%
  anti_join(FSR_low_high)

# Now calcualte FSRs at the protein level and generate summary statistics 

FSR_summary<-FSR %>%
  group_by(Protein.Description,Protein.Accession,Protein.Gene,Subject) %>%
  summarise(mean_FSR_day=mean(fsr_day),SD_FSR_day=sd(fsr_day),median_FSR_day=median(fsr_day),IQR_FSR_day=IQR(fsr_day),peptides=n(),
            unique_peptides=length(unique(Peptide.Sequence)))

FSR_Z_score<-FSR %>%
  left_join(FSR_summary) %>%
  mutate(Z_score=(fsr_day-mean_FSR_day)/SD_FSR_day)

# Set Z score to zero where there are 5 or fewer peptides used to calculate the mean and SD from which the Z score is calculated

FSR_Z_score$Z_score[FSR_Z_score$peptides<=5]<-0

# Subject summary comprising: number of proteins identified per subject, number of proteins with more than 5 peptides (per subject), number of peptides (per subject), number of Unique peptides (per subject)

subject_summary<-FSR_summary %>%
  mutate(more_than_5_peptides=case_when(peptides>5 ~ 1,
         peptides<=5 ~ 0)) %>%
  group_by(Subject) %>%
  summarise(No_proteins=n(), No_proteins_over_5_peptides=sum(more_than_5_peptides), No_peptides=sum(peptides),No_unique_peptides=sum(unique_peptides))

# Export the peptide FSR data, protein FSR_summary and Z-scores data

write_csv(FSR,"FSR.csv")
write_csv(FSR_summary,"FSR_summary.csv")
write_csv(FSR_Z_score,"FSR_Z_score.csv")
write_csv(subject_summary,"subject_summary.csv")
write_csv(phase1,"phase1.csv")


### B) Post FSR calculation data filtering ###

### Load libraries ###
library(tidyverse)
library(lme4)
library(afex)

### Set working directory ###

# In R Studio, this line of code sets the working directory to the same path as this script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Import FSR.csv and remove FSR calculations where the Fit parameter is more than that specified in the QC_criteria file 
# (QC_criteria file is reloaded to allow the FSR_fit parameter to be adjusted)

FSR<-read.csv("FSR.csv")

QC_criteria<-read_tsv("QC_filters.txt")
QC_poor_fit<-QC_criteria$Value[QC_criteria$Parameter=="FSR_fit"]

FSR_poor_fit<-FSR %>%
  filter(Fit>QC_poor_fit | Fit<0)

write_csv(FSR_poor_fit,"FSR_poor_fit.csv")

FSR<-FSR %>%
  anti_join(FSR_poor_fit)

rm(FSR_poor_fit)

# Also remove all FSRs where baseline_actual_over_theory is outwith the the tolerance set in the QC_criteria file.
baseline_QC_top<-QC_criteria$Value[QC_criteria$Parameter=="baseline_vs_theory_top"]
baseline_QC_bottom<-QC_criteria$Value[QC_criteria$Parameter=="baseline_vs_theory_bottom"]

baseline_filter<-FSR %>%
  filter(!between(baseline_actual_over_theory,1-baseline_QC_bottom,1+baseline_QC_top))

FSR<-FSR %>%
  anti_join(baseline_filter)

# Export the results to a CSV file

write_csv(baseline_filter,"basline_filter.csv")

rm(baseline_filter)

# Now calcualte FSRs at the protein level and generate summary statistics using these "cleaned up" (final) FSRs

FSR_summary<-FSR %>%
  group_by(Protein.Description,Protein.Accession,Protein.Gene,Subject) %>%
  summarise(mean_FSR_day=mean(fsr_day),SD_FSR_day=sd(fsr_day),median_FSR_day=median(fsr_day),IQR_FSR_day=IQR(fsr_day),peptides=n(),
            unique_peptides=length(unique(Peptide.Sequence)))

FSR_Z_score<-FSR %>%
  left_join(FSR_summary) %>%
  mutate(Z_score=(fsr_day-mean_FSR_day)/SD_FSR_day)

# Set Z score to zero where there are 5 or fewer peptides used to calculate the mean and SD from which the Z score is calculated

FSR_Z_score$Z_score[FSR_Z_score$peptides<=5]<-0

# Subject summary comprising: number of proteins identified per subject, number of proteins with more than 5 peptides (per subject), number of peptides (per subject), number of Unique peptides (per subject)

subject_summary<-FSR_summary %>%
  mutate(more_than_5_peptides=case_when(peptides>5 ~ 1,
                                        peptides<=5 ~ 0)) %>%
  group_by(Subject) %>%
  summarise(No_proteins=n(), No_proteins_over_5_peptides=sum(more_than_5_peptides), No_peptides=sum(peptides),No_unique_peptides=sum(unique_peptides))

# Export the "final" peptide FSR data, protein FSR_summary and Z-scores data

write_csv(FSR,"FSR_final.csv")
write_csv(FSR_summary,"FSR_summary_final.csv")
write_csv(FSR_Z_score,"FSR_Z_score_final.csv")
write_csv(subject_summary,"subject_summary_final.csv")

### C) MPE and FSR plotting ###

### Load libraries ###
library(tidyverse)

### Set working directory ###

# In R Studio, this line of code sets the working directory to the same path as this script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

### Plot mpe vs. interval for all peptides

FSR<-read.csv("FSR_final.csv")
phase1<-read.csv("phase1.csv")
FSR_low_high<-read.csv("FSR_low_high.csv")

pdf("mpe_plots.pdf",16,8)

for(i in unique(na.omit(FSR$Peptide.Sequence))){
  
  protein<-unique(phase1 %>%
                    filter(Peptide.Sequence==i) %>%
                    select(Protein.Description))
  
  curve_values<-FSR %>%
    filter(Peptide.Sequence==i) %>% 
    slice(rep(1:n(), each = 100))
  
  curve_values<-cbind(curve_values,interval=rep(seq(0,as.numeric(phase1 %>%
                                                                   filter(Peptide.Sequence==i) %>%
                                                                   summarise(time=max(calc_interval))),length.out=100),nrow(FSR %>%
                                                                                                                              filter(Peptide.Sequence==i)))) %>%
    mutate(Fit=Deut_labs*exp(Intercept)*exp(Gradient*(interval))*(1-exp(-(fsr_hr/100)*(interval))))
  
  p<-ggplot(phase1 %>%
              filter(Peptide.Sequence==i)) +
    geom_point(aes(x=calc_interval,y=mpe)) +
    geom_line(data=curve_values,aes(x=interval,y=Fit)) +
    xlab("Time (hours)") +
    ylab("mpe (%)") +
    ggtitle(paste(protein,i)) +
    facet_grid(Product.Charge~Subject)
  
  print(p)
  
}

dev.off()

pdf("mpe_plots_low_high_FSR.pdf",16,8)

for(i in unique(na.omit(FSR_low_high$Peptide.Sequence))){
  
  protein<-unique(phase1 %>%
                    filter(Peptide.Sequence==i) %>%
                    select(Protein.Description))
  
  curve_values<-FSR_low_high %>%
    filter(Peptide.Sequence==i) %>% 
    slice(rep(1:n(), each = 100))
  
  curve_values<-cbind(curve_values,interval=rep(seq(0,as.numeric(phase1 %>%
                                                                   filter(Peptide.Sequence==i) %>%
                                                                   summarise(time=max(calc_interval))),length.out=100),nrow(FSR_low_high %>%
                                                                                                                              filter(Peptide.Sequence==i)))) %>%
    mutate(Fit=Deut_labs*exp(Intercept)*exp(Gradient*(interval))*(1-exp(-(fsr_hr/100)*(interval))))
  
  p<-ggplot(phase1 %>%
              filter(Peptide.Sequence==i)) +
    geom_point(aes(x=calc_interval,y=mpe)) +
    geom_line(data=curve_values,aes(x=interval,y=Fit)) +
    xlab("Time (hours)") +
    ylab("mpe (%)") +
    ggtitle(paste(protein,i)) +
    facet_grid(Product.Charge~Subject)
  
  print(p)
  
}

dev.off()

### Plot the FSRs
library(gridExtra)

FSR<-FSR %>%
  mutate(Intervention=str_sub(str_remove_all(Subject,"[0-9]+"),
                              regexpr(" ",str_remove_all(Subject,"[0-9]+")),
                              str_length(str_remove_all(Subject,"[0-9]+"))),
         Experiment=str_sub(Subject,0,regexpr(" ",Subject)))

pdf("FSR_plots_colour_by_peptide.pdf",12,8)

for(i in unique(FSR$Protein.Description)){
  
  p2<-ggplot(FSR %>% filter(Protein.Description==i)) +
    geom_jitter(aes(x=Intervention,y=fsr_day,colour=Peptide.Sequence),width=0.2) +
    ggtitle(i) +
    facet_wrap(~Experiment) +
    theme(legend.position = "none") +
    xlab("")
  
  p3<-ggplot(FSR %>% filter(Protein.Description==i)) +
    geom_boxplot(aes(x=Intervention,y=fsr_day,colour=Peptide.Sequence)) +
    facet_wrap(~Experiment) +
    theme(legend.position = "none") +
    xlab("")
  
  grid.arrange(p2,p3,ncol=1)
  
}

dev.off()

pdf("FSR_plots_colour_by_peptide_replot.pdf",12,8)

for(i in unique(FSR$Protein.Description)){
  
  p1<-ggplot(FSR %>% filter(Protein.Description==i & Experiment == "Experiment1 ")) +
    geom_jitter(aes(x=Intervention,y=fsr_day,colour=Peptide.Sequence),width=0.2) +
    ggtitle(i) +
    theme(legend.position = "none") +
    xlab("")
  
  p2<-ggplot(FSR %>% filter(Protein.Description==i & Experiment == "Experiment2 ")) +
    geom_jitter(aes(x=Intervention,y=fsr_day,colour=Peptide.Sequence),width=0.2) +
    theme(legend.position = "none") +
    xlab("")
  
  p3<-ggplot(FSR %>% filter(Protein.Description==i & Experiment == "Experiment1 ")) +
    geom_boxplot(aes(x=Intervention,y=fsr_day,colour=Peptide.Sequence)) +
    theme(legend.position = "none") +
    xlab("")
  
  p4<-ggplot(FSR %>% filter(Protein.Description==i & Experiment == "Experiment2 ")) +
    geom_boxplot(aes(x=Intervention,y=fsr_day,colour=Peptide.Sequence)) +
    theme(legend.position = "none") +
    xlab("")
  
  grid.arrange(p1,p2,p3,p4,ncol=2)
  
}

dev.off()




pdf("FSR_plots_colour_by_subject.pdf",12,8)

for(i in unique(FSR$Protein.Description)){
  
  p5<-ggplot(FSR %>% filter(Protein.Description==i)) +
    geom_jitter(aes(x=Intervention,y=fsr_day,colour=Subject),width=0.2) +
    ggtitle(i) +
    facet_wrap(~Experiment) +
    theme(legend.position = "none") +
    xlab("")
  
  p6<-ggplot(FSR %>% filter(Protein.Description==i)) +
    geom_boxplot(aes(x=Intervention,y=fsr_day,colour=Subject)) +
    facet_wrap(~Experiment) +
    theme(legend.position = "none") +
    xlab("")
  
  grid.arrange(p5,p6,ncol=1)
  
  
}

dev.off()



pdf("FSR_plots_colour_by_subject_replot.pdf",12,8)

for(i in unique(FSR$Protein.Description)){
  
  p7<-ggplot(FSR %>% filter(Protein.Description==i & Experiment == "Experiment1 ")) +
    geom_jitter(aes(x=Intervention,y=fsr_day,colour=Subject),width=0.2) +
    ggtitle(i) +
    theme(legend.position = "none") +
    xlab("")
  
  p8<-ggplot(FSR %>% filter(Protein.Description==i & Experiment == "Experiment2 ")) +
    geom_jitter(aes(x=Intervention,y=fsr_day,colour=Subject),width=0.2) +
    theme(legend.position = "none") +
    xlab("")
  
  p9<-ggplot(FSR %>% filter(Protein.Description==i & Experiment == "Experiment1 ")) +
    geom_boxplot(aes(x=Intervention,y=fsr_day,colour=Subject)) +
    theme(legend.position = "none") +
    xlab("")
  
  p10<-ggplot(FSR %>% filter(Protein.Description==i & Experiment == "Experiment2 ")) +
    geom_boxplot(aes(x=Intervention,y=fsr_day,colour=Subject)) +
    theme(legend.position = "none") +
    xlab("")
  
  grid.arrange(p7,p8,p9,p10,ncol=2)
  
  
}

dev.off()
