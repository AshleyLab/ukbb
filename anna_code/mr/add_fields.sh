#Add age, age^2, primary cause of death 
#convert categorical primary cause of death variable to an ordinal 
python add_fields.py --main_f covariates.txt \
--source_f augmented_ages_recruitment_death.txt primary_cause_of_death_collapsed.txt \ 
--source_fields AgeAtRecruitment AgeAtRecruitmentSquared UnderlyingPrimaryCauseOfDeath \ 
--outf covariates_augmented.txt

#add survival time as an outcome 
python add_fields.py --main_f mr_outcome.euro.txt \
--source_f augmented_ages_recruitment_death.txt \ 
--source_fields SurvivalTime \ 
--outf mr_outcome.euro.augmented.txt

