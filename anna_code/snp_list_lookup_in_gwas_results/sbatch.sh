#for phenotype in 0_1 15_16 21_22 6_7  DurationVigorousActivity JobInvolvesManualOrPhysicalWork TimeSpentOutdoorsSummer X12pm_6pm X9pm_5am 10_11 16_17 22_23 7_8  DurationWalkingForPleasure JobInvolvesShiftWork TimeSpentOutdoorsWinter X5am_9pm 11_12 17_18 2_3 8_9  DWT_SMV  NumberDaysWalked10Minutes Transition10 X5pm_12am 1_2 18_19 23_24 9_10 DWT_SMV1  NumberOfDaysModeratePhysicalActivity Transition25 X6am_12pm 12_13 19_20 3_4 DurationModerateActivity FrequencyStrenuousSportsLast4Weeks NumberOfDaysVigorousPhysicalActivity UsualWalkingPace X6pm_12am 13_14 2000mg 4_5 DurationOfWalks FrequencyWalkingForPleasure OverallAccelerationAverage X12am_6am X9am_5pm 14_15 20_21 5_6 DurationStrenuousSports JobInvolvesMainlyWalkingOrStanding StandardDeviationOfAcceleration X12am_9am X9mg 
#for phenotype in JobInvolvesManualOrPhysicalWork JobInvolvesShiftWork JobInvolvesMainlyWalkingOrStanding
for phenotype in 9mg
do
    for chrom in `seq 1 22`
    do
        sbatch -J "$phenotype.$chrom.continuous" -o logs/$phenotype.$chrom.o -e logs/$phenotype.$chrom.e -p akundaje,euan,owners  lookup_snp_list.sh $phenotype $chrom
	#sbatch -J "$phenotype.$chrom.continuous" -o logs/$phenotype.$chrom.o -e logs/$phenotype.$chrom.e -p akundaje,euan,owners  lookup_snp_list.categorical.sh $phenotype $chrom
    done
done
