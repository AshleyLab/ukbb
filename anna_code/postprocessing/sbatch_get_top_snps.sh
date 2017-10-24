for task in 0_1 14_15 2000mg 3_4 9_10 DurationWalkingForPleasure JobInvolvesHeavyManualOrPhysicalWork JobInvolvesShiftWork_linear TimeSpentOutdoorsSummer 10_11 15_16 20_21 4_5 9mg DWT_SMV JobInvolvesHeavyManualOrPhysicalWork_linear NumberDaysWalked10Minutes TimeSpentOutdoorsWinter 11_12 16_17 21_22 5_6 DurationModerateActivity DWT_SMV1 JobInvolvesMainlyWalkingOrStanding NumberOfDaysModeratePhysicalActivity Transition10 1_2 17_18 22_23 6_7 DurationOfWalks FrequencyStrenuousSportsLast4Weeks JobInvolvesMainlyWalkingOrStanding_linear NumberOfDaysVigorousPhysicalActivity Transition25 12_13 18_19 2_3 7_8 DurationStrenuousSports FrequencyWalkingForPleasure JobInvolvesManualOrPhysicalWork OverallAccelerationAverage UsualWalkingPace 13_14 19_20 23_24 8_9 DurationVigorousActivity 
do
    sbatch -J "top_snps_for_task.$task" -o logs/top_snps_for_task.$task.o -e logs/top_snps_for_task.$task.e -p akundaje,owners get_top_snps_per_task.sh $task 
done


for task in JobInvolesMainlyWalkingOrStanding JobInvolvesShiftWork StandardDeviationOfAcceleration
do
    sbatch -J "top_snps_for_task.$task" -o logs/top_snps_for_task.$task.o -e logs/top_snps_for_task.$task.e -p akundaje,owners get_top_snps_per_task_categorical.sh $task 
done
