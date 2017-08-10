for task in 0_1 13_14 18_19 22_23 5_6 9mg DWT_SMV JobInvolvesMainlyWalkingOrStanding OverallAccelerationAverage Transition25 10_11 14_15 19_20 2_3 6_7 DurationModerateActivity DWT_SMV1 JobInvolvesShiftWork StandardDeviationOfAcceleration UsualWalkingPace 11_12 15_16 2000mg 23_24 7_8 DurationOfWalks FrequencyStrenuousSportsLast4Weeks NumberDaysWalked10Minutes TimeSpentOutdoorsSummer 1_2 16_17 20_21 3_4 8_9 DurationStrenuousSports FrequencyWalkingForPleasure NumberOfDaysModeratePhysicalActivity TimeSpentOutdoorsWinter 12_13 17_18 21_22 4_5 9_10 DurationWalkingForPleasure JobInvolvesHeavyManualOrPhysicalWork NumberOfDaysVigorousPhysicalActivity Transition10 
do
    sbatch -J "top_snps_for_task.$task" -o logs/top_snps_for_task.$task.o -e logs/top_snps_for_task.$task.e -p euan,owners get_top_snps_per_task.sh $task 
done


for task in JobInvolvesHeavyManualOrPhysicalWork JobInvolvesMainlyWalkingOrStanding JobInvolvesShiftWork
do
    sbatch -J "top_snps_for_task.$task" -o logs/top_snps_for_task.$task.o -e logs/top_snps_for_task.$task.e -p euan,owners get_top_snps_per_task_categorical.sh $task 
done
