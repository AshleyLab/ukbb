#for task in 0_1 10_11 11_12 1_2 12_13 13_14 14_15 15_16 16_17 17_18 18_19 19_20 2000mg 20_21 21_22 22_23 2_3 23_24 3_4 4_5 5_6 6_7 7_8 8_9 9_10 9mg DurationModerateActivity DurationOfWalks DurationStrenuousSports DurationVigorousActivity DurationWalkingForPleasure DWT_SMV DWT_SMV1 FrequencyStrenuousSportsLast4Weeks FrequencyWalkingForPleasure JobInvolvesMainlyWalkingOrStanding JobInvolvesManualOrPhysicalWork JobInvolvesShiftWork NumberDaysWalked10Minutes NumberOfDaysModeratePhysicalActivity NumberOfDaysVigorousPhysicalActivity OverallAccelerationAverage StandardDeviationOfAcceleration TimeSpentOutdoorsSummer TimeSpentOutdoorsWinter Transition10 Transition25 UsualWalkingPace
for task in TimeSpentOutdoorsWinter
#for task in X0_1 X1_2 X2_3 X3_4 X4_5 X5_6
do
    sbatch -J "clump.$task" -o logs/clump.$task.o -e logs/clump.$task.e -p euan,owners --mem=20000 aggregate_clump.sh $task 
done
