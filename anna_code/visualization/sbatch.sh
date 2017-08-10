#for task in 10_11 11_12 1_2 12_13 13_14 14_15 15_16 16_17 17_18 18_19 19_20 2000mg 20_21 21_22 22_23 2_3 23_24 3_4 4_5 5_6 6_7 7_8 8_9 9_10 9mg DurationModerateActivity DurationOfWalks DurationStrenuousSports DurationWalkingForPleasure DWT_SMV DWT_SMV1 FrequencyStrenuousSportsLast4Weeks FrequencyWalkingForPleasure NumberDaysWalked10Minutes NumberOfDaysModeratePhysicalActivity NumberOfDaysVigorousPhysicalActivity OverallAccelerationAverage StandardDeviationOfAcceleration TimeSpentOutdoorsSummer TimeSpentOutdoorsWinter Transition10 Transition25 UsualWalkingPace JobInvolvesHeavyManualOrPhysicalWork JobInvolvesMainlyWalkingOrStanding JobInvolvesShiftWork
for task in 0_1
do
    #sbatch -J "vis.$task" -o logs/vis.$task.o -e logs/vis.$task.e -p euan,owners --mem=40000 make_inputs_for_manhattan_and_qq.sh $task
    sbatch -J "vis.$task" -o logs/vis.$task.o -e logs/vis.$task.e -p euan,owners --mem=30000 filter_maf.sh $task
done

#for task in JobInvolvesHeavyManualOrPhysicalWork JobInvolvesMainlyWalkingOrStanding JobInvolvesShiftWork
#do
#    sbatch -J "vis.$task" -o logs/vis.$task.o -e logs/vis.$task.e -p euan,owners --mem=40000 make_inputs_for_manhattan_and_qq.categorical.sh $task
#done
