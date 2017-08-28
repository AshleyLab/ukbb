for i in 0_1 2000mg 8_9 JobInvolvesShiftWork 10_11 20_21 9_10 NumberDaysWalked10Minutes 11_12 21_22 9mg NumberOfDaysModeratePhysicalActivity 12_13 21 DurationModerateActivity NumberOfDaysVigorousPhysicalActivity 1_2 22_23 DurationOfWalks OverallAccelerationAverage 13_14 23_24 DurationStrenuousSports StandardDeviationOfAcceleration 14_15 2_3 DurationWalkingForPleasure TimeSpentOutdoorsSummer 15_16 3_4 DWT_SMV TimeSpentOutdoorsWinter 16_17 4_5 FrequencyStrenuousSportsLast4Weeks Transition10 17_18 5_6 FrequencyWalkingForPleasure Transition25 18_19 6_7 JobInvolvesHeavyManualOrPhysicalWork UsualWalkingPace 19_20 7_8 JobInvolvesMainlyWalkingOrStanding DurationVigorousActivity
do
#filter by maf
python filter_fuma_input.py --fuma_input $i.top.snps.tsv.fuma.input.txt \
       --maf_filter_upper 0.01 \
       --maf_filter_lower 0.001 \
       --metadata ../../postprocessing/top_snps_activity_merged.txt \
       --outf  $i.maf_filtered_0.001_0.01
#filter by info 
python filter_fuma_input.py --fuma_input $i.top.snps.tsv.fuma.input.txt \
       --info_filter 0.8 \
       --maf_filter_upper 0.01 \
       --maf_filter_lower 0.001 \
       --metadata ../../postprocessing/top_snps_activity_merged.txt \
       --outf $i.info_filtered_0.001_0.01

done
