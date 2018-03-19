#for exposure in 9mg DurationVigorousActivity JobInvolvesMainlyWalkingOrStanding NumberOfDaysModeratePhysicalActivity Transition10 X12am_6am X5am_9pm X6pm_12am DurationModerateActivity DurationWalkingForPleasure JobInvolvesManualOrPhysicalWork NumberOfDaysVigorousPhysicalActivity Transition25 X12am_9am X5pm_12am X9pm_5am DurationOfWalks FrequencyWalkingForPleasure NumberDaysWalked10Minutes OverallAccelerationAverage UsualWalkingPace X12pm_6pm X6am_12pm
for exposure in X9mg
do
    for outcome in CvdStatus AliveAt65 AliveAt70
    do
	#sbatch -J "mr$exposure.$outcome" -o logs/$exposure.$outcome.o -e logs/$exposure.$outcome.e -p euan,owners --mem=5000 mr_funcs.sh $exposure $outcome
	sbatch -J "mr$exposure.$outcome" -o logs2/$exposure.$outcome.o -e logs2/$exposure.$outcome.e -p akundaje --mem=10000 mr_funcs.sh $exposure $outcome
    done
done
