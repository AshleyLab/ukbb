outf=open('all_chroms.txt','w')
for i in range(1,23):
    outf.write('chrom'+str(i)+'.bed\t'+'chrom'+str(i)+'.bim\t'+'chrom'+str(i)+'.fam\n')
outf.write('chromX.bed\tchromX.bim\tchromX.fam\n')
outf.write('chromY.bed\tchromY.bim\tchromY.fam\n')
