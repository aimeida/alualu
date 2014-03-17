fin=$1
fout=$2

for ci in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y 
#for ci in 1 
do
    chrn=chr${ci} 
    less ${fin}.${chrn} | sort -k4n > ${fout}.${chrn}
done