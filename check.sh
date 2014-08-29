## sh check.sh 0 1473 chr1 109573232
## sh check.sh 1 1473 chr1 109573232

version=$1
family=$2
chrn=$3      ##"chr1"
pos1=$4

if [ $version = 0 ]; then
    path0="/home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/fixed_delete0/" 
    fout0="/home/qianyuxx/faststorage/Alu/tmpi"
fi

if [ $version = 1 ]; then
    path0="/home/qianyuxx/faststorage/AluDK/outputs/delete_alu0"
    fout0="/home/qianyuxx/faststorage/Alu/tmpd"
fi

echo "############# checking " $path0   > $fout0
pos2=$(( ${pos1} -1 ))

cd $path0
fn1=`ls *.vcf`

echo `grep ${pos1} ${fn1}`              >> $fout0

echo "############# checking tmp1"      >> $fout0

for fn in `ls tmp1s/${family}*.tmp1`
do
    echo $fn                            >> $fout0
    echo `grep ${pos2} ${fn}`           >> $fout0
done


echo "############# checking tmp2"      >> $fout0

for fn in `ls tmp2s/${family}*.tmp2`
do
    echo $fn                            >> $fout0
    echo `grep ${pos2} ${fn}`           >> $fout0
done

echo "output to" ${fout0}