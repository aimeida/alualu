## g++ -I../libraries trybam.cpp -lz -o trybam
##../debug/trybam /nfs/gpfs/data/Results/GWS/AODNCUD/AODNCUD.bam outputs/AODNCUD

binpath='/nfs_mount/bioinfo/users/yuq/work/Alu/stat/alu/debug/'

input1='/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/'
input2='/nfs/gpfs/data/Results/GWS/'
out1='/nfs_mount/bioinfo/users/yuq/work/Alu/stat/alu/output_test/'
out2='/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/insert_len/'
out3='/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/20_chr0/'

file1=${input2}'AODNCUD/AODNCUD.bam'
pn_file=${input1}PN_noP

if [ $1 -eq "1" ]
then
    debug/utils_debug 0 chr1 |grep -v "^@"
elif [ $1 -eq "2" ]
then
    ${binpath}alu_delete 1 ${input1}config.properties 0 chr0 
elif [ $1 -eq "22" ]
then
    ${binpath}alu_delete 2 ${input1}config.properties
elif [ $1 -eq "23" ]
then
    ${binpath}alu_delete 3 ${input1}config.properties 0
elif [ $1 -eq "3" ]
then
    ${binpath}alu_insert 1 ${input1}config.properties 0 ## step one, print out locations 
elif [ $1 -eq "32" ]
then
    ${binpath}alu_insert 2 ${input1}config.properties 0 ## step one, print out locations 
elif [ $1 -eq "4" ]
then
    echo ${binpath}alu_now 2 ${input1}config.properties 1
elif [ $1 -eq "5" ]
then
    debug/alu_multi test/config.test
elif [ $1 -eq "91" ]
then
    fin='/nfs/gpfs/data/pbstmp/jonsv/alu_multi/results3/BJDGHGJ/InsertLengths.txt'
    fout='/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/grocery/InsertLengths.100_1000_5'
    ${binpath}build_dist ${fin} ${fout} 100 1000 5
    Rscript /nfs_mount/bioinfo/users/yuq/work/Alu/cmds/InsertLenDist.R ${fout}
    # fout='/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/grocery/InsertLengths.100_1000_10'
    # ${binpath}build_dist ${fin} ${fout} 100 1000 10
    # Rscript /nfs_mount/bioinfo/users/yuq/work/Alu/cmds/InsertLenDist.R ${fout}
fi  

