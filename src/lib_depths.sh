#!/bin/bash

project_dir="/home/kh593/project/nfkb_seq"
atac_table="${project_dir}/data/atac_libs.tsv"
mint_table="${project_dir}/data/mint_libs.tsv"
log_dir="/home/kh593/scratch60/nfkb_seq/logs"

echo -n > ${project_dir}/data/atac_depths.csv

count=1
while read donor expt stim lib
do
if [ $count == 1 ]
then
echo -e "lib\tdepth" >> ${project_dir}/data/atac_depths.csv
else
cat ${log_dir}/callpeaks_atac/peakcall_atac${count}.out | awk -v library=$lib \
'BEGIN{ OFS="\t"; depth=0 }
/^Fullset size: [[:digit:]]+/ { depth=$3 }
END{ print library,depth }' >> ${project_dir}/data/atac_depths.csv
fi
    
(( count++ ))
done < ${atac_table}



## The same thing for mintchip reads 
echo -n > ${project_dir}/data/mint_depths.csv

count=1
while read donor expt stim lib
do
if [ $count == 1 ]
then
echo -e "lib\tdepth" >> ${project_dir}/data/mint_depths.csv
else
cat ${log_dir}/callpeaks_mint/peakcall_mint${count}.out | awk -v library=$lib \
'BEGIN{ OFS="\t"; depth=0 }
/^Fullset size: [[:digit:]]+/ { depth=$3 }
END{ print library,depth }' >> ${project_dir}/data/mint_depths.csv
fi
    
(( count++ ))
done < ${mint_table}
