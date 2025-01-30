#ID
head -n 1 "/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/81499/released/2022-06-07/data/data.51913.csv" | sed 's/,/\n/g' | cat -n | grep 'eid'
#Sex
head -n 1 "/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/81499/released/2022-06-07/data/data.51913.csv" | sed 's/,/\n/g' | cat -n | grep '31-0.0'
#Age
head -n 1 "/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/81499/released/2022-06-07/data/data.51913.csv" | sed 's/,/\n/g' | cat -n | grep '21022'
#Smokignstatus 
head -n 1 "/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/81499/released/2022-06-07/data/data.51913.csv" | sed 's/,/\n/g' | cat -n | grep '20160-0.0'
#FEV1
head -n 1 "/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/81499/released/2022-06-07/data/data.51913.csv" | sed 's/,/\n/g' | cat -n | grep '20256-0.0'


cut -d',' -f 1,23,24,9034,9603 "/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/81499/released/2022-06-07/data/data.51913.csv" > /mnt/storage/scratch/bu19525/BristolPhD/MVMR_follow_up/Data/SmoStatFEV1.txt
