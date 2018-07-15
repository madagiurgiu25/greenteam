#!/bin/sh

python3 insertJsonIntoCollection.py "mirbase_hg38" "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/mirbase/mirbase_hg38.json"
python3 insertJsonIntoCollection.py "mirbase_mm10" "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/mirbase/mirbase_mm10.json"
python3 insertJsonIntoCollection.py "gencode_lnc_hg38" "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/GENCODE/hg38/hg38_long_noncoding.json"
python3 insertJsonIntoCollection.py "gencode_mi_hg38" "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/GENCODE/hg38/hg38_miRNA.json"
python3 insertJsonIntoCollection.py "gencode_lnc_mm10" "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/GENCODE/mm10/mm10_long_noncoding.json"
python3 insertJsonIntoCollection.py "gencode_mi_mm10" "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/GENCODE/mm10/mm10_miRNA.json"
python3 insertJsonIntoCollection.py "noncode_hg38" "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/Noncode/hg38_long_noncoding_noncode.json"
python3 insertJsonIntoCollection.py "noncode_mm10" "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/Noncode/mm10_long_noncoding_noncode.json"
python3 insertJsonIntoCollection.py "lncipedia" "/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/lncipedia/lncipedia_hg38.json"

# db.getCollectionNames().forEach(function(collection) {     print("Collection '" + collection + "' documents: " + db[collection].count());   });
