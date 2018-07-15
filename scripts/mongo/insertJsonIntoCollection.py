import json
import sys
import pprint
from mongodb.pymongodb import MongoDB

def loadJson(df_file):
    with open(df_file) as df_input:    
        df = json.load(df_input)
    return df

if __name__ == '__main__':
    collection_name = sys.argv[1]
    fin_json = sys.argv[2]

    mongo = MongoDB('minglerna')

    collection = mongo.getCollection(collection_name)
    json_bulk = loadJson(fin_json)

    bulkInsert = mongo.insertJSONIntoCollection(collection, json_bulk)

# python3 insertJsonIntoCollection.py mm10_interactions ./mm10/mm10_interactions_allDBs_and_pc.json
# python3 insertJsonIntoCollection.py mm10_assembly ./mm10/mm10_primary_assembly_and_lncRNA.json

# python3 insertJsonIntoCollection.py hg38_interactions ./hg38/Mirbase_human_gencode_lncRNA_filtered_new.json
# python3 insertJsonIntoCollection.py hg38_interactions ./hg38/Mirbase_human_gencode_pc_filtered_new.json
# python3 insertJsonIntoCollection.py hg38_interactions ./hg38/Mirbase_human_lncipedia_filtered.json
# python3 insertJsonIntoCollection.py hg38_interactions ./hg38/Mirbase_human_noncode_filtered.json
# python3 insertJsonIntoCollection.py hg38_assembly ./hg38/hg38_primary_assembly_and_lncRNA.json

# python3 insertJsonIntoCollection.py interactions ./mm10/mm10_interactions_allDBs_and_pc.json
# python3 insertJsonIntoCollection.py interactions ./hg38/Mirbase_human_gencode_lncRNA_filtered_new.json
# python3 insertJsonIntoCollection.py interactions ./hg38/Mirbase_human_gencode_pc_filtered_new.json
# python3 insertJsonIntoCollection.py interactions ./hg38/Mirbase_human_lncipedia_filtered.json
# python3 insertJsonIntoCollection.py interactions ./hg38/Mirbase_human_noncode_filtered.json
# python3 insertJsonIntoCollection.py assembly ./hg38/hg38_primary_assembly_and_lncRNA.json
# python3 insertJsonIntoCollection.py assembly ./mm10/mm10_primary_assembly_and_lncRNA.json