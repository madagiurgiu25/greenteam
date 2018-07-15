import json
import sys
import pprint
from mongodb.pymongodb import MongoDB

if __name__ == '__main__':
    collection_name = sys.argv[1]
    gene_id = sys.argv[2]

    mongo = MongoDB('minglerna')

    collection = mongo.getCollection(collection_name)
    interactions = mongo.getJsonObjectForGeneID(collection, gene_id)

    read_interactions = json.loads(interactions)
    #print(read_interactions)

    for interaction in read_interactions:
    	print(interaction['transcript_list'])

# python3 getJsonForGeneID.py interactions LNC_GE_mm10_00500789
