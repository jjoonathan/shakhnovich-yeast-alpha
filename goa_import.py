import sys, os
import pymongo
if len(sys.argv) != 2:
    print "Usage: %s goa_yeast.file"
conn = pymongo.Connection()
goa_yeast = conn['cerevisiae-paradoxus']['goa-yeast']
goa_yeast.ensure_index('uniprot_id')
for line in open(sys.argv[1]):
    line = line.split('\t')
    uniprot_id, go_tag = line[1], line[4]
    goa_yeast.update({'uniprot_id':uniprot_id}, {'go_tag':go_tag}, upsert=True)
    
