print "Loading libraries..."
import os, sys, zlib, multiprocessing, threading, signal, re, traceback, urllib, eventlet, httplib

# Requires installation
import lxml.etree
import pymongo
from eventlet.green import urllib2

species1 = 'cerevisiae'
species2 = 'paradoxus'
run_name = species1 + '-' + species2

sys.stdout.write("Connecting to mongodb...\n")
sys.stdout.flush()
connection = pymongo.Connection('localhost',27017)
genes = connection[run_name]['genes']
ortholog_groups = connection[run_name]['ortholog_groups']
# genes.ensure_index("ensembl_name", unique=True, sparse=True)
genes.ensure_index("local_name")
genes.ensure_index("uniprot_id")
ortholog_groups.ensure_index("cerevisiae")

cell_idx_to_species = ['cerevisiae', 'paradoxus', 'mikatae', 'bayanus']
def update_orthology():
	print "Updating orthology..."
	for line in open('data/listing.txt'):
		cells = line.split("\t")
		for i in range(len(cells)):
			try:
				local_name, ensembl_name = cells[i].split(' ')
			except:
				continue
			species = cell_idx_to_species[i]
			row = {'local_name':local_name, 'species':species, 'ensembl_name':ensembl_name}
			genes.update(row, {'$set':row}, upsert=True)
# update_orthology()

def archive_sequences():
	sys.stdout.write("Archiving sequences...\n\n")
	sys.stdout.flush()
	for species in cell_idx_to_species:
		gene_name = None
		gene_seq = ''
		for line in open('data/%s.fasta'%species):
			if line[0]=='>' and gene_seq!='':
				row  = {'local_name':gene_name}
				gene_name = re.match('>(\S+)\s+(\S+)',line).groups()[0]
				m = re.match('^([atgcN.]*)(.*?)([atgcN.]*)$',gene_seq)
				if m==None:
					print "MATCH FAILED: %s vs %s"%(gene_name,gene_seq)
					gene_seq=''
					continue
				gene_seq=''
				pre, orf, post = m.groups()
				genes.update(row, {'$set':{'pre':pre,'post':post,'orf':orf}})
			elif line[0]=='>' and gene_seq=='':
				gene_name = re.match('>(\w+)\s+(\w+)',line).groups()[0]
			else:
				gene_seq = gene_seq + line.strip()
# archive_sequences()

q = multiprocessing.JoinableQueue()
number_enqueued_tasks = 0 # We have to count this ourselves, WTF
number_tasks_done = 0
def start_workers(num_workers):
	sys.stdout.write("Starting workers...\n")
	sys.stdout.flush()
	def worker():
		try:
			global number_enqueued_tasks
			global number_tasks_done
			while True:
				fractn = number_tasks_done * 100. / (number_enqueued_tasks+.01)
				print "\x1B[1FProcessing queue: %i%% (%i/%i)"%(int(fractn), number_tasks_done, number_enqueued_tasks)
				item = q.get()
				success = perform_uniprot_search(item)
				q.task_done()
				if success:
					number_tasks_done += 1
				else:
					q.put(item)
		except EOFError:
			return
	for i in range(num_workers):
		t = threading.Thread(target=worker)
		t.daemon = True
		t.start()
def joinq():
	try:
		q.join()
	except KeyboardInterrupt as e:
		print "Ctrl-C Pressed."
		traceback.print_exc(10)
		exit(0)

def perform_uniprot_search(row):
	try:
		url = 'http://www.uniprot.org/uniprot/?query=%s&sort=score&format=xml'%urllib.quote(row['ensembl_name'])
		xmldoc = urllib2.urlopen(url).read()
		try:
			xmltree = lxml.etree.fromstring(xmldoc)
		except Exception as e:
			print "Parsing %s failed.\n"%url
			return True
		uni_ids = xmltree.xpath('//up:entry[@dataset="Swiss-Prot"]/up:accession/text()',namespaces={'up':'http://uniprot.org/uniprot'})
		row = {'_id':row['_id']}
		newkeys = {'uniprot_id':uni_ids[0], 'all_uniprot_ids':' '.join(uni_ids)}
		genes.update(row, {'$set':newkeys}, multi=True)
	except Exception as e:
		print "Exception searching for %s: %s"%(ensembl_name,e)
		traceback.print_exc(10)
		sys.stdout.flush()
		return False
	sys.stdout.flush()
	return True
def convert_ensembl_to_uniprot():
	print "Converting ensembl to uniprot..."
	for row in genes.find({'uniprot_id':None}):
		global number_enqueued_tasks
		number_enqueued_tasks += 1
		q.put(row)
	start_workers(9)
	joinq()
# convert_ensembl_to_uniprot()



def import_roundup_clusters():
	print "Importing roundup clusters..."
	xmldoc = open('data/roundup_cerevisiae_clusters.xml').read()
	xmltree = lxml.etree.fromstring(xmldoc)
	xmlns = {'ox':'http://orthoXML.org/2011/'}
	gene_id_to_uniprot_id = {}
	print "\tCreating local id map..."
	for species in xmltree.xpath('//ox:species',namespaces=xmlns):
		species_name = species.attrib['name']
		for gene in species.xpath('ox:database/ox:genes/ox:gene',namespaces=xmlns):
			gene_id_to_uniprot_id[int(gene.attrib['id'])] = gene.attrib['protId']
	print "\tProcessing ortholog groups..."
	for orthologGroup in xmltree.xpath('//ox:orthologGroup',namespaces=xmlns):
		gene_ref_to_uniprot_id = lambda geneRef: gene_id_to_uniprot_id[int(geneRef.attrib['id'])]
		member_uniprot_ids = map(gene_ref_to_uniprot_id, orthologGroup.xpath('ox:geneRef',namespaces=xmlns))
		a_member = member_uniprot_ids[0]
		newrow = {'cerevisiae':member_uniprot_ids}
		ortholog_groups.update({'cerevisiae':a_member}, {'$set':newrow}, upsert=True)
		for uniprot_id in member_uniprot_ids:
			newrow = {'uniprot_id':uniprot_id, 'species':'cerevisiae'}
			genes.update({'uniprot_id':uniprot_id}, {'$set':newrow}, upsert=True)
# import_roundup_clusters()



def fetch_gene(g):  # g is a uniprot identifier
    fetchurl = 'http://www.uniprot.org/uniprot/%s.xml'%g
    try:
	response = urllib2.urlopen(fetchurl).read()
	xmldoc = lxml.etree.fromstring(response)
	xmlns = {'up':'http://uniprot.org/uniprot'}
	seq_refs = xmldoc.xpath('//up:dbReference[@type="EMBL"]/up:property[@type="protein sequence ID"]/@value',namespaces=xmlns)
	if len(seq_refs)==0:
	    print "Error: %s had no associated sequence references."%g
	    return "FAIL"
    except (urllib2.HTTPError, httplib.BadStatusLine) as e:
	print "Error '%s'->%s"%(fetchurl,e)
	return "RETRY"
    try:
	url = "http://www.ebi.ac.uk/ena/data/view/%s&display=fasta"%(seq_refs[0],)
	return urllib2.urlopen(url).read()
    except (urllib2.HTTPError, httplib.BadStatusLine) as e:
	print "Error2 %s %s"%(url, e)
	return "RETRY"
def run_gene_fetch(gene):  # gene is a db row
	fasta = "RETRY"
	for i in range(5):
		fasta = fetch_gene(gene['uniprot_id'])
		if fasta != "RETRY":
			break
	if fasta == "FAIL":
		return False
	genes.update({'_id':gene['_id']}, {'$set':{'fasta':fasta}})
	return True
def fetch_new_genes():
	print "Fetching new genes...\n"
	e, dne = {'$exists':True}, {'$exists':False}
	genesToFetch = genes.find({'uniprot_id':e, 'orf':dne, 'fasta':dne}, {'uniprot_id':1})
	num = genesToFetch.count()
	downloader_pool = eventlet.greenpool.GreenPool(size=8)
	i=0
	for result in downloader_pool.imap(run_gene_fetch, genesToFetch):
		i += 1
		print "\x1B[1F%i%% (%i/%i)"%(int(i*100./num),i,num)
# fetch_new_genes()
	
def create_fasta_for_old_genes():  # SMART-RESTART
	print "Adding fasta to old genes..."
	e, dne = {'$exists':True}, {'$exists':False}	
	genes_to_update = genes.find({'uniprot_id':e, 'orf':e, 'fasta':dne}, {'orf':1, 'local_name':1, 'uniprot_id':1})
	num = genes_to_update.count()
	i = 0
	for gene in genes_to_update:
		fasta = "> %s %s\n%s"%(gene['local_name'], gene['uniprot_id'], gene['orf'])
		genes.update({"_id":gene['_id']}, {'$set':{'fasta':fasta}})
#create_fasta_for_old_genes()

def add_non_cerevisiae_genes_to_ortholog_groups():
	print "Adding non-cerevisiae genes to ortholog groups..."
	for gene in genes.find({'species':{'$ne':'cerevisiae'}}, {'uniprot_id':1,'local_name':1,'species':1}):
		if 'uniprot_id' not in gene:
			print "No uniprot for %s, aborting.\n"%gene['local_name']
			continue
		ortholog_groups.update({'cerevisiae':gene['uniprot_id']}, {'$set':{gene['species']:gene['local_name']}})
#add_non_cerevisiae_genes_to_ortholog_groups()

def create_clustal_ids():
	if ortholog_groups.find({'clustal_id':{'$exists':False}}).count() == 0:
		return
	print "Creating clustal ids..."
	i = 1
	for group in ortholog_groups.find():
		ortholog_groups.update({"_id":group['_id']}, {'$set':{'clustal_id':i}})
		i += 1
#create_clustal_ids()

def create_clustal_in():
	print "Writing clustal input files..."
	for species2 in ['paradoxus', 'mikatae', 'bayanus']:
		indir = 'clustalin-%s'%species2
		if not os.path.isdir(indir):
			os.mkdir(indir)
		for orthogrp in ortholog_groups.find({species2:{'$exists':True}}):
			def get_fasta(uniprot_id):
				gene = genes.find_one({'uniprot_id':uniprot_id})
				if gene == None:
					return None
				return gene['fasta']
			cerevisiae_genes = map(get_fasta, orthogrp['cerevisiae'])
			cerevisiae_genes = [x for x in cerevisiae_genes if x != None]
			outgroup_gene = genes.find_one({'local_name':orthogrp[species2]})['fasta']
			if outgroup_gene==None or len(cerevisiae_genes)==0:
				print "Incomplete gene set for %s\n"%orthogrp['clustal_id']
				continue
			ingenes = list(cerevisiae_genes)
			ingenes.append(outgroup_gene)
			infile = open(indir+'/%i.fasta'%orthogrp['clustal_id'],'w')
			infile.write('\n'.join(ingenes))
			infile.close()
#create_clustal_in()

def read_clustalout():
	for species2 in ['paradoxus', 'mikatae', 'bayanus']:
		outdir = 'clustalout-'+species2
		mkout
		for fname in os.listdir(outdir):
			
