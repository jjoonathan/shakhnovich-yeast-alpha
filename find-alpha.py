print "Loading libraries..."
import os, sys, zlib, multiprocessing, urllib2, threading, signal, re, traceback, urllib

# Requires installation
import lxml.etree
import pymongo

species1 = 'cerevisiae'
species2 = 'paradoxus'
run_name = species1 + '-' + species2

sys.stdout.write("Connecting to mongodb...\n")
sys.stdout.flush()
connection = pymongo.Connection('localhost',27017)
genes = connection[run_name]['genes']
genes.ensure_index("ensembl_name", unique=True, sparse=True)
genes.ensure_index("local_name", unique=True)

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
start_workers(8)

def perform_uniprot_search(ensembl_name):
	try:
		url = 'http://www.uniprot.org/uniprot/?query=%s&sort=score&format=xml'%urllib.quote(ensembl_name)
		xmldoc = urllib2.urlopen(url).read()
		try:
			xmltree = lxml.etree.fromstring(xmldoc)
		except Exception as e:
			print "Parsing %s failed.\n"%url
			return True
		uni_ids = xmltree.xpath('//up:entry[@dataset="Swiss-Prot"]/up:accession/text()',namespaces={'up':'http://uniprot.org/uniprot'})
		row = {'ensembl_name':ensembl_name}
		newkeys = {'uniprot_id':uni_ids[0], 'all_uniprot_ids':' '.join(uni_ids)}
		genes.update(row,{'$set':newkeys})
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
		q.put(row['ensembl_name'])
convert_ensembl_to_uniprot()

try:
	q.join()
except KeyboardInterrupt as e:
	print "Ctrl-C Pressed."
	traceback.print_exc(10)
	exit(0)

def create_clustal_in():
	if not os.path.isdir('clustalin'):
		os.mkdir('clustalin')
# create_clustal_in()
