print "Loading libraries..."
import os, sys, zlib, Queue, urllib2, threading, signal, re
signal.signal(signal.SIGINT, lambda a,b: exit(0))

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

q = Queue.Queue()
number_enqueued_tasks = 0 # We have to count this ourselves, WTF
number_tasks_done = 0
def start_workers(num_workers):
	sys.stdout.write("Starting workers...\n")
	sys.stdout.flush()
	def worker():
		global number_enqueued_tasks
		global number_tasks_done
		while True:
			fractn = number_tasks_done * 100. / (number_enqueued_tasks+.01)
			print "\x1B[1FProcessing queue: %i%% (%i/%i)"%(int(fractn), number_tasks_done, number_enqueued_tasks)
			item = q.get()
			success = item()
			q.task_done()
			if success:
				number_tasks_done += 1
			else:
				q.put(item)
	for i in range(num_workers):
		t = threading.Thread(target=worker)
		t.daemon = True
		t.start()
start_workers(8)

def convert_ensembl_to_uniprot():
	def perform_uniprot_search(ensembl_name):
		try:
			print "a"
			xmldoc = urllib2.urlopen('http://www.uniprot.org/uniprot/?query=%s&sort=score&format=xml'%ensembl_name).read()
			print "b"
			xmltree = lxml.etree.fromstring(xmldoc)
			print "c"
			uni_ids = xmltree.xpath('//up:entry[@dataset="Swiss-Prot"]/up:accession/text()',namespaces={'up':'http://uniprot.org/uniprot'})
			print "d"
			row = {'ensembl_name':ensembl_name}
			print "upd"
			genes.update(row,{'$set':{'uniprot_id':uni_ids[0], 'all_uniprot_ids':' '.join(uni_ids)}})
		except Exception as e:
			print "Exception searching for %s: %s"%(ensembl_name,e)
			sys.stdout.flush()
			return False
		sys.stdout.flush()
		return True
	print "Converting ensembl to uniprot..."
	for row in genes.find({'uniprot_id':None}):
		global number_enqueued_tasks
		number_enqueued_tasks += 1
		q.put(lambda: perform_uniprot_search(row['ensembl_name']))
convert_ensembl_to_uniprot()
q.join()

def create_clustal_in():
	if not os.path.isdir('clustalin'):
		os.mkdir('clustalin')
# create_clustal_in()
