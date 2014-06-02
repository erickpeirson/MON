# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from nltk.corpus import stopwords
import string
import rdflib
import urllib2
import xml.etree.ElementTree as ET
import os
import subprocess
import pickle
import matplotlib.pyplot as plt

# <markdowncell>

# ## Relevant paths & services

# <codecell>

orgpath = '/Users/erickpeirson/Model Organisms in Neuroscience/data'

# <codecell>

datapath = orgpath + '2014-01-18 RDF'

# <codecell>

taxonomy_api = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy"
taxonomy_api_get = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy"

# <markdowncell>

# ## Methods

# <codecell>

def find_names(paper, linnaeus_path="/Applications/linnaeus/bin/linnaeus-2.0.jar", tmp_path="/tmp/"):
    """
    Yields a NCBI Taxonomy id for organisms, via LINNAEUS NER.

    Calls the LINNAEUS NER engine to find organism names (common & scientific) in the
    paper title and abstract.
    """

    with open(tmp_path+"text.txt", 'w') as f:
        f.write(str(strip_non_ascii(paper.title)) + "\n")
        f.write(str(strip_non_ascii(paper.abstract)) + "\n")
    
    subprocess.call(['java', '-jar', linnaeus_path, '--text', tmp_path+"text.txt", '--out', tmp_path+"results.tsv" ])
    
    with open(tmp_path+"results.tsv", 'r') as f:
        results = [ line.strip('\n').split('\t') for line in f.readlines() ]
    
    return set([ r[0].split(':')[-1] for r in results[1:] ])

# <codecell>

# Cache for classify()
classified = {}

# <codecell>

def classify(taxon):
    """
    Yields the higher-level taxonomic group to which taxon belongs.

    Retrieves data for taxon from NCBI Taxonomy, then passes the result
    to merge() to classify. Uses classified (dict) to cache results.
    """
    try:
        return classified[taxon]
    except KeyError:
        response = urllib2.urlopen(taxonomy_api_get + "&id=" + taxon).read()
        xml = ET.fromstring(response)
        merged = merge(xml)
        classified[taxon] = merged 
        return merged

# <codecell>

def merge(result):
    """
    Uses NCBI Taxonomy XML response to decide to which higher-level
    taxonomic group a taxon belongs.
    """
    
    taxa = [ tax.text for tax in result.findall('.//TaxId') ]
    if '89593' in taxa: # Craniata
        if '7898' in taxa or '7894' in taxa or '7878' in taxa: # Fishes
            return 'Fish'
        elif '8292' in taxa: # Amphibeans
            return 'Amphibean'
        elif '40674' in taxa: # Mammals
            if '9681' in taxa: # Cats
                return 'Cat'
            elif '10066' in taxa: # Mice & rats
                return 'Mouse or Rat'
            elif '9443' in taxa: # Primates
                if '9605' in taxa: # Human
                    return 'Human'
                elif '9539' in taxa: # Macaque
                    return 'Macaque'
                else:
                    return 'Other Primate'
            else: # Other mammals
                return 'Other Mammal'
        # Remainder
        return 'Other Craniata'
    else: # Not in Craniata
        return 'Non-Craniata'

# <codecell>

class paper:
    """
    Describes a minimal PubMed entry.
    """
    def __init__(self, uri, abstract, journal, title):
        self.uri = uri
        self.abstract = abstract
        self.journal = journal
        self.title = title

# <markdowncell>

# ## General utility methods

# <codecell>

def tokenize(s):
    s = s.translate(string.maketrans("",""), string.punctuation)
    return s.lower().split()

# <codecell>

def remove_stopwords(t):
    return [ word for word in t if word not in stopwords.words() and word != 'ab' ]

# <codecell>

def strip_non_ascii(string):
    ''' Returns the string without non ASCII characters'''
    stripped = (c for c in string if 0 < ord(c) < 127)
    return ''.join(stripped)

# <codecell>

def remove_integers(t):
    new_t = []
    for word in t:
        try:
            int(word)
        except ValueError:
            new_t.append(word)
    return new_t

# <codecell>

class dictionary:
    """
    A two-way index for integer/string pairs.
    """
    def __init__(self):
        self.by_str = {}
        self.by_int = {}
        
    def __setitem__(self, key, value):
        if type(key) == str:
            self.by_str[key] = value
            self.by_int[value] = key
        if type(key) == int:
            self.by_int[key] = value
            self.by_str[value] = key
    
    def __getitem__(self, key):
        if type(key) == str:
            return self.by_str[key]
        if type(key) == int:
            return self.by_int[key]

# <markdowncell>

# ## Parse RDF

# <codecell>

g = rdflib.Graph()
for filename in os.listdir(datapath):
    if filename.split('.')[-1] == 'rdf':
        g.load(datapath+"/"+filename)

# <markdowncell>

# ### For convenience, index RDF triples.

# <codecell>

asdict = {}    # Index by subject.
for s,p,o in g:  # Subject, predicate, and object in RDF triple.
    try:
        asdict[s][p] = o
    except KeyError:
        asdict[s] = {}
        asdict[s][p] = o

# <markdowncell>

# ### Generate a list of papers (see class `paper`, above)

# <codecell>

aslist = []
for key, value in asdict.iteritems():
    type_key = rdflib.term.URIRef(u'http://www.w3.org/1999/02/22-rdf-syntax-ns#type')
    article_key = rdflib.term.URIRef(u'http://purl.org/net/biblio#Article')
    if value[type_key] == article_key:
        try:
            p = paper(key, value[rdflib.term.URIRef(u'http://purl.org/dc/terms/abstract')], value[rdflib.term.URIRef(u'http://purl.org/dc/terms/isPartOf')], value[rdflib.term.URIRef(u'http://purl.org/dc/elements/1.1/title')])
            aslist.append(p)
        except KeyError:
            pass

# <codecell>

len(aslist) # number of abstracts

# <markdowncell>

# ## Generate a complete vocabulary for the set of `paper` abstracts
# * Remove all integers
# * Remove stopwords (NLTK stopwords)
# * Calculate overall term frequencies ( `termcounts` )
# * Index the vocabulary ( `mo_dict` )

# <codecell>

vocabulary = set([])
for p in aslist:
    s = str(strip_non_ascii(p.abstract))
    words = set(tokenize(s))
    vocabulary = vocabulary | words

# <codecell>

vocabulary = remove_integers(list(vocabulary))

# <codecell>

vocabulary = set(vocabulary) - set([word for word in stopwords.words()])

# <codecell>

len(vocabulary)

# <codecell>

termcounts = {}
for p in aslist:
    for term in tokenize(str(strip_non_ascii(p.abstract))):
        try:
            termcounts[term] += 1
        except KeyError:
            termcounts[term] = 1

# <codecell>

mo_dict = dictionary()
i = 0
for word in list(vocabulary):
    mo_dict[word] = i
    i += 1

# <markdowncell>

# ## First attempt: dictionary look-up for each term, using NCBI Taxonomy (slow)

# <codecell>

species = []
for word in list(vocabulary):
    found = False
    response = urllib2.urlopen(taxonomy_api + "&term=" + word).read()
    xml = ET.fromstring(response)
    for e in xml:
        if e.tag == 'Count':
            if int(e.text) > 0:
                found = True
                name = word
        if e.tag == 'IdList' and len(e) > 0:
            ident = e[0].text
    if found:
        response = urllib2.urlopen(taxonomy_api_get + "&id=" + ident).read()
        xml = ET.fromstring(response)
        for e in xml[0]:
            if e.tag == 'Rank':
                rank = e.text
        species.append( (name, rank, ident) )

# <codecell>

species_set = set([ s for s,r,i in species ])

# <codecell>

counts = {}
for org,rank,ident in species:
    counts[org] = 0
    for p in aslist:
        if org in p.abstract:
            counts[org] += 1

# <codecell>

fig = plt.figure(figsize(40, 10))
ax = fig.add_subplot(1, 1, 1)

names = [ k for k,v in counts.iteritems() if v > 10 ]
values = [ v/float(len(aslist)) for v in counts.values() if v > 10 ]
N = len(names)
ind = range(N)
ax.bar(ind, values, align='center')
ax.set_xticks(ind)
ax.set_xticklabels(names)
fig.autofmt_xdate()
plt.show()

# <markdowncell>

# ## Second attempt: LINNAEUS NER & merging via NCBI Taxonomy

# <codecell>

paper_orgs = {}
journal_orgs = {}
for p in aslist:
    r = find_names(p)
    paper_orgs[p.uri] = r
    try:
        journal_orgs[p.journal]
    except KeyError:
        journal_orgs[p.journal] = {}
    for i in r:
        try:
            journal_orgs[p.journal][i] += 1
        except KeyError:
            journal_orgs[p.journal][i] = 1

# <codecell>

by_journal = {}
for uri, values in journal_orgs.iteritems():
    rels = { p:v for p,v in g[uri] }
    title = str(rels[rdflib.term.URIRef(u'http://purl.org/dc/elements/1.1/title')])
    try:
        by_journal[title]
    except KeyError:
        by_journal[title] = {}
    for k,v in values.iteritems():
        try:
            by_journal[title][k] += v
        except KeyError:
            by_journal[title][k] = v

# <codecell>

overall_taxa = {}
for paper,orgs in paper_orgs.iteritems():
    for org in orgs:
        taxon = classify(org)
        try:
            overall_taxa[taxon] += 1
        except KeyError:
            overall_taxa[taxon] = 1

# <codecell>

journal_taxa = {}
for journal, orgs in by_journal.iteritems():
    journal_taxa[journal] = {}
    N = np.sum([c for c in orgs.values()])
    for org,count in orgs.iteritems():
        taxon = classify(org)
        try:
            journal_taxa[journal][taxon] += float(count)/N
        except KeyError:
            journal_taxa[journal][taxon] = float(count)/N

# <codecell>

for journal, values in journal_taxa.iteritems():
    font = { 'size'   : 16}
    matplotlib.rc('font', **font)    
    
    fig = plt.figure(figsize=(20,10), dpi=300)
    ind = np.array(range(len(values)))
    plt.bar(ind, np.array(values.values())*100)
    plt.title(journal)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xticks(ind+0.5)
    ax.set_xticklabels(values.keys())
    plt.ylabel("Percent of Papers")
    fig.autofmt_xdate(rotation=90)
    plt.savefig("/Users/erickpeirson/Model Organisms in Neuroscience/figures/"+journal.lower().replace(' ','')+".png")

# <codecell>

fig = plt.figure(figsize=(20,10), dpi=300)
ind = np.array(range(len(overall_taxa)))
plt.bar(ind, overall_taxa.values())
ax = fig.add_subplot(1, 1, 1)
ax.set_xticks(ind+0.5)
ax.set_xticklabels(overall_taxa.keys())
fig.autofmt_xdate(rotation=90)
plt.title("All Journals")
plt.ylabel("Number of Papers")
plt.savefig("/Users/erickpeirson/Model Organisms in Neuroscience/figures/alljournals.png")

# <markdowncell>

# ### Persistency

# <codecell>

with open(orgpath+"/organisms.txt", "w") as f:
    for n,r,i in species:
        f.write(n + "\t" + r + "\t" + i + "\n")

# <codecell>

with open(orgpath+"/organisms.txt", "r") as f:
    species = [ tuple(line.strip('\n').split('\t')) for line in f.readlines() ]

# <codecell>

with open("/Users/erickpeirson/Model Organisms in Neuroscience/data/paper_orgs.pickle", "w") as f:
    pickle.dump(paper_orgs, f)
with open("/Users/erickpeirson/Model Organisms in Neuroscience/data/journal_orgs.pickle", "w") as f:
    pickle.dump(journal_orgs, f)
with open("/Users/erickpeirson/Model Organisms in Neuroscience/data/by_journal.pickle", "w") as f:
    pickle.dump(by_journal, f)
with open("/Users/erickpeirson/Model Organisms in Neuroscience/data/overall_taxa.pickle", "w") as f:
    pickle.dump(overall_taxa, f)
with open("/Users/erickpeirson/Model Organisms in Neuroscience/data/journal_taxa.pickle", "w") as f:
    pickle.dump(journal_taxa, f)

# <markdowncell>

# ## Modeling (ignore for now)

# <codecell>

from vsm.corpus import Corpus
from vsm.corpus.util import corpusbuilders as cb
from vsm.model import tf, tfidf, lsa, ldagibbs, beaglecomposite, beaglecontext, beagleenvironment, beagleorder, ldacgsmulti
from vsm.viewer import ldagibbsviewer, tfidfviewer
from vsm.viewer import beagleviewer
import pickle

# <codecell>

for i in xrange(len(aslist)):
    aslist[i].abstract = tokenize(str(strip_non_ascii(aslist[i].abstract)))

# <codecell>

abstracts = [ tokenize(str(strip_non_ascii(p.abstract))) for p in aslist]

# <codecell>

model_path = "/Users/erickpeirson/Desktop/MO"
filename = "ACSChemicalNeuroscience"

# <codecell>

ec = cb.corpus_fromlist(abstracts, context_type='sentence')
#pickle.dump(ec, open(model_path+filename+"_ec.model", "w"))
print str(len(ec.words)) + " in unfiltered vocabulary"

cc = ec.apply_stoplist(stopwords.words())
#pickle.dump(cc, open(model_path+filename+"_cc.model", "w"))
print str(len(cc.words)) + " after stoplist applied"
    
be = beagleenvironment.BeagleEnvironment(ec)
be.train()
#pickle.dump(be, open(model_path+filename+"_be.model", "w"))
print "environment done"
    
ms = beaglecontext.BeagleContextSeq(cc, ec, be.matrix)
ms.train()
#pickle.dump(ms, open(model_path+filename+"_ms.model", "w"))
print "context done"
    
#od = beagleorder.BeagleOrderSeq(ec, be.matrix)
#od.train()
#pickle.dump(od, open(model_path+filename+"_od.model", "w"))
#print "order done"

#bc = beaglecomposite.BeagleComposite(cc, ms.matrix, ec, od.matrix)
#bc.train()
#pickle.dump(bc, open(model_path+filename+"_bc.model", "w"))
    
#bv = beagleviewer.BeagleViewer(cc, bc)
bv_x = beagleviewer.BeagleViewer(cc, ms)
#bv_o = beagleviewer.BeagleViewer(ec, od)

# <codecell>

bv_x.sim_word_word(['mice'])

# <codecell>

print bv_x.simmat_words(['mouse', 'rat', 'human', 'alzheimers', 'schizophrenia', 'model'])

# <codecell>

'alzheimers' in vocabulary

# <codecell>

bv_o.sim_word_word(['mice'])

# <codecell>

bv.sim_word_word(['mice'])

# <codecell>

bv_x.sim_word_word(['rabbits'])

# <codecell>

bv_x.sim_word_word(['monkeys'])

# <codecell>

bv_x.sim_word_word(['mouse'])

# <codecell>

l = ldacgsmulti.LdaCgsMulti(cc, 'sentence', K=20)

# <codecell>

l.train(itr=100, n_proc=8)

# <codecell>


