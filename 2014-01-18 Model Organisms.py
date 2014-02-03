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

# <codecell>

def find_names(paper, linnaeus_path="/Applications/linnaeus/bin/linnaeus-2.0.jar", tmp_path="/tmp/"):
    with open(tmp_path+"text.txt", 'w') as f:
        f.write(str(strip_non_ascii(paper.title)) + "\n")
        f.write(str(strip_non_ascii(paper.abstract)) + "\n")
    
    subprocess.call(['java', '-jar', linnaeus_path, '--text', tmp_path+"text.txt", '--out', tmp_path+"results.tsv" ])
    
    with open(tmp_path+"results.tsv", 'r') as f:
        results = [ line.strip('\n').split('\t') for line in f.readlines() ]
    
    return set([ r[0].split(':')[-1] for r in results[1:] ])

# <codecell>

print '.'

# <codecell>

def merge(result):
    """
    result : ET xml root
    """
    
    taxa = [ child.text for taxon in ET.findall('.//Taxon') for child in taxon if child.tag == 'TaxId' ]
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
    def __init__(self, uri, abstract, journal, title):
        self.uri = uri
        self.abstract = abstract
        self.journal = journal
        self.title = title

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

# <codecell>

taxonomy_api = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy"
taxonomy_api_get = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy"

# <codecell>

datapath = "/Users/erickpeirson/Desktop/MO/2014-01-18 RDF"

# <codecell>

p = rdflib.parser.Parser()

# <codecell>

g = rdflib.Graph()

# <codecell>

for filename in os.listdir(datapath):
    if filename.split('.')[-1] == 'rdf':
        g.load(datapath+"/"+filename)

# <codecell>

len(g)

# <codecell>

asdict = {}
aslist = []

# <codecell>

for s,p,o in g:
    try:
        asdict[s][p] = o
    except KeyError:
        asdict[s] = {}
        asdict[s][p] = o

# <codecell>

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

str(strip_non_ascii(aslist[600].abstract))

# <codecell>

'human' in species_set

# <codecell>

len(aslist) # number of abstracts

# <codecell>

vocabulary = set([])
for p in aslist:
    s = str(strip_non_ascii(p.abstract))
    words = set(tokenize(s))
    vocabulary = vocabulary | words

# <codecell>

len(vocabulary)

# <codecell>

vocabulary = remove_integers(list(vocabulary))

# <codecell>

len(vocabulary)

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

# <codecell>

orgpath = '/Users/erickpeirson/Desktop/MO'

# <codecell>

with open(orgpath+"/organisms.txt", "r") as f:
    species = [ tuple(line.strip('\n').split('\t')) for line in f.readlines() ]

# <codecell>

species_set = set([ s for s,r,i in species ])

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

journal_orgs.keys()

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
                print word + ' found'
        if e.tag == 'IdList' and len(e) > 0:
            ident = e[0].text
    if found:
        response = urllib2.urlopen(taxonomy_api_get + "&id=" + ident).read()
        xml = ET.fromstring(response)
        for e in xml[0]:
            if e.tag == 'Rank':
                rank = e.text
                print rank
        species.append( (name, rank, ident) )

# <codecell>

len(species)

# <codecell>

with open(datapath+"/organisms.txt", "w") as f:
    for n,r,i in species:
        f.write(n + "\t" + r + "\t" + i + "\n")

# <codecell>

for i in xrange(len(aslist)):
    aslist[i].abstract = tokenize(str(strip_non_ascii(aslist[i].abstract)))

# <codecell>

counts = {}
for org,rank,ident in species:
    counts[org] = 0
    for p in aslist:
        if org in p.abstract:
            counts[org] += 1

# <codecell>

import matplotlib.pyplot as plt

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

# <codecell>

from vsm.corpus import Corpus
from vsm.corpus.util import corpusbuilders as cb
from vsm.model import tf, tfidf, lsa, ldagibbs, beaglecomposite, beaglecontext, beagleenvironment, beagleorder, ldacgsmulti
from vsm.viewer import ldagibbsviewer, tfidfviewer
from vsm.viewer import beagleviewer
import pickle

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


