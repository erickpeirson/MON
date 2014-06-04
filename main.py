# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>



from nltk.corpus import stopwords
import string
import rdflib
import urllib2
import xml.etree.ElementTree as ET
import os
import subprocess
import pickle
import matplotlib.pyplot as plt
import tempfile
from collections import Counter
import csv

def tokenize(s):
    s = s.translate(string.maketrans("",""), string.punctuation)
    return s.lower().split()

def remove_stopwords(t):
    """
    Removes words found in the NLTK ``stopwords`` corpus.
    
    Parameters
    ----------
    t : list
        A list of strings.
    
    Returns
    -------
    list
        A somewhat shorter list of strings.
    """

    return [ word for word in t if word not in stopwords.words() and word != 'ab' ]



def strip_non_ascii(string):
    """
    Beligerently removes non-ASCII characters
    
    Parameters
    ----------
    string : str
    
    Returns
    -------
    stripped_string : str
    """
    
    stripped = (c for c in string if 0 < ord(c) < 127)
    stripped_string = ''.join(stripped)
    return stripped_string

def remove_integers(t):
    """
    Removes words that cannot be cast as integers.
    
    Parameters
    ----------
    t : list
    
    Returns
    -------
    new_t : list
    """

    new_t = []
    for word in t:
        try:
            int(word)
        except ValueError:
            new_t.append(word)
    return new_t



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

class paper:
    """
    Describes a minimal PubMed entry.
    """
    def __init__(self, uri, abstract, journal, title):
        self.uri = uri
        self.abstract = abstract
        self.journal = journal
        self.title = title

def to_path(uri):
    """
    Converts a URI to something that can be used as a filepath.
    
    Parameters
    ----------
    uri : str
    
    Returns
    -------
    path : str
    """
    path =  str(uri).replace('/', '___')
    return path



def to_uri(path):
    return str(path).replace('___', '/')


def merge(result):
    """
    Uses NCBI Taxonomy XML response to decide to which higher-level
    taxonomic group a taxon belongs.

    Parameters
    ----------
    result : ET XML object
        NCBI Taxonomy response.
    
    Returns
    -------
    str
        Taxon category (e.g. 'Human', 'Mouse or Rat', etc).
    """
    
    # All higher-level taxonomic groups to which this taxon belongs.
    taxa = [ tax.text for tax in result.findall('.//TaxId') ]

    if '33208' in taxa and '7742' not in taxa:    # Any metazoan that is not a vertebrate.
        return 'Invertebrate'
    if '89593' in taxa: # Craniata
        if '7898' in taxa or '7894' in taxa or '7878' in taxa: # Fishes
            return 'Fish'
        elif '8782' in taxa:    # Aves
            return 'Bird'
        elif '1294634' in taxa or '8505' in taxa or '8509' in taxa or '8459' in taxa: # Crocodylia (in Amniota); Sphenodontia; Squamata; Testudines
            return 'Reptile'
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
        return 'Other Craniate'
    else: # Not in Craniata
        return 'Non-Craniate'



def classify(taxon):
    """
    Yields the higher-level taxonomic group to which taxon belongs.

    Retrieves data for taxon from NCBI Taxonomy, then passes the result
    to merge() to classify. Uses classified (dict) to cache results.

    Parameters
    ----------
    taxon : str
        NCBI Taxonomy ID.

    Returns
    -------
    merged : str
        Taxon category (e.g. 'Human', 'Mouse or Rat', etc).
    """
    try:
        return classified[taxon]
    except KeyError:
        response = urllib2.urlopen(taxonomy_api_get + "&id=" + taxon).read()
        xml = ET.fromstring(response)
        classified[taxon] = merge(xml)
        return merged



def find_names_batch(papers, linnaeus_path="/Applications/linnaeus/bin/linnaeus-2.0.jar", outpath='./'):
    """
    Named Entity Recognition step.

    Generates a temporary directory containing one text file per paper (title+abstract),
    and throws LINNAEUS at it.

    Parameters
    ----------
    papers : list
        A list of :class:`paper` objects.
    linnaeus_path : str
        Path to LINNAEUS .jar file.
    outpath : str
        Target directory for NER results.

    Returns
    -------
    resultspath : str
        Path to CSV-formated NER results file.
    """
    temp = tempfile.mkdtemp()
    for paper in papers:
        filepath = '{0}/{1}.txt'.format(temp, to_path(paper.uri))
        with open(filepath, 'w') as f:
            f.write(str(strip_non_ascii(paper.title)) + "\n")
            f.write(str(strip_non_ascii(paper.abstract)) + "\n")

    resultspath = '{0}/orgnames.csv'.format(outpath)    # CSV results will go here.
    
    # Call LINAEUS.
    subprocess.call(['java', '-jar', linnaeus_path, '--textDir', temp, '--out', resultspath ])

    del temp    # Cleanup.
    
    return resultspath



def handle_result_taxon(result):
    """
    Parse the taxon ID porition of LINNAEUS NER result.

    If multiple possible taxa are suggested, chooses the one with the greater
    confidence level.

    Parameters
    ----------
    result : list
        A CSV-parsed row from the LINNAEUS NER results file.

    Returns
    -------
    str
        NCBI Taxonomy ID for NER-matched taxon.
    """
    taxa = [ r.split(':')[-1] for r in result[0].split('|') ]
    if len(taxa) > 1:   
        confidence = []
        ids = []        
        for t in taxa:
            t_c = t.split('?')
            confidence.append(t_c[1])
            ids.append(t_c[0])
    
        return ids[confidence.index(max(confidence))]
    else:
        return taxa[0].split('?')[0]



def handle_result(result):
    """
    Assigns taxon matches from NER to each respective paper.

    Parameters
    ----------
    result : list
        A CSV-parsed row from the LINNAEUS NER results file.        
    """    
    taxon = handle_result_taxon(result)
    uri = rdflib.term.URIRef(to_uri(result[1]))
    orgs_in_paper[uri].append(taxon)



def handle_paper_taxa(uri):
    """
    Applies logic for assigning taxon-group hits to journals.

    Humans are treated slightly differently than the rest. To accrue a hit
    for 'Human', it must be the only hit for that paper, or there must have
    been at least two hits for that paper.

    Parameters
    ----------
    uri : str
        URI belonging to a :class:`.paper`.

    Returns
    -------
    to_increment : dict
        Keys are taxon-groups to increment for corresponding journal,
        values are the amount to increment each taxon-group.
    """
    to_increment = {}
    taxa = taxa_in_paper[uri]
    taxon_hits = taxa_hits_in_paper[uri]
    
    for taxon, count in taxa.iteritems():
        if taxon == 'Human':
            if len(taxa) == 1 or taxon_hits['Human'] > 1:
                to_increment['Human'] = 1
        else:
            to_increment[taxon] = count
    return to_increment
    
if __name__ == '__main__':


    orgpath = '/Users/erickpeirson/Model-Organisms-in-Neuroscience/data'



    datapath = '/Users/erickpeirson/Model-Organisms-in-Neuroscience/data/RDF'



    taxonomy_api = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy"
    taxonomy_api_get = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy"


    # ## Load Data



    g = rdflib.Graph()
    for filename in os.listdir(datapath):
        if filename.split('.')[-1] == 'rdf':
            g.load(datapath+"/"+filename)


    # ### `asdict`: Index by subject



    asdict = {}    # Index by subject.
    for s,p,o in g:  # Subject, predicate, and object in RDF triple.
        try:
            asdict[s][p] = o
        except KeyError:
            asdict[s] = {}
            asdict[s][p] = o


    # ### `aslist`: Generate a list of `paper`s



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



    len(aslist) # number of abstracts


    # ## LINNAEUS NER & classification via NCBI Taxonomy
    # 
    # ### Named Entity Recognition
    # * [LINNAEUS](http://linnaeus.sourceforge.net/)
    # 
    # ### Classification logic:
    # * For each paper in a journal:
    #     * For each unique taxon in that paper,
    #         * If the taxon group is human:
    #             * If human is the only taxon group in that paper, or if the number of hits for human is greater than 1:
    #                 * Increment the human group
    #         * Otherwise,
    #             * Increment the taxon group


    # ### Globals



    classified = {}    # For caching classification results.



    orgs_in_paper = { p.uri:[] for p in aslist }



    taxa_in_paper = { p.uri:Counter() for p in aslist }
    taxa_hits_in_paper = { p.uri:Counter() for p in aslist }



    taxa_in_journal = {}


    # ### Methods






    # ### Execution



    resultpath = find_names_batch(aslist)    # <-- The action.



    # Read the CSV-formatted NER results file from LINNAEUS.
    results = []
    with open(resultpath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            results.append(row)



    # Assign taxon matches to papers.
    for result in results[1:]:
        handle_result(result)



    # Classify taxon matches, and accrue to papers.
    for uri,ts in orgs_in_paper.iteritems():
        for t in ts:    # Count number of taxon-group hits in title and abstract.
            c = classify(t)
            taxa_hits_in_paper[uri][c] += 1
        for t in set(ts):    # Count number of unique taxon-group hits.
            c = classify(t)
            taxa_in_paper[uri][c] += 1



    # Index journal names by the URIs of their papers.
    journal = {}
    for pap in aslist:
        # Get object of dc.title triple for paper's journal (pap.journal is a journal URI).
        title = str(asdict[pap.journal][rdflib.term.URIRef(u'http://purl.org/dc/elements/1.1/title')])
        journal[pap.uri] = title.lower()



    # Assign taxon-group hits to journals.
    for pap in aslist:
        ti = handle_paper_taxa(pap.uri)
        jtitle = journal[pap.uri]
        for k,v in ti.iteritems():
            try:
                taxa_in_journal[jtitle][k] += v
            except KeyError:
                taxa_in_journal[jtitle] = Counter()
                taxa_in_journal[jtitle][k] += v


    # ## Generate figures



    tkeys = ['Human', 
             'Macaque', 
             'Other Primate',
             'Mouse or Rat',
             'Cat',
             'Other Mammal',
             'Bird',
             'Reptile',
             'Amphibean',
             'Fish',
             'Invertebrate']



    plt.figure(figsize=(8,50))
    J = len(taxa_in_journal)
    sums = Counter()
    for j in xrange(J):
        key = taxa_in_journal.keys()[j]
        counts = []
        for t in tkeys:    # Use the same keys for all figures, so that taxon-groups are aligned.
            try:
                c = taxa_in_journal[key][t]
            except KeyError:
                c = 0
            counts.append(c)
            sums[t] += c

        ax = plt.subplot(J, 1, j+1)    # Subplots are arranged vertically.
        T = len(tkeys)

        ax.bar(arange(T), np.array(counts)+0.1, align='center', color='green', alpha=0.8, lw=2)
        
        # TODO: wrap tick label text?
        plt.xticks(arange(T))
        ax.set_xticklabels(tkeys, rotation=35, horizontalalignment='right')

        # TODO: clean this up....
        if key == 'the neuroscientist: a review journal bringing neurobiology, neurology and psychiatry':
            title = 'The Neuroscientist'
        elif key == 'cerebral cortex (new york, n.y.: 1991)':
            title = 'Cerebral Cortex'
        elif key == 'acs chemical neuroscience':
            title = 'ACS Chemical Neuroscience'
        else:
            title = key.title()
            
        plt.title('{0}, 2005-2010'.format(title))
    plt.tight_layout()    # Fixes subplot spacing.
    plt.show()



    plt.figure( figsize=(10,7))
    plt.bar(arange(len(sums)), [ sums[t] for t in tkeys ], align='center', color='green', alpha=0.8, lw=2)
    plt.xticks(arange(len(sums)), tkeys, rotation=35, horizontalalignment='right')
    plt.title('All Journals, 2005-2010')
    plt.show()

