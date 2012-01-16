#!/usr/bin/python
import sys, os, random, getopt
from pythonlib import matching

# Copyright 2007, 2008, 2009
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Moritz Gerstung,
# Lukas Geyrhofer,
# Osvaldo Zagordi,
# ETH Zurich

# This file is part of ShoRAH.
# ShoRAH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ShoRAH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ShoRAH.  If not, see <http://www.gnu.org/licenses/>.

'''


usage: mm.py filename [max_haplos]

------------
Input: filename should consist of lines like
0 ACAGACGA
1 CAGACGA
4 AAAAAAAAAA

where the first column is the start position and the second column is the read.
The start positions must be indexed from 0.

Any alphabet is allowed for the reads.  

Duplicate reads or reads which are redundant might cause unexpected results.
For example, the reads
0 ACGT
1 CGT
are redundant (since the first contains the second exactly).

------------
Output:
	the last 4 characters of filename are deleted and replaced with geno
	so if input is file.read then output is file.geno
------------
parameters (see the beginning of main() in the code to change these)

MAXGENO = 10000 (can be changed with max_haplos on command line)
 - maximum number of haplotypes we output

MAXITER = 500
 - maximum number of times we rerun the matching algorithm to try to find new haplotypes

ONEPATH = True
 - extend chains to just one path in the graph?  If false, chains get extended to all possible paths.
   This can be a problem if the chain is small, you can get too many paths.


'''



''' 
TODO
 use getopt to parse options
 don't just strip off last 4 characters of filename to get outputfile
 ONEPATH should have an option to output up to MAXPATHS paths
'''

#import psyco
#psyco.full()
#psyco.log()

def permuteList(l,perm):
	#apply a permutation to a list
	return [perm[i] for i in l]

def permuteGraph(graph,perm):
	#apply a permutation to a graph
	return dict([(perm[g],permuteList(child,perm)) for (g,child) in graph.iteritems()])

def permute(graph, reads, descList):
	''' pick a random elt of S_{len(graph) - 2} and apply the permutation
	    (fixing the source and sink) to graph, reads, and descList
		'''
	#build permutation in S_n
	n = len(graph) - 2
	perm = range(n)
	random.shuffle(perm)
	#and then make sure to fix the source and sink
	perm.append(source)
	perm.append(sink)
	if verbose: print 'permuting with', perm 
	newGraph = permuteGraph(graph,perm)
	#newReads = [reads[i] for i in perm[:-2]] # THIS DOESN'T WORK --  NEED INVERSE OF PERM
	newReads = [[] for i in range(n)]
	for i in range(n):
		newReads[perm[i]] = reads[i]
	#permute labels in descList
	tmpdescList = map(lambda l: permuteList(l, perm), map(list,descList))
	#also permute places ### VIA INVERSE
	newdescList = [[] for i in tmpdescList]
	for i in range(n + 2):
		newdescList[perm[i]] = tmpdescList[i]
	return newGraph, newReads, newdescList

def readFasta(filename):
	try:
		f=open(filename, 'r')
	except IOError:
		return [], 0
	ans = []
	maxEnd = 0
	for line in f:
		pos, seq = line.split()
		pos = int(pos)
		end = pos + len(seq)
		if end > maxEnd: maxEnd = end
		ans.append([pos,seq])
	return ans, maxEnd

def printFasta(seqs, lineLen = 80):
	start = 0
	totalLen = len(seqs[0])

	while start < totalLen:
		length = min(lineLen, totalLen - start)
		print
		for i in seqs:
			print i[start:start + length]
		start += length
	return 0

def addSourceSink(graph, descList,reads, seqLen):
	""" 
	adds a source and sink onto the graph
	and into descList
	"""
	source = len(reads)
	sink = len(reads) + 1
	maxstart = 0
	for r in reads:
		if r[0] > maxstart: maxstart = r[0]
	#reads.append([-1,''])
	#reads.append([seqLen + 1,''])

	### need edges source -> r whenever r has no parents
	### and r -> sink whenever it has no children

	#find everything which is not a child

	allChildren = set([])
	map(allChildren.update, [graph[i] for i in graph])
	initialReads = set(range(len(graph))).difference(allChildren)
	graph.setdefault(source,list(initialReads))

	#find everything without children and create edge to sink
	for r in range(len(graph)):
		if graph[r] == []: 
			graph[r].append(sink)
	graph.setdefault(sink,[])

	#print graph
	if OUTPUTGRAPH:
		print 'BEGIN'
		for r in graph:
			#print 'node name', r
			#print 'node start', reads[r][0]
			#print 'children', graph[r]
			if r == sink:
				print 1000, ':', r, ':', []
			if r == source:
				print -1, ':', r, ':', graph[r]
			if (r != sink and r != source):
				print reads[r][0], ':', r, ':', graph[r]
		print 'END'
		print graph
	#fix up descList
	#everything as sink as desc
	for i in descList:
		i.add(sink)
	#push source and then sink on descList
	#source has everthing as desc
	sourceDesc = set(range(len(descList)))
	sourceDesc.add(sink)
	sinkDesc = set([])
	descList.append(sourceDesc)
	descList.append(sinkDesc)


def makeGraph(reads):
	def readsMatch(r,s):
		# check if two reads agree on their overlap
		# called with r starting before s
		overlap = r[0] + len(r[1]) - s[0]
		if overlap < MINOVERLAP: return 0
		if overlap > len(s[1]): overlap = len(s[1])
		#if overlap < 0: overlap = 0
		return r[1][-overlap:] == s[1][:overlap]

	graph = {}
	# initialize graph
	for r in range(len(reads)):
		graph.setdefault(r,[])

	# set edges between reads
	for r in range(len(reads)):
		if r%100 == 0: print r
		for s in range(len(reads)):
			if reads[r][0] < reads[s][0]:
				if readsMatch(reads[r],reads[s]):
					#edge from r to s
					graph.setdefault(r,[]).append(s)
	return graph

def findDescendants(graph):
	descList = []
	for i in range(len(graph)):
		descList.append(set([]))

	def getDesc(node):
		#print 'getDesc called on', node
		#node is an int,
		#descList is a list of lists of ints
		#set descList[node] by:
		# union of all descList[x] for x in graph[node]
		if len(descList[node]) > 0:
			#print node, descList[node]
			return descList[node]
		else:
			#descList[node] = descList[node].union(set(graph[node]))
			descList[node].update(set(graph[node]))
			for x in graph[node]:
				#doesn't speed it up
				#if len(descList[x]) > 0:
				#	descList[node] = descList[node].union(descList[x])
				#else:
				#descList[node] = descList[node].union(getDesc(x))
				descList[node].update(getDesc(x))
			return descList[node]

	#must call getDesc on the entire graph
	# but most of the calls are trivial
	for i in range(len(graph)):
		getDesc(i)
		#print 'desc of', i, 'are', getDesc(i)
	return descList

def removeCycles(graph, descList):
	for i in graph.keys():
		##print 'checking the',len(graph[i]),'children of node',i
		#python is a pain if we modify graph[i] while looping over it
		#so need a temporary variable 
		tmp = []
		tmp.extend(graph[i])
		for j in tmp:
			for k in tmp:
				if j != k:
					#print 'checking',j,k
					if j in descList[k]:
						#print 'removing node',j
						#print 'graph[i] =',graph[i]
						graph[i].remove(j);
						break;
		###print 'reduced to',len(graph[i]), 'children'
	return graph

def matchToChain(match, n):
	''' takes a matching and creates the 
	corresponding chain in the poset 
	Who knows how it works...'''
	startpoints = set(range(n)).difference(match.values())
	allChains = []
	for i in startpoints:
		curChain = []
		curChain.append(i)
		tmp = i
		while True:
			if tmp in match.keys():
				tmp = match[tmp]
				curChain.append(tmp)
			else:
				break
		allChains.append(curChain)
	return allChains

def onePath(x,y):
	''' find one path between x and y 
		using descList
		currently unused, but might be nice
	'''
	path = [x]
	if y not in descList[x]: return None
	if x == y: return path
	done = 0
	while not done:
		for child in graph[path[-1]]:
			if y in descList[child]:
				path.append(child)
				break
			if y == child:
				path.append(y)
				done = 1
				break
	return path

def find_all_paths(GRAPH, start, end, path=[]):
	# from Guido's essay
	path = path + [start]
	if start == end:
		return [path]
	if not GRAPH.has_key(start):
		return []
	paths = []
	for node in GRAPH[start]:
		if node not in path:
			newpaths = find_all_paths(GRAPH, node, end, path)
			for newpath in newpaths:
				paths.append(newpath)
	return paths

def allPaths(x,y):
	##print 'allPaths',x,y
	# find all paths from x to y in the graph
	paths = [[x]]
	if y not in descList[x]: return None
	if x == y: return paths
	finishedpaths = []
	while len(paths) > 0:
		#for last element of each list in paths
		# loop over children, making new list of paths
		newpaths = []
		#for every path
		for i in range(len(paths)):
			#try to extend it
			for child in graph[paths[i][-1]]:
				#if we can extend along this child
				if y in descList[child]:
					#y is later, so create a new, longer path
#FIXME: test this, possibly just use find_all_paths
					#tmppath = []
					#tmppath.extend(paths[i])
					#tmppath.append(child)
					newpaths.append(paths[i] + [child])
				if y == child:
					#done, push on final
					tmppath = []
					tmppath.extend(paths[i])
					#don't include y
					#tmppath.append(child)
					finishedpaths.append(tmppath)
		#print 'newpaths = ',newpaths
		paths = newpaths
	#print 'returns', len(finishedpaths),'paths'
	return finishedpaths

def pathToGeno(path):
	geno = ''
	'''
	first get rid of sink and source 
	(this has the side effect of modifying the paths)
	'''
	if len(path) == 0: return None
	if path[0] == source:
		path.pop(0)
	if len(path) == 0: return None
	if path[-1] == sink:
		path.pop(-1)
	for i in range(len(path) - 1):
		r = path[i]
		s = path[i+1]
		#nonoverlap is how much of read r occurs before read s starts
		nonoverlap = reads[s][0] - reads[r][0]
		geno = geno + reads[r][1][:nonoverlap]
		#geno.append(reads[r][1][:nonoverlap]
	#add the last read
	geno = geno + reads[path[-1]][1]
	#geno.append(reads[path[-1]][1])

	#pad with .'s at beginning and end to make right length
	countBegin = 0
	countEnd = 0

	for i in range(reads[path[0]][0]):
		#geno.append('.')
		geno = '.' + geno
		countBegin+=1
	
	offsetEnd = seqLen - reads[path[-1]][0] - len(reads[path[-1]][1])
	for i in range(offsetEnd):
		geno = geno + '.'
		countEnd+=1

	'''
	print 'path of length',len(path), 'has', countBegin, '. and', countEnd
	if countBegin > .5 * len(geno):
		geno = None
	if countEnd > .5 * len(geno):
		geno = None
	'''

	return geno



def extendChainToOne(chain):
	#extend chain to ONE path in the graph
	#input: a chain (list of ints)
	#output: a list containing a path that extends the chain
	ans = []
	if len(chain) == 0: return ans
	if chain[0] == source:
		chain.pop(0)
	if len(chain) == 0: return None
	if chain[-1] == sink:
		chain.pop(-1)
	#deal with paths from source to chain[0] but
	# if chain started as [source, sink] then it is of
	# length 0 now
	if len(chain) == 0:
		ans = onePath(source,sink)
	else:
		ans = onePath(source,chain[0])
	for i in range(len(chain)):
		if i != len(chain) - 1: nxt = chain[i+1]
		else: nxt = sink
		if nxt not in graph[chain[i]]:
			ans.extend(onePath(chain[i],nxt))
		else:
			ans.append(chain[i])
	return [ans]


def extendChain(chain):
	#extend chain to all possible paths in graph
	#input: a chain (list of ints)
	#output: a list of all paths that extend the chain
	#print 'input:',chain
	ans = [[]]
	newans = []
	if len(chain) == 0: return ans
	if chain[0] == source:
		chain.pop(0)
	if len(chain) == 0: return None
	if chain[-1] == sink:
		chain.pop(-1)
	#deal with paths from source to chain[0] but
	# if chain started as [source, sink] then it is of
	# length 0 now
	if len(chain) == 0:
		ans = allPaths(source,sink)
	else:
		ans = allPaths(source,chain[0])
	for i in range(len(chain)):
		#print 'currently have', len(ans),'chains'
		if i != len(chain) - 1: nxt = chain[i+1]
		else: nxt = sink
		if nxt not in graph[chain[i]]:
			for p in allPaths(chain[i],nxt):
				for j in ans:
					newans.append(j + p)
			ans = newans
			newans = []
		else:
			for j in ans:
				j.append(chain[i])
	##print 'output:',ans
	return ans


def minimalCover(graph, reads, descList):
	""" the main code: takes a graph and creates a matching
	    problem.  then solves the matching to get chains in 
		the graph, which are extended to paths in all ways """

	#print 'making bipartite graph'
	#make bipartiteGraph
	biGraph = {}
	for i in range(len(graph)):
		biGraph.setdefault(i,list(descList[i]))
	
	print '...finding matching'
	#find matching
	a = matching.bipartiteMatch(biGraph)
	# a[0] is the matching, a[1] and a[2] the parts of
	# the maximal independent sets in the two parts of biGraph
	print '...Reversing matching'
	#reverse the matching
	match = dict([(v,k) for (k,v) in a[0].iteritems()])
	#free some memory, perhaps
	del a
	#print 'matching =', match

	#transform to chains in poset
	chainList = matchToChain(match,len(graph))
	if OUTPUTGRAPH:
		print 'chains=',chainList
	maximalAntiChain = len(chainList)
	print '...largest antichain of size', maximalAntiChain
	if verbose:
		for zz in range(len(chainList)):
			print zz, '==', chainList[zz]

	#now that we're done with the matching, add the source and sink into graph
	#this is needed for the extendChain and pathToGeno functions
	# now done before the matching
	#addSourceSink(graph, descList, reads, seqLen)
	
	print '...transforming chains to paths'
	##transform chains to paths in all possible ways
	if ONEPATH:
		paths = map(extendChainToOne,chainList)
	else:
		paths = map(extendChain,chainList)
	
	print len(chainList), 'chains give', sum(map(len, paths)), 'paths (', map(len,paths), ')'
	if verbose:
		print 'paths=',paths

	print '...translating paths->genotypes'
	genos = []
	for i in range(len(paths)):
		genos.extend(map(pathToGeno,paths[i]))
	
	return genos, maximalAntiChain

def countPaths(graph,source,sink):
	#numPaths[i] = # of paths from i to the sink in graph
	numPaths = [-1 for i in range(len(graph))]
	numPaths[sink] = 1
	def getPaths(node):
		''' count paths from node to sink '''
		if numPaths[node] < 0:
			numPaths[node] = sum(map(getPaths, graph[node]))
		#print 'node', node, 'has', numPaths[node], 'paths to sink'
		return numPaths[node]
	getPaths(source)
	#print 'paths =', numPaths
	return numPaths[source]


###########main code#########


def main(opts, args):
	
	global verbose, seqLen, descList, graph, reads, source, sink
	global ONEPATH
	global MINOVERLAP
	global OUTPUTGRAPH
	
	''' new options parsing 
	usage:
	-h 
	-m maxgeno
	-i maxiter
	-o (boolean, onepath)
	-a (boolean, onepath)
	-l minoverlap
	-g (output graph)
	'''
	
	usage = "usage: %s [-h -g -m MAXGENO -i MAXITER -l MINOVERLAP --one --all ] filename [max_haplos]" % os.path.basename(sys.argv[0])

	#defaults:
	# keep searching until we find at least MAXGENO haplotypes
	MAXGENO = 10000
	# but only look up to MAXITER rounds
	MAXITER = 50 
	# extend chains to ONEPATH or many ?
	ONEPATH = True
	#minimum overlap to include in graph
	MINOVERLAP = 1
	#output a graph?
	OUTPUTGRAPH = False
	for o, a in opts:
		if o in ("-h", "--help"):
			print usage
			sys.exit()
		if o in ("-o", "--one"):
			ONEPATH = True
		if o in ("-m", "--maxhaplo"):
			MAXGENO = int(a)
		if o in ("-i", "--maxiter"):
			MAXITER = int(a)
		if o in ("-a", "--all"):
			ONEPATH = False
		if o in ("-l", "--length"):
			MINOVERLAP = int(a)
		if o in ("-g", "--graph"):
			OUTPUTGRAPH = True
	''' Parse remaining arguments '''
	#print 'arguments:', args
	if len(args) < 1:
		print usage
		sys.exit(2)
	else:
		readfile = args[0]
		if len(args) >= 2:
			MAXGENO = int(args[1])

	''' want to do outfile = readfile; outfile ~= s/rest/geno$/ '''
	if readfile[-4:] == 'rest' or readfile[-4:] == 'read':
	#outFileName = readfile[:-4] + str(MAXGENO) + '.geno'
		outFileName = readfile[:-4] + 'geno'
	else:
		outFileName = readfile + '.geno'



	print 'inputfile =', readfile
	print 'outputfile=', outFileName
	print 'Creating up to', MAXGENO, 'haplotypes'
	print 'Permuting up to', MAXITER, 'times'
	print 'ONEPATH = ', ONEPATH
	print 'minimum overlap in graph:', MINOVERLAP

	''' get reads '''
	reads, seqLen = readFasta(readfile)
	if seqLen == 0:
		print 'No reads in file (or couldn\'t open file)'
		return 0
	source = len(reads)
	sink = len(reads) + 1
	if len(reads) < 51: verbose = 1
	else: verbose = 0
	print len(reads),'reads, sequence length =', seqLen
	if verbose: 
		for r in reads: 
			print reads.index(r), r

	''' Build read graph '''
	graph = makeGraph(reads)
	if verbose: print 'graph = ', graph
	print 'finished graph'

	''' take transitive closure of graph '''
	descList = findDescendants(graph)
	if verbose: print 'desc =', map(list,descList)
	print '...done with transitive closure'

	print '...removing cycles'
	graph = removeCycles(graph,descList)
	if verbose: print 'graph \ cycles = ', graph

	''' add on a source and sink to the graph '''
	addSourceSink(graph, descList, reads, seqLen)
	if verbose:
		print 'graph with SS', graph
		print 'new Desclist', map(list,descList)
	''' count paths in the graph '''
	
	totalPaths = countPaths(graph, source, sink)
	print 'There are', totalPaths, 'paths in the graph'


	#global maximalAntiChain

	allGenos = set([])
	DONE = False


	iter = 0
	if totalPaths < MAXGENO:
		print 'This is less than MAXGENO =', MAXGENO
		print 'So matching is not run'
		paths = extendChain([source,sink])
		genos = map(pathToGeno,paths)
		allGenos.update(genos)
		DONE = True
		genos, maximalAntiChain = minimalCover(graph,reads, descList)
	else:
		print 'starting matching'
		print
		while (not DONE) and (len(allGenos) < MAXGENO) and (iter < MAXITER):
			''' Done with preproccessing. Now permute
		    graph, reads, and descList and get a new set of genotypes
			then find a minimal chain decomposition of the graph
		    and return all genotypes obtained by extending
			these chains '''
			if iter > 1:
				print '...permuting'
				graph, reads, descList = permute(graph,reads,descList)
			genos, maximalAntiChain = minimalCover(graph, reads, descList)
			sizeBefore = len(allGenos)
			allGenos.update(genos)
			sizeAfter = len(allGenos)
			print 'found', len(allGenos), 'genotypes out of', MAXGENO
			print 'gained', sizeAfter - sizeBefore, 'in round', iter
			print
			#printFasta(genos,1000)
			iter+=1


	outfile = open(outFileName, 'w')

	for i in list(allGenos):
		#check if i == None
		if i:
			outfile.write(i)
			outfile.write('\n')
	outfile.close()
	if verbose:
		printFasta(list(allGenos),1000)
	
	#summary
	print 'graph has', totalPaths, 'total paths'
	print 'minimal path cover is of size', maximalAntiChain
	print 'ran matching algorithm', iter, 'times'
	print 'output', len(allGenos), 'haplotypes'

	print '% TOTALPATHS =', totalPaths
	print '% ANTICHAIN =', maximalAntiChain
	print '% HAPLO =', len(allGenos)


################

'''
WARNINGS:
	only works if reads are indexed from 0 (e.g., start positions)
	this is pretty fundamental to the code, unfortunately
'''

if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ghm:i:oal:", ["graph", "help", "maxhaplo=", "maxiter=", "one", "all", "length"])
	except getopt.GetoptError:
		# print help information and exit:
		print usage
		sys.exit(2)
	main(opts, args)

#	import hotshot
#	prof = hotshot.Profile("hotshot_edi_stats")
#	prof.runcall(main)
#	prof.close()

