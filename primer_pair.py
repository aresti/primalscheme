
class PrimerPair():

	def __init__(self, left, right):

		self.left = left
		self.right = right
		self.total = 0

class Primer():
	
	def __init__(self, scheme, region, output, n, direction):

		self.direction = direction
		self.name = '%i_%i_left_%i' %(scheme, region, n)
		self.start = output['PRIMER_%s_%i' %(direction, n)][0]
		self.end = output['PRIMER_%s_%i' %(direction, n)][0] + (
			output['PRIMER_%s_%i' %(direction, n)][1])
		self.length = output['PRIMER_%s_%i' %(direction, n)][1]
		self.seq = output['PRIMER_%s_%i_SEQUENCE' %(direction, n)]
		self.gc = output['PRIMER_%s_%i_GC_PERCENT' %(direction, n)]
		self.tm = output['PRIMER_%s_%i_TM' %(direction, n)]
		self.sub_total = 0

class Region():

	def __init__(self, scheme, region, output):

		self.scheme = scheme
		self.region = region
		self.pool = '2' if self.region%2==0 else '1'
		self.pairs = []
		for n in range(5):
			left = Primer(scheme, region, output, n, 'LEFT')
			right = Primer(scheme, region, output, n, 'RIGHT')
			self.pairs.append(PrimerPair(left, right))

	def sort_pairs(self):
		
		self.pairs.sort(key=lambda x: x.total, reverse=True)
