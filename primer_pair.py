
class primer_pair():

	def __init__(self, scheme, region, output, n):

		#stop using setattr
		setattr(self, 'left_name', '%i_%i_left_%i' %(scheme, region, n))
		setattr(self, 'left_start', output['PRIMER_LEFT_%i' %n][0])
		setattr(self, 'left_end', output['PRIMER_LEFT_%i' %n][0] + output['PRIMER_LEFT_%i' %n][1])
		setattr(self, 'left_length', output['PRIMER_LEFT_%i' %n][1])
		setattr(self, 'left_seq', output['PRIMER_LEFT_%i_SEQUENCE' %n])
		setattr(self, 'left_gc', output['PRIMER_LEFT_%i_GC_PERCENT' %n])
		setattr(self, 'left_tm', output['PRIMER_LEFT_%i_TM' %n])
		
		setattr(self, 'right_name', '%i_%i_right_%i' %(scheme, region, 0))
		setattr(self, 'right_start', output['PRIMER_RIGHT_%i' %n][0])
		setattr(self, 'right_end', output['PRIMER_RIGHT_%i' %n][0] + output['PRIMER_RIGHT_%i' %n][1])
		setattr(self, 'right_length', output['PRIMER_RIGHT_%i' %n][1])
		setattr(self, 'right_seq', output['PRIMER_RIGHT_%i_SEQUENCE' %n])
		setattr(self, 'right_gc', output['PRIMER_RIGHT_%i_GC_PERCENT' %n])
		setattr(self, 'right_tm', output['PRIMER_RIGHT_%i_TM' %n])


class region_object():

	def __init__(self, scheme, region, output):

		self.pairs = []
		for n in range(5):
			self.pairs.append(primer_pair(scheme, region, output, n))

