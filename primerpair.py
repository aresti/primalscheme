
class primerpair():

	def __init__(self, scheme, region, output):
		
		self.left_0_name = '%i_%i_left_%s' %(scheme, region, 0)
		self.left_0_start = output['PRIMER_LEFT_0'][0]
		self.left_0_end = output['PRIMER_LEFT_0'][0] + output['PRIMER_LEFT_0'][1] 
		self.left_0_length = output['PRIMER_LEFT_0'][1]
		self.left_0_seq = output['PRIMER_LEFT_0_SEQUENCE']
		self.left_0_gc = output['PRIMER_LEFT_0_GC_PERCENT']
		self.left_0_tm = output['PRIMER_LEFT_0_TM']
		
		self.right_0_name = '%i_%i_right_%s' %(scheme, region, 0)
		self.right_0_start = output['PRIMER_RIGHT_0'][0]
		self.right_0_end = output['PRIMER_RIGHT_0'][0] + output['PRIMER_RIGHT_0'][1]
		self.right_0_length = output['PRIMER_RIGHT_0'][1]
		self.right_0_seq = output['PRIMER_RIGHT_0_SEQUENCE']
		self.right_0_gc = output['PRIMER_RIGHT_0_GC_PERCENT']
		self.right_0_tm = output['PRIMER_RIGHT_0_TM']

		self.left_1_name = '%i_%i_left_%s' %(scheme, region, 1)
		self.left_1_start = output['PRIMER_LEFT_1'][0]
		self.left_1_end = output['PRIMER_LEFT_1'][0] + output['PRIMER_LEFT_1'][1] 
		self.left_1_length = output['PRIMER_LEFT_1'][1]
		self.left_1_seq = output['PRIMER_LEFT_1_SEQUENCE']
		self.left_1_gc = output['PRIMER_LEFT_1_GC_PERCENT']
		self.left_1_tm = output['PRIMER_LEFT_1_TM']
		
		self.right_1_name = '%i_%i_right_%s' %(scheme, region, 1)
		self.right_1_start = output['PRIMER_RIGHT_1'][0]
		self.right_1_end = output['PRIMER_RIGHT_1'][0] + output['PRIMER_RIGHT_1'][1]
		self.right_1_length = output['PRIMER_RIGHT_1'][1]
		self.right_1_seq = output['PRIMER_RIGHT_1_SEQUENCE']
		self.right_1_gc = output['PRIMER_RIGHT_1_GC_PERCENT']
		self.right_1_tm = output['PRIMER_RIGHT_1_TM']		
