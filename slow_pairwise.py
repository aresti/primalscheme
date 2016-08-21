def slow_pairwise(record, pair):
	print
	print 'record           ', record.id

	left = Seq.Seq(pair.left_0_seq)
	aln = pairwise2.align.localms(left, record, 2, -1, -1, -1)[0]
	left_start = aln[3]
	left_end = aln[4]
	left_score = aln[2]

	right = Seq.Seq(pair.right_0_seq)
	aln = pairwise2.align.localms(right, record, 2, -1, -1, -1)[0]
	right_start = aln[3]
	right_end = aln[4]
	right_score = aln[2]

	print 'primer           ', left
	print 'start            ', left_start
	print 'end              ', left_end
	print 'score            ', left_score	
	print 'primer           ', right
	print 'start            ', right_start
	print 'end              ', right_end
	print 'score            ', right_score

	return (left_start, left_end, left_score, right_start, right_end, right_end)
