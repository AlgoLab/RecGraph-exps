def main():
	bam = pysam.AlignmentFile(sys.argv[1])
	latest_start = -1
	earliest_end = float("inf")
	for read in bam.fetch():
		if read.flag in [0, 16]:
			if (x:= read.reference_start) > latest_start:
				latest_start = x #0-based inclused
			if (x:= read.reference_end) < earliest_end:
				earliest_end = x #0-based excluded

	# print(latest_start+1, earliest_end)

	bam = pysam.AlignmentFile(sys.argv[1])
	r_id=1
	for read in bam.fetch():
		if read.flag in [0, 16]:
			p = read.get_aligned_pairs(matches_only=True)
			ix = 0
			while p[ix][1] < latest_start:
				ix+=1
			rs = p[ix][0]
			ix = len(p)-1
			while p[ix][1] >= earliest_end:
				ix-=1
			re = p[ix][0]		
			# print(f'>{read.query_name}.{rs+1}-{re+1}')
			print(f'>S_{r_id}')
			r_id+=1
			print(str(read.query_sequence[rs:re+1]))

if __name__ == '__main__':
	import pysam
	import sys
	main()