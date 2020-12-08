import argparse
import fileinput


parser = argparse.ArgumentParser(
    description='Update GFA from minigraph before vg conversion')
parser.add_argument('-r', default='hg38_chr20',
                    help='name of reference contig')

args = parser.parse_args()

print('H\tVN:Z:1.0')

ref_segs = []

for line in fileinput.input(files='-'):
    line = line.rstrip().split('\t')
    if line[0] == 'S':  # segment
        # change segment names to numeric
        line[1] = line[1][1:]
        # uppercase sequences
        line[2] = line[2].upper()
        # save segment if on reference path
        if line[4] == 'SN:Z:' + args.r:
            ref_segs.append(line[1])
    elif line[0] == 'L':
        # change segment names to numeric
        line[1] = line[1][1:]
        # change segment names to numeric
        line[3] = line[3][1:]
    print('\t'.join(line))

# write path
ref_path = ['P', args.r]
ref_path_segs = []
ref_path_ols = []
for seg in ref_segs:
    ref_path_segs.append(seg + '+')
    ref_path_ols.append('*')
ref_path.append(','.join(ref_path_segs))
ref_path.append(','.join(ref_path_ols[1:]))
# ref_path.append('*')
print('\t'.join(ref_path))
