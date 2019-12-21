import sys
import pysam

if len(sys.argv) != 4:
    sys.exit('Usage: python sc_atac_samespecies_individual_readcounter.py [Input Bam file] [Input Index table (NA if no table)] [Out file]')

inputbam = sys.argv[1]
indextable = sys.argv[2]
outfile = sys.argv[3]

totalct = {}
descriptor = {}
descdic = {}
if indextable != 'NA':
    descer = open(indextable,'r')
    for line in descer:
        liner = line.strip().split()
        descdic[liner[0]] = liner[1]

print("Counting total reads...")
bamfile = pysam.Samfile(inputbam,'rb')
for read in  bamfile.fetch():
    tagger = read.qname.split(':')[0]
    try:
        totalct[tagger] += 1
    except KeyError:
        totalct[tagger] = 1
        try:
            descriptor[tagger] = descdic[tagger]
        except KeyError:
            descriptor[tagger] = 'bkgd'

bamfile.close()

# outter = open(outfile,'w')
# print >> outter, 'Barcode\tExperiment\tReadCount'
with open(outfile, 'w') as w:
    w.write('Barcode\tExperiment\tReadCount\n')
    for tag in sorted(totalct.keys()):
        # print >> outter, tag + '\t' + descriptor[tag] + '\t' + str(totalct[tag])
        w.write(tag + '\t' + descriptor[tag] + '\t' + str(totalct[tag]) + '\n')

# outter.close()
