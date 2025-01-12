import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='A program to convert a deduplicated bam file to a matrix of sites by cells for scATAC-seq analysis.')
parser.add_argument('-B','--inbam', help='Deduplicated bam file',dest='inbam',required=True)
parser.add_argument('-I','--indextable', help='File listing acceptable barcodes',dest='indextable',required=True)
parser.add_argument('-O','--outdir', help='Output directory',dest='outdir',required=True)
parser.add_argument('-P','--prefix',help='Output file prefix',dest='prefix',required=True)
parser.add_argument('-C','--cutoff',help='Read depth cutoff for cells',dest='cutoff',required=False,default="auto")
parser.add_argument('-W','--windows',help='Bed file containing locations of windows reads should be counted in',dest='windows',required=True)
args = parser.parse_args()

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

if args.outdir[-1] != '/':
    args.outdir = args.outdir + '/'

try:
    os.makedirs(args.outdir)
except OSError:
    print('Outdir already exists...')

print("Counting reads assigned to barcodes...")
scriptdir = os.path.dirname(os.path.abspath(__file__))
# readcounter = 'python ' + scriptdir + '/sc_atac_barcode_read_counter.py ' + args.inbam + ' ' + args.indextable + ' ' + args.outdir + args.prefix + '.readcount.report.txt'
read_rpt = os.path.join(args.outdir, args.prefix + '.readcount.report.txt')
readcounter = 'python {} {} {} {}'.format(
	os.path.join(scriptdir, 'sc_atac_barcode_read_counter.py'),
	args.inbam,
	args.indextable,
	read_rpt)
submitter(readcounter)

print("Determining read depth cutoff...")
# cutoffer = 'Rscript ' + scriptdir + '/sc_atac_cell_cutoff.R ' + args.outdir + args.prefix + ' ' + args.cutoff
cutoffer = '/usr/bin/Rscript {} {} {}'.format(
	os.path.join(scriptdir, 'sc_atac_cell_cutoff.R'),
	read_rpt,
	args.cutoff)
submitter(cutoffer)

print("Building cell matrix...")
# builder = 'python ' + scriptdir + '/sc_atac_window_counter.py ' + args.inbam + ' ' + args.outdir + args.prefix + '.readdepth.cells.indextable.txt ' + args.windows + ' ' + args.outdir + args.prefix + '.cellmatrix.txt True'
builder = 'python {} {} {} {} {} {}'.format(
	os.path.join(scriptdir, 'sc_atac_window_counter.py'),
	args.inbam,
	os.path.join(scriptdir, '2to4.readdepth.cells.indextable.txt'),
	args.windows,
	os.path.join(args.outdir, args.prefix + '.cellmatrix.txt'),
	'True')
# submitter(builder)
print(builder)
