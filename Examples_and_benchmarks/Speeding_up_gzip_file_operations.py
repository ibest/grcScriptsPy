# Demonstrate and benchmark options for reading from and writing to a gzip file using
# pipes.

from subprocess import Popen, PIPE, STDOUT
import time
import gzip
import sys


def sp_gzip_read(file, bufsize=-1):
    p = Popen('gzip --decompress --to-stdout'.split() + [file], stdout=PIPE, stderr=STDOUT, bufsize=bufsize)
    return p.stdout


def sp_gzip_write(file, bufsize=-1):
    filep = open(file, 'wb')
    p = Popen('gzip', stdin=PIPE, stdout=filep, shell=True, bufsize=bufsize)
    return p.stdin


class fastqIter:
    " A simple file iterator that returns 4 lines for fast fastq iteration. "
    def __init__(self, handle):
        self.inf = handle

    def __iter__(self):
        return self

    def next(self):
        lines = {'id': self.inf.readline().strip()[1:],
                 'seq': self.inf.readline().strip(),
                 '+': self.inf.readline().strip(),
                 'qual': self.inf.readline().strip()}
        assert(len(lines['seq']) == len(lines['qual']))
        if lines['id'] == '' or lines['seq'] == '' or lines['+'] == '' or lines['qual'] == '':
            raise StopIteration
        else:
            return lines

    @staticmethod
    def parse(handle):
        return fastqIter(handle)

    def close(self):
        self.inf.close()


def writeFastq(handle, fq):
    handle.write('>' + fq['id'] + '\n')
    handle.write(fq['seq'] + '\n')
    handle.write(fq['+'] + '\n')
    handle.write(fq['qual'] + '\n')


#Test reading in from file
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    filename = 'L5I3_CAGATC_L001_R1_001.fastq.gz'

############### fastqIter+gzip_pipe ##############
t = time.time()
i = 0
for record in fastqIter(sp_gzip_read(filename, bufsize=-1)):
    i += 1
print "Reading" + "-"*100
print "Reading: fastqIter+gzip_pipe bufsize=-1:"
print "\ttotal records: %s" % i
print "\trecords per second: %s" % (i/(time.time() - t))
print "-"*100

############### fastqIter+gzip_pipe unbuffered ##############
t = time.time()
i = 0
for record in fastqIter(sp_gzip_read(filename, bufsize=0)):
    i += 1
print "Reading" + "-"*100
print "Reading: fastqIter+gzip_pipe bufsize=0 (unbuffered):"
print "\ttotal records: %s" % i
print "\trecords per second: %s" % (i/(time.time() - t))
print "-"*100


############### fastqIter + gzip.open ##############
t = time.time()
i = 0
for record in fastqIter(gzip.open(filename, 'rb')):
    i += 1

print "Reading" + "-"*100
print "fastqIter + Python gzip module:"
print "\ttotal records: %s" % i
print "\trecords per second: %s" % (i/(time.time() - t))
print "-"*100


############### buffered gzip_pipe with no parsing  ##############
t = time.time()
i = 0
for record in sp_gzip_read(filename):
    i += 1

print "Reading" + "-" * 100
print "Buffered gzip_pipe with no parsing  "
print "\ttotal lines: %s, total records: %s " % (i, i/4)
print "\trecords per second: %s " % ((i/4)/(time.time() - t))
print "-"*100


###------------------------------------------------------------
#         Writing tests                    #
##################################################
print "\n\n\n\n"

outfname = "test.fastq.gz"
############### buffered gzip_pipe with parsing  ##############
t = time.time()
i = 0
with sp_gzip_write(outfname, bufsize=-1) as oh:
    for record in fastqIter(sp_gzip_read(filename, bufsize=-1)):
        i += 1
        writeFastq(oh, record)

    print "Writing" + "-" * 100
    print "Buffered input and output gzip_pipe"
    print "\ttotal records: %s" % i
    print "\trecords per second: %s" % (i/(time.time() - t))
    print "-"*100


############### buffered gzip_pipe input, unbuffered output with parsing  ##############
t = time.time()
i = 0

with sp_gzip_write(outfname, bufsize=0) as oh:
    for record in fastqIter(sp_gzip_read(filename, bufsize=-1)):
        i += 1
        writeFastq(oh, record)

    print "Writing" + "-" * 100
    print "Gzip_pipe, buffered input and unbuffered output"
    print "\ttotal records: %s" % (i)
    print "\trecords per second: %s " % ((i)/(time.time() - t))
    print "-"*100

############### buffered gzip_pipe input with python_gzip output and parsing  ##############
t = time.time()
i = 0

with gzip.open(outfname, 'wb') as oh:
    for record in fastqIter(sp_gzip_read(filename, bufsize=-1)):
        i += 1
        writeFastq(oh, record)

    print "Writing" + "-" * 100
    print "Buffered input and unbuffered output gzip_pipe"
    print "\ttotal records: %s" % (i)
    print "\trecords per second: %s " % ((i)/(time.time() - t))
    print "-"*100

