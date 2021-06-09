import csv
import os
import gzip


class GenomicPosition:

    def __init__(self, chrom, pos):
        self.chrom = chrom
        self.pos = pos

    def __repr__(self):
        return self.chrom + "," + str(self.pos)

    def __lt__(self, other):
        if self.chrom != other.chrom:
            return self.chrom < other.chrom
        else:
            return self.pos < other.pos

    def __eq__(self, other):
        return self.chrom == other.chrom and self.pos == other.pos

    def __ne__(self, other):
        return self.chrom != other.chrom or self.pos != other.pos

def make_pos(d):
    d["gp"] = GenomicPosition( d["chr"], d["pos"] )
    del d["chr"], d["pos"]
    return d

def main():
    path="/Users/marti/UNI/bismark_5_cells/eb/"
    files = [ gzip.open(path+f, mode = 'rt') for f in os.listdir(path) if f.endswith(".tsv.gz")]
    print()
    readers = [ csv.DictReader(f, delimiter='\t') for f in files ]

    currec = [ make_pos(next(r)) for r in readers ]

    writerv = open("values.txt", "w")
    writerp = open("positions.txt", "w")
    line_num = 0

    while True:
        curpos = min( rc["gp"] for rc in currec )
        if curpos.chrom == "ZZZ":
            break
        writerp.write(str(curpos) + "," + str(line_num) + "\n" )

        for i in range(len(readers)):
            if currec[i]["gp"] == curpos:
                writerv.write(str(i) + "," + str(currec[i]["met_reads"]) + "," +
                      str(currec[i]["nonmet_reads"]) + "\n")
                line_num += 1
                try:
                    currec[i] = make_pos( next( readers[i] ) )
                except StopIteration:
                    currec[i] = { "gp": GenomicPosition( "ZZZ", 0 ) }
        print()

    writerv.close()
    writerp.close()
    for f in files:
        f.close()

if __name__ == '__main__':
    main()
