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

def read_bismarck_files(files):
    readers = [ csv.DictReader(f, delimiter='\t') for f in files ]
    currec = [ make_pos(next(r)) for r in readers ]

    line_num = 0

    while True:
        curpos = min( rc["gp"] for rc in currec )
        if curpos.chrom == "ZZZ":
            break

        data = []
        for i in range(len(readers)):
            if currec[i]["gp"] == curpos:
                data.append( { "cell": i, "count_unmeth": currec[i]["nonmet_reads"], "count_meth": currec[i]["met_reads"] } )
                try:
                    currec[i] = make_pos( next( readers[i] ) )
                except StopIteration:
                    currec[i] = { "gp": GenomicPosition( "ZZZ", 0 ) }
        yield curpos, data

def write_text_files( pos_file, val_file, data ):
    line_num = 0
    for pos, data in data:
        pos_file.write(str(curpos) + "," + str(line_num) + "\n" )

        for a in data:
            val_file.write( str(a["cell"]) + "," + a["count_unmeth"] + "," + a["count_meth"] + "\n")
            line_num += 1

def main():
    path="eb"
    files = [ gzip.open(path+"/"+f, mode = 'rt') for f in os.listdir(path) if f.endswith(".tsv.gz") ]
    writerv = open("values.txt", "w")
    writerp = open("positions.txt", "w")

    write_text_files( writerp, writerv, read_bismarck_files(files) )

    writerv.close()
    writerp.close()
    for f in files:
        f.close()

if __name__ == '__main__':
    main()
