# This file takes an input directory (IN_DIR) with output file from 
# Bismarck (one for each cell) and produces intermediate files in an 
# output directory (OUT_DIR), which can then be consolidated into a
# HDF5 file by the "part2"  script.

const IN_DIR = "eb/"  # do not omit trailing slash
const OUT_DIR = "bstmp/"   # ditto; must exist and should be empty

import GZip
import CSV

const eof_marker = "zzzz"
# This eof_marker should sort later than all chromosome names.

struct GenomicPosition
	chrom ::String
	pos ::UInt32
end

function Base.isless( gp1 ::GenomicPosition, gp2 ::GenomicPosition )
	if gp1.chrom == gp2.chrom
		return gp1.pos < gp2.pos
	else
		return gp1.chrom < gp2.chrom
	end
end

function Base.show( io ::IO, gp ::GenomicPosition )
	print( io, gp.chrom, ":", gp.pos )
end

struct BismarckRecord
	gp ::GenomicPosition
	count_unmeth ::UInt8
	count_meth ::UInt8
end

struct BismarckInput
	io ::IO
	_prev_chrom ::Base.RefValue{String}
	function BismarckInput( io ::IO )
		firstline = readline( io )
		if firstline != "chr\tpos\tmet_reads\tnonmet_reads\trate\n" 
			error( "Wrong file formar" )
		end
		return new( io, Ref("") )
	end
end

function next( b ::BismarckInput )
	if eof(b.io)
		return BismarckRecord( GenomicPosition(eof_marker, 0), 0, 0 ) 
	end
	line = readline( b.io )
	if line == "" || line == "\n"
		return BismarckRecord( GenomicPosition(eof_marker, 0), 0, 0 ) 
	end	
	fields = split( line, "\t" )
	if length(fields) != 5
		println( ">",line,"<")
		error( "Input line does not contain 5 fields." )
	end
	gp = GenomicPosition( fields[1], parse( Int, fields[2] ) )
	br = BismarckRecord( gp, parse( UInt8, fields[4] ), parse( UInt8, fields[3] ) )
	if gp.chrom != b._prev_chrom[]
		# New chrosome starts, do a few checks
		if gp.chrom > eof_marker
			error( "Chromosome name sorts after EOF marker." )
		end
		if gp.chrom < b._prev_chrom[]
			error( "Chromosomes are not lexically sorted in input file." )
		end
		b._prev_chrom[] = gp.chrom
	end
	return br
end

function main()
	# Open input files
	input_filenames = readdir( IN_DIR )
	input_files = [ GZip.open( IN_DIR*fn ) for fn in input_filenames ]
	inputs = [ BismarckInput(f) for f in input_files ]

	# Write cell list
	cell_list = open( OUT_DIR*"cells.lst", "w" )
	for fn in input_filenames
		println( cell_list, fn )
	end
	close( cell_list )

	# Open other output files
	values_out = open( OUT_DIR*"values.bin", "w" )
	index_out = nothing  # will be opened later
	cur_idx::UInt32 = 0

	chrom_list_out = open( OUT_DIR*"chromosomes.csv", "w" )

	currec = [ next(bi) for bi in inputs ]
	prev_chrom = ""

	while true
		curpos = minimum( rc.gp for rc in currec )
		if curpos.chrom == eof_marker
			break
		end	
		if curpos.chrom != prev_chrom
			# We start a new chromosome
			if index_out != nothing
				close( index_out )
				println( chrom_list_out, prev_chrom, ",", cur_idx )
				println( "Finished processing chromosome ", prev_chrom)
			end
			prev_chrom = curpos.chrom
			index_out = open( OUT_DIR*curpos.chrom*".idx", "w" )
		end
		#println( curpos )
		write( index_out, UInt32(curpos.pos), UInt32(cur_idx) )
		for i in 1:length(inputs)
			if currec[i].gp == curpos
				#println( "   cell ", i, ": ", currec[i].count_unmeth, "+", currec[i].count_meth )
				currec[i] = next( inputs[i] )
				write( values_out, UInt16(i), UInt8(currec[i].count_unmeth), UInt8(currec[i].count_meth) )
				cur_idx += 1
			end
		end
	end

	for f in input_files
		close(f)
	end

	close( index_out )
	close( values_out )
	
	println( chrom_list_out, prev_chrom, ",", cur_idx )
	println( "Finished processing chromosome ", prev_chrom )
	close( chrom_list_out )


end

main()


# Description of output format:
#
# values.bin
#    each record has 4 Bytes, as follows:
#      Bytes 0,1: cell number, zero-based, UInt16
#      Byte 2: count unmethylated reads
#      Byte 3: count methylated reads
#
# 