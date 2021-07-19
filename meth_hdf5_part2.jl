const IN_DIR = "bstmp/"  # do not omit trailing slash
const OUT_FILE = "foo.hdf5"  

using HDF5

function main()

	hf = h5open( OUT_FILE, "w" )

	# Write list of cell names
	open( IN_DIR * "cells.lst" ) do f
		hf["cells"] = readlines(f)
	end

	# Create groups
	grp_val = create_group( hf, "values" )
	grp_idx = create_group( hf, "indices" )


	# Read chromosome lengths
	chroms = []
	open( IN_DIR * "chromosomes.csv" ) do f
		prev_pos = 0
		while !eof(f)
			fields = split( readline(f), "," )
			pos = parse( Int, fields[2] )
			push!( chroms, ( fields[1], pos - prev_pos ) )
			prev_pos = pos
		end
	end
	
	# Write indices

	for (chrom, len) in chroms
		ds_chrom = create_dataset( grp_idx, chrom, datatype(UInt32), (len,2) )
		open( IN_DIR * "values.bin" ) do f
			for i in 1:len
				ds_chrom[ i, 1 ] = read( f, UInt32 )
				ds_chrom[ i, 2 ] = read( f, UInt32 )
			end
		end
		println( "Index for chromosome ", chrom, " written." )		
	end


	# Write values

	num_vals = Int( stat( IN_DIR * "values.bin" ).size / 4 )

	ds_cell = create_dataset( grp_val, "cell", datatype(UInt16), (num_vals,) )
	ds_umeth = create_dataset( grp_val, "count_unmeth", datatype(UInt8), (num_vals,) )
	ds_meth = create_dataset( grp_val, "count_meth", datatype(UInt8), (num_vals,) )

	open( IN_DIR * "values.bin" ) do f
		for i in 1:num_vals
			ds_cell[i] = read( f, UInt16 )
			ds_umeth[i] = read( f, UInt8 )
			ds_meth[i] = read( f, UInt8 )
		end
	end
	# The abvove loop is slow; should be changed to work with chunks

	close( hf )

	println( "All values written." )		


end

main()