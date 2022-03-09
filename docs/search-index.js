var searchIndex = JSON.parse('{\
"gfatk":{"doc":"","t":[0,0,0,0,0,0,0,0,0,0,5,5,5,5,0,0,0,12,3,3,3,11,11,11,11,11,11,11,11,11,11,12,12,11,11,11,11,11,11,11,11,12,12,12,12,11,11,11,11,11,12,12,11,11,11,11,11,11,11,11,11,12,12,3,3,5,11,11,11,11,11,11,11,11,11,11,11,11,11,5,11,11,11,11,11,11,11,5,5,5,5,5,5,12,3,3,11,11,11,11,11,11,12,11,11,11,12,12,11,11,11,12,5,11,12,11,11,11,11,11,11,12,3,3,11,11,11,11,11,11,11,11,11,5,11,11,5,5,5,11,11,5,11,12,11,5,11,5,12,11,11,11,11,11,11,11,11,11,11],"n":["dot","extract","extract_mito","fasta","gfa","linear","load","overlap","stats","utils","dot","extract","extract_mito","fasta","gfa","graph","writer","0","GFAtk","Overlap","Overlaps","borrow","borrow","borrow","borrow_mut","borrow_mut","borrow_mut","determine_path_overlaps","from","from","from","from_orient","from_segment","gen_cov_hash","into","into","into","into_digraph","into_ungraph","make_overlaps","node_seq_len_and_cov","overlap_str_from_f","overlap_str_from_r","overlap_str_to_f","overlap_str_to_r","print","print_extract","print_path_to_fasta","print_sequences","sequence_stats","to_orient","to_segment","try_from","try_from","try_from","try_into","try_into","try_into","type_id","type_id","type_id","0","0","GFAdigraph","GFAungraph","all_paths","all_paths_all_node_pairs","borrow","borrow","borrow_mut","borrow_mut","dot","edge_count","from","from","into","into","node_count","recursive_search","segments_subgraph","try_from","try_from","try_into","try_into","type_id","type_id","weakly_connected_components","gfa_string","force_linear","byte_lines_iter","load_gfa","load_gfa_stdin","overlap","0","Stat","Stats","borrow","borrow","borrow_mut","borrow_mut","clone","clone_into","cov","extract_mito","from","from","gc","index","into","into","push","segments","stats","to_owned","total_sequence_length","try_from","try_from","try_into","try_into","type_id","type_id","0","GFAGraphLookups","GFAGraphPair","borrow","borrow","borrow_mut","borrow_mut","clone","clone","clone_into","clone_into","fmt","format_usize_to_kb","from","from","gc_content","get_edge_coverage","get_option_string","into","into","is_stdin","new","node_index","node_index_to_seg_id","parse_cigar","push","reverse_complement","seg_id","seg_id_to_node_index","to_owned","to_owned","to_string","try_from","try_from","try_into","try_into","type_id","type_id"],"q":["gfatk","","","","","","","","","","gfatk::dot","gfatk::extract","gfatk::extract_mito","gfatk::fasta","gfatk::gfa","","","gfatk::gfa::gfa","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","gfatk::gfa::graph","","","","","","","","","","","","","","","","","","","","","","","","","","gfatk::gfa::writer","gfatk::linear","gfatk::load","","","gfatk::overlap","gfatk::stats","","","","","","","","","","","","","","","","","","","","","","","","","","","","gfatk::utils","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"d":["Make a DOT language representation of a GFA.","Extract a subgraph from a GFA.","Extract the putative mitochondrial subgraph in a GFA.","Print all the sequences in a GFA to fasta format.","A module with all the methods to manipulate GFA’s in.","Coerce a GFA into a fasta, finding the longest path …","Helper functions to load a GFA from a file, or read from …","Generate overlapping sequences between segments in a GFA.","Generate statistics about the input GFA file.","Utility functions used throughout.","Make a DOT (https://graphviz.org/doc/info/lang.html) …","Supply a sequence/segment ID from the GFA, and extract the …","Using a combination of GC% of the segments, relative …","Print a fasta representation of the sequences in a GFA.","A module where all the methods of GFA manipulations are …","A module where a GFA is coerced into a petgraph <code>Graph</code> …","Simple modified wrapper for a function in the gfa crate.","","A wrapper around GFA from the gfa crate","Overlap from one segment to another.","A vector of <code>Overlap</code> structs.","","","","","","","Computes the overlaps between segments as a vector of: …","","","","Orientation of from segment.","ID of from segment.","Returns a <code>HashMap</code> of relative coverage of each node …","","","","Returns a tuple of GFAGraphLookups (a struct of …","Returns a tuple of GFAGraphLookups (a struct of …","Returns the overlaps between all the segments in a GFA.","Return the coverage and sequence length for a segment, …","From segment forward.","From segment reverse.","To segment forward.","To segment reverse.","Print overlaps to STDOUT.","A method to print a GFA to STDOUT, given a vector of …","A method to print to STDOUT a fasta, given a path through …","The internal function called when <code>gfatk fasta</code> is called.","The internal function called in <code>gfatk stats</code>.","Orientation of to segment.","ID of to segment.","","","","","","","","","","","","A wrapper of petgraph’s directed <code>Graph</code> struct, applied …","A wrapper of petgraph’s undirected <code>Graph</code> struct, applied …","A function generic over certain types of <code>Directed</code> petgraph …","","","","","","","","","","","","","","Returns a subgraph GFA that only contains elements with …","","","","","","","","Writes a GFA to a string.","Force a linear representation of the GFA.","Iterate over the byte lines of a file.","Given a path, load the GFA into a <code>GFA</code> struct.","If the file is coming from STDIN, this function reads a …","Generate overlaps between segments, with an optional …","","The statistics associated with a subgraph in a GFA.","A vector of <code>Stat</code>.","","","","","","","The average coverage across a subgraph.","The function called from <code>gfatk extract-mito</code>.","","","The average GC% across a subgraph.","Arbitrary index of the subgraph(s).","","","Add a new <code>Stat</code> to <code>Stats</code>.","Names of the segments.","Internal function called in <code>gfatk stats</code>.","","Total sequence length of all the segments.","","","","","","","","A vector of <code>GFAGraphPair</code>’s.","A pair consisting of a node index and a segment ID.","","","","","","","","","","Format a sequence length (<code>usize</code>) to kilobases.","","","Calculate the GC content of a string slice.","Get the coverage associated with an edge (<code>ec</code> tag in the …","Format a GFA option field into a string.","","","Check if there is anything coming from STDIN.","Create a new GFAGraphLookups","The node index (petgraph’s <code>NodeIndex</code>).","Return segment ID from a node index.","Parse a CIGAR string slice into an overlap length.","Push a new <code>GFAGraphPair</code> to the end.","Reverse complement a string slice.","The segment ID.","Return a node index from a segment ID.","","","","","","","","",""],"i":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,2,3,1,2,3,1,1,2,3,2,2,1,1,2,3,1,1,1,1,2,2,2,2,3,1,1,1,1,2,2,1,2,3,1,2,3,1,2,3,4,5,0,0,0,5,4,5,4,5,5,5,4,5,4,5,5,4,0,4,5,4,5,4,5,5,0,0,0,0,0,0,6,0,0,6,7,6,7,7,7,7,6,6,7,7,7,6,7,6,7,0,7,7,6,7,6,7,6,7,8,0,0,9,8,9,8,9,8,9,8,8,0,9,8,0,0,0,9,8,0,8,9,8,0,8,0,9,8,9,8,8,9,8,9,8,9,8],"f":[null,null,null,null,null,null,null,null,null,null,[[["argmatches",3]],["result",6]],[[["argmatches",3]],["result",6]],[[["argmatches",3]],["result",6]],[[["argmatches",3]],["result",6]],null,null,null,null,null,null,null,[[]],[[]],[[]],[[]],[[]],[[]],[[["vec",3],["gfagraphlookups",3],["vec",3]],["result",6,[["vec",3]]]],[[]],[[]],[[]],null,null,[[["gfagraphlookups",3]],["result",6,[["hashmap",3,[["nodeindex",3],["usize",15]]]]]],[[]],[[]],[[]],[[],["result",6]],[[],["result",6]],[[["usize",15]],["result",6,[["overlaps",3]]]],[[["usize",15]],["result",6]],null,null,null,null,[[["usize",15]]],[[["vec",3,[["usize",15]]]]],[[["indexmap",3,[["usize",15],["vec",3]]],["str",15],["vec",3,[["usize",15]]]],["result",6]],[[],["result",6]],[[["bool",15]],["result",6]],null,null,[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["typeid",3]],[[],["typeid",3]],[[],["typeid",3]],null,null,null,null,[[["graph",3],["nodeindex",3,[["indextype",8]]],["nodeindex",3,[["indextype",8]]],["option",4,[["hashmap",3]]]],["result",6,[["vec",3,[["vec",3,[["nodeindex",3,[["indextype",8]]]]]]]]]],[[["gfagraphlookups",3],["option",4,[["hashmap",3]]]],["result",6]],[[]],[[]],[[]],[[]],[[["gfatk",3]],["result",6]],[[],["usize",15]],[[]],[[]],[[]],[[]],[[],["usize",15]],[[["usize",15],["i32",15],["vec",3,[["nodeindex",3]]],["gfagraphlookups",3]],["result",6,[["vec",3,[["usize",15]]]]]],[[["gfa",3],["vec",3,[["usize",15]]]],["gfa",3,[["usize",15],["",26,[["optfields",8],["clone",8]]]]]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["typeid",3]],[[],["typeid",3]],[[["gfagraphlookups",3]],["result",6,[["vec",3,[["vec",3,[["usize",15]]]]]]]],[[["gfa",3]],["string",3]],[[["argmatches",3]],["result",6]],[[["read",8]],["box",3,[["iterator",8]]]],[[],["result",6,[["gfa",3]]]],[[["stdinlock",3]],["result",6,[["gfa",3],["parseerror",4]]]],[[["argmatches",3]],["result",6]],null,null,null,[[]],[[]],[[]],[[]],[[],["stat",3]],[[]],null,[[],["vec",3,[["usize",15]]]],[[]],[[]],null,null,[[]],[[]],[[["stat",3]]],null,[[["argmatches",3],["bool",15]],["result",6,[["option",4]]]],[[]],null,[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["typeid",3]],[[],["typeid",3]],null,null,null,[[]],[[]],[[]],[[]],[[],["gfagraphpair",3]],[[],["gfagraphlookups",3]],[[]],[[]],[[["formatter",3]],["result",6]],[[["usize",15]],["string",3]],[[]],[[]],[[],["f32",15]],[[["vec",3]],["result",6,[["i64",15]]]],[[["vec",3,[["optfield",3]]]],["result",6,[["string",3]]]],[[]],[[]],[[],["bool",15]],[[]],null,[[["nodeindex",3]],["result",6,[["usize",15]]]],[[],["result",6,[["usize",15]]]],[[["gfagraphpair",3]]],[[],["vec",3,[["u8",15]]]],null,[[["usize",15]],["result",6,[["nodeindex",3]]]],[[]],[[]],[[],["string",3]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["typeid",3]],[[],["typeid",3]]],"p":[[3,"GFAtk"],[3,"Overlap"],[3,"Overlaps"],[3,"GFAungraph"],[3,"GFAdigraph"],[3,"Stats"],[3,"Stat"],[3,"GFAGraphLookups"],[3,"GFAGraphPair"]]}\
}');
if (window.initSearch) {window.initSearch(searchIndex)};