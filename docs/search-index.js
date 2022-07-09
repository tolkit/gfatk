var searchIndex = JSON.parse('{\
"gfatk":{"doc":"<code>gfatk</code> is a tool for Graphical Fragment Assembly (GFA) …","t":[0,0,0,0,0,0,0,0,0,0,0,0,0,5,5,5,5,5,0,5,0,12,12,3,3,3,11,11,11,11,11,11,11,11,11,11,11,12,11,12,11,11,11,11,11,11,11,11,11,11,12,12,12,12,11,11,11,11,11,11,12,11,12,11,11,11,11,11,11,11,11,11,12,12,3,3,17,5,11,11,11,11,11,11,11,11,11,11,11,11,5,5,11,5,11,11,11,11,11,11,11,11,5,5,5,5,5,5,4,13,3,3,13,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,11,11,11,11,12,5,5,5,11,12,11,11,11,11,11,11,11,11,11,11,11,11,12,13,4,13,13,3,3,11,11,11,11,11,11,11,11,11,11,12,12,11,11,11,11,11,11,12,12,12,11,11,11,12,12,11,11,12,5,11,11,12,11,11,11,11,11,11,11,11,11,5,12,3,3,11,11,11,11,11,11,11,11,11,11,11,5,11,11,5,5,5,11,11,5,11,12,11,5,5,11,5,12,11,5,11,11,11,11,11,11,11,11,11],"n":["dot","extract","extract_chloro","extract_mito","fasta","gfa","linear","load","overlap","path","stats","trim","utils","dot","extract","extract_chloro","extract_mito","fasta","gfa","gfa_string","graph","0","0","GFAtk","Overlap","Overlaps","borrow","borrow","borrow","borrow_mut","borrow_mut","borrow_mut","clone","clone_into","from","from","from","from_orient","from_path_cli","from_segment","gen_cov_hash","get_coverage","into","into","into","into_digraph","into_ungraph","make_overlaps","new","node_seq_len_and_cov","overlap_str_from_f","overlap_str_from_r","overlap_str_to_f","overlap_str_to_r","parse_coverage_opt","print","print_extract","print_sequences","push","sequence_stats","to_orient","to_owned","to_segment","try_from","try_from","try_from","try_into","try_into","try_into","type_id","type_id","type_id","0","0","GFAdigraph","GFAungraph","MAX_RECURSION_DEPTH","all_paths","all_paths_all_node_pairs","borrow","borrow","borrow_mut","borrow_mut","dot","edge_count","from","from","into","into","node_count","recursive_path_finder_incl_coverage","recursive_path_finder_no_coverage","recursive_search","segments_subgraph","trim","try_from","try_from","try_into","try_into","type_id","type_id","weakly_connected_components","linear","linear_inner","byte_lines_iter","load_gfa","load_gfa_stdin","overlap","CLIOpt","File","GFAPath","GFAPathElement","String","borrow","borrow","borrow","borrow_mut","borrow_mut","borrow_mut","clone","clone","clone_into","clone_into","fmt","fmt","from","from","from","index","inner","into","into","into","new","orientation","parse_path","parse_path_string","path","push","segment_id","to_fasta_header","to_owned","to_owned","try_from","try_from","try_from","try_into","try_into","try_into","type_id","type_id","type_id","0","Chloroplast","GenomeType","Mitochondria","None","Stat","Stats","borrow","borrow","borrow","borrow_mut","borrow_mut","borrow_mut","clone","clone","clone_into","clone_into","cov","edge_count","eq","extract_organelle","fmt","from","from","from","gc","graph_indices_subgraph","index","into","into","into","is_circular","node_count","print_tabular","push","segments","stats","to_owned","to_owned","total_sequence_length","try_from","try_from","try_from","try_into","try_into","try_into","type_id","type_id","type_id","trim","0","GFAGraphLookups","GFAGraphPair","borrow","borrow","borrow_mut","borrow_mut","clone","clone","clone_into","clone_into","fmt","fmt","fmt","format_usize_to_kb","from","from","gc_content","get_edge_coverage","get_option_string","into","into","is_stdin","new","node_index","node_index_to_seg_id","nucleotide_counts","parse_cigar","push","reverse_complement","seg_id","seg_id_to_node_index","switch_base","to_owned","to_owned","to_string","try_from","try_from","try_into","try_into","type_id","type_id"],"q":["gfatk","","","","","","","","","","","","","gfatk::dot","gfatk::extract","gfatk::extract_chloro","gfatk::extract_mito","gfatk::fasta","gfatk::gfa","","","gfatk::gfa::gfa","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","gfatk::gfa::graph","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","gfatk::linear","","gfatk::load","","","gfatk::overlap","gfatk::path","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","gfatk::stats","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","gfatk::trim","gfatk::utils","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"d":["Make a DOT language representation of a GFA.","Extract a subgraph from a GFA.","Extract the putative plastid subgraph in a GFA.","Extract the putative mitochondrial subgraph in a GFA.","Print all the sequences in a GFA to fasta format.","A module with all the methods to manipulate GFA’s in.","Coerce a GFA into a fasta, finding the longest path …","Helper functions to load a GFA from a file, or read from …","Generate overlapping sequences between segments in a GFA.","Extract a fasta given a path.","Generate statistics about the input GFA file.","Utility to trim a GFA of isolated nodes.","Utility functions used throughout.","Make a DOT (https://graphviz.org/doc/info/lang.html) …","Supply a sequence/segment ID from the GFA, and extract the …","Using a combination of GC% of the segments, relative …","Using a combination of GC% of the segments, relative …","Print a fasta representation of the sequences in a GFA.","A module where all the methods of GFA manipulations are …","Writes a GFA to a string.","A module where a GFA is coerced into a petgraph <code>Graph</code> …","","","A wrapper around GFA from the gfa crate","Overlap from one segment to another.","A vector of <code>Overlap</code> structs.","","","","","","","","","Returns the argument unchanged.","Returns the argument unchanged.","Returns the argument unchanged.","Orientation of from segment.","Take a <code>GFAPath</code> and print out the path from a GFA.","ID of from segment.","Returns a <code>HashMap</code> of relative coverage of each node …","","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Returns a tuple of GFAGraphLookups (a struct of …","Returns a tuple of GFAGraphLookups (a struct of …","Returns the overlaps between all the segments in a GFA.","Create a new instance of <code>Overlaps</code>.","Return the coverage and sequence length for a segment, …","From segment forward.","From segment reverse.","To segment forward.","To segment reverse.","Two internal functions below to parse coverage of a GFA …","Print overlaps to STDOUT.","A method to print a GFA to STDOUT, given a vector of …","The internal function called when <code>gfatk fasta</code> is called.","Append to <code>Overlaps</code>, adding another <code>Overlap</code>.","The internal function called in <code>gfatk stats</code>.","Orientation of to segment.","","ID of to segment.","","","","","","","","","","","","A wrapper of petgraph’s directed <code>Graph</code> struct, applied …","A wrapper of petgraph’s undirected <code>Graph</code> struct, applied …","A recursion depth limit, so we don’t hit a stack overflow","A function generic over certain types of <code>Directed</code> petgraph …","The main function called from <code>gfatk linear</code>.","","","","","The main function called from <code>gfatk dot</code>.","Simple wrapper of <code>Graph.edge_count()</code> in petgraph.","Returns the argument unchanged.","Returns the argument unchanged.","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Simple wrapper of <code>Graph.node_count()</code> in petgraph.","Function called by <code>all_paths</code> where a <code>HashMap</code> is supplied …","A safer and more reliable alternative to …","The algorithm called in <code>gfatk extract</code>.","Returns a subgraph GFA that only contains elements with …","Trim a graph to include only nodes connected to two or …","","","","","","","Split the GFA digraph into subgraphs which are the weakly …","Force a linear representation of the GFA.","Reusable function to call on subgraphs in a GFA if …","Iterate over the byte lines of a file.","Given a path, load the GFA into a <code>GFA</code> struct.","If the file is coming from STDIN, this function reads a …","Generate overlaps between segments, with an optional …","Which option is used on the CLI, either a string or a file …","From a file.","A series of GFA path elements.","A GFA path element. Of the form <code>&lt;segment ID&gt;&lt;+/-&gt;</code>","From the command line.","","","","","","","","","","","","","Returns the argument unchanged.","Returns the argument unchanged.","Returns the argument unchanged.","The index of the struct.","","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Constructor for <code>GFAPath</code>.","The orientation of the segment.","Parse either a string, or a file, containing the path.","Parses a path string to a <code>GFAPath</code> object.","A function to generate a path through a GFA file given an …","Push an element to the end of the <code>GFAPath</code> buffer.","The ID of the GFA path element. It must match that of the …","Convert to a string for inclusion in the fasta header.","","","","","","","","","","","","","The chloroplast/plastid genome","Enumeration of the genomes we are interested in.","The mitochondrial genome","This will process the stats and return nothing","The statistics associated with a subgraph in a GFA.","A vector of <code>Stat</code>.","","","","","","","","","","","The average coverage across a subgraph.","The edge count of the graph.","","Extract the putative mitochondrial/chloroplast genome from …","","Returns the argument unchanged.","Returns the argument unchanged.","Returns the argument unchanged.","The average GC% across a subgraph.","The segments of the (sub)graph.","Arbitrary index of the subgraph(s).","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Whether the subgraph is circular (only applies to …","The node count of the graph.","Print tabular form of <code>Stats</code> to STDOUT.","Add a new <code>Stat</code> to <code>Stats</code>.","Names of the segments.","Internal function called in <code>gfatk stats</code>.","","","Total sequence length of all the segments.","","","","","","","","","","Trim a GFA file of segments which are connected only to …","","A vector of <code>GFAGraphPair</code>’s.","A pair consisting of a node index and a segment ID.","","","","","","","","","","","","Format a sequence length (<code>usize</code>) to kilobases.","Returns the argument unchanged.","Returns the argument unchanged.","Calculate the GC content of a string slice.","Get the coverage associated with an edge (<code>ec</code> tag in the …","Format a GFA option field into a string.","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Check if there is anything coming from STDIN.","Create a new GFAGraphLookups","The node index (petgraph’s <code>NodeIndex</code>).","Return segment ID from a node index.","Collect nucleotide counts into a <code>HashMap</code>.","Parse a CIGAR string slice into an overlap length.","Push a new <code>GFAGraphPair</code> to the end.","Reverse complement a string slice.","The segment ID.","Return a node index from a segment ID.","Used in <code>reverse_complement</code> to switch to a complementary …","","","","","","","","",""],"i":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0,3,1,2,3,1,2,2,2,3,1,2,3,2,3,2,2,3,1,2,2,2,2,1,2,3,3,3,3,2,1,2,2,1,2,3,2,3,3,1,2,3,1,2,3,1,2,4,5,0,0,0,0,5,4,5,4,5,5,5,4,5,4,5,5,0,0,4,0,5,4,5,4,5,4,5,5,0,0,0,0,0,0,0,6,0,0,6,6,7,8,6,7,8,7,8,7,8,7,8,6,7,8,7,8,6,7,8,8,7,0,0,0,8,7,8,7,8,6,7,8,6,7,8,6,7,8,9,10,0,10,10,0,0,9,10,11,9,10,11,10,11,10,11,11,11,10,9,11,9,10,11,11,11,11,9,10,11,11,11,9,9,11,0,10,11,11,9,10,11,9,10,11,9,10,11,0,12,0,0,13,12,13,12,13,12,13,12,13,12,12,0,13,12,0,0,0,13,12,0,12,13,12,0,0,12,0,13,12,0,13,12,12,13,12,13,12,13,12],"f":[null,null,null,null,null,null,null,null,null,null,null,null,null,[[["argmatches",3]],["result",6]],[[["argmatches",3]],["result",6]],[[["argmatches",3],["genometype",4]],["result",6]],[[["argmatches",3],["genometype",4]],["result",6]],[[["argmatches",3]],["result",6]],null,[[["gfa",3]],["string",3]],null,null,null,null,null,null,[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["gfatk",3]],[[["",0],["",0]]],[[]],[[]],[[]],null,[[["",0],["gfapath",3],["hashmap",3,[["string",3],["usize",0]]],["str",0],["option",4,[["str",0]]]],["result",6]],null,[[["",0],["gfagraphlookups",3]],["result",6,[["hashmap",3,[["nodeindex",3],["usize",0]]]]]],[[["",0]],["result",6,[["f32",0]]]],[[]],[[]],[[]],[[["",0]],["result",6]],[[["",0]],["result",6]],[[["",0],["usize",0]],["result",6,[["overlaps",3]]]],[[]],[[["",0],["usize",0]],["result",6]],null,null,null,null,[[["optfieldval",4]],["result",6,[["f32",0]]]],[[["usize",0]]],[[["",0],["vec",3,[["usize",0]]]]],[[["",0],["option",4,[["string",3]]]],["result",6]],[[["",0],["overlap",3]]],[[["",0],["genometype",4],["bool",0]],["result",6]],null,[[["",0]]],null,[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[["",0]],["typeid",3]],[[["",0]],["typeid",3]],[[["",0]],["typeid",3]],null,null,null,null,null,[[["graph",3],["nodeindex",3,[["indextype",8]]],["nodeindex",3,[["indextype",8]]],["option",4,[["hashmap",3]]],["usize",0]],["result",6,[["vec",3,[["vec",3,[["nodeindex",3,[["indextype",8]]]]]]]]]],[[["",0],["gfagraphlookups",3],["option",4,[["hashmap",3]]]],["result",6]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0],["gfatk",3]],["result",6]],[[["",0]],["usize",0]],[[]],[[]],[[]],[[]],[[["",0]],["usize",0]],[[["graph",3],["nodeindex",3,[["indextype",8]]],["nodeindex",3,[["indextype",8]]],["hashmap",3],["hashmap",3],["usize",0],["usize",0]],["option",4,[["result",6,[["vec",3,[["vec",3,[["nodeindex",3,[["indextype",8]]]]]]]]]]]],[[["graph",3],["nodeindex",3,[["indextype",8]]],["nodeindex",3,[["indextype",8]]],["hashset",3]],["result",6,[["vec",3,[["vec",3,[["nodeindex",3,[["indextype",8]]]]]]]]]],[[["",0],["vec",3,[["usize",0]]],["i32",0],["vec",3,[["nodeindex",3]]],["gfagraphlookups",3]],["result",6,[["vec",3,[["usize",0]]]]]],[[["gfa",3],["vec",3,[["usize",0]]]],["gfa",3,[["usize",0],["",26,[["optfields",8],["clone",8]]]]]],[[["",0],["gfagraphlookups",3]],["vec",3,[["usize",0]]]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[["",0]],["typeid",3]],[[["",0]],["typeid",3]],[[["",0],["gfagraphlookups",3]],["result",6,[["vec",3,[["vec",3,[["usize",0]]]]]]]],[[["argmatches",3]],["result",6]],[[["gfatk",3],["bool",0],["gfagraphlookups",3],["gfadigraph",3],["option",4,[["string",3]]]],["result",6]],[[["read",8]],["box",3,[["iterator",8]]]],[[],["result",6,[["gfa",3]]]],[[["stdinlock",3]],["result",6,[["gfa",3],["parseerror",4]]]],[[["argmatches",3]],["result",6]],null,null,null,null,null,[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["gfapathelement",3]],[[["",0]],["gfapath",3]],[[["",0],["",0]]],[[["",0],["",0]]],[[["",0],["formatter",3]],["result",6]],[[["",0],["formatter",3]],["result",6]],[[]],[[]],[[]],null,null,[[]],[[]],[[]],[[]],null,[[["str",0],["cliopt",4],["gfatk",3]],["result",6]],[[["str",0],["gfatk",3]],["result",6]],[[["argmatches",3]],["result",6]],[[["",0],["gfapathelement",3]]],null,[[["",0]],["string",3]],[[["",0]]],[[["",0]]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[["",0]],["typeid",3]],[[["",0]],["typeid",3]],[[["",0]],["typeid",3]],null,null,null,null,null,null,null,[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["genometype",4]],[[["",0]],["stat",3]],[[["",0],["",0]]],[[["",0],["",0]]],null,null,[[["",0],["genometype",4]],["bool",0]],[[["",0],["usize",0],["usize",0],["f32",0],["f32",0]],["result",6,[["vec",3,[["usize",0]]]]]],[[["",0],["formatter",3]],["result",6]],[[]],[[]],[[]],null,null,null,[[]],[[]],[[]],null,null,[[["",0]]],[[["",0],["stat",3]]],null,[[["argmatches",3],["genometype",4]],["result",6,[["option",4]]]],[[["",0]]],[[["",0]]],null,[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[["",0]],["typeid",3]],[[["",0]],["typeid",3]],[[["",0]],["typeid",3]],[[["argmatches",3]],["result",6]],null,null,null,[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["",0]],[[["",0]],["gfagraphpair",3]],[[["",0]],["gfagraphlookups",3]],[[["",0],["",0]]],[[["",0],["",0]]],[[["",0],["formatter",3]],["result",6]],[[["",0],["formatter",3]],["result",6]],[[["",0],["formatter",3]],["result",6]],[[["usize",0]],["string",3]],[[]],[[]],[[],["f32",0]],[[],["result",6,[["i64",0]]]],[[["vec",3,[["optfield",3]]]],["result",6,[["string",3]]]],[[]],[[]],[[],["bool",0]],[[]],null,[[["",0],["nodeindex",3]],["result",6,[["usize",0]]]],[[],["hashmap",3,[["u8",0],["i32",0]]]],[[],["result",6,[["usize",0]]]],[[["",0],["gfagraphpair",3]]],[[],["vec",3,[["u8",0]]]],null,[[["",0],["usize",0]],["result",6,[["nodeindex",3]]]],[[["u8",0]],["u8",0]],[[["",0]]],[[["",0]]],[[["",0]],["string",3]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[["",0]],["typeid",3]],[[["",0]],["typeid",3]]],"p":[[3,"Overlaps"],[3,"GFAtk"],[3,"Overlap"],[3,"GFAungraph"],[3,"GFAdigraph"],[4,"CLIOpt"],[3,"GFAPathElement"],[3,"GFAPath"],[3,"Stats"],[4,"GenomeType"],[3,"Stat"],[3,"GFAGraphLookups"],[3,"GFAGraphPair"]]}\
}');
if (window.initSearch) {window.initSearch(searchIndex)};