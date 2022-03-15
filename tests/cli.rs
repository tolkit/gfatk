// Command line tests for `gfatk`
// run with `cargo test --release`

use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;

// test `gfatk linear`

// H	VN:Z:1.0
// S	11	ACCTT	ll:f:30.0
// S	12	TCAAGG	ll:f:60.0
// S	13	CTTGATT	ll:f:30.0
// L	11	+	12	-	4M	ec:i:1
// L	12	-	13	+	5M	ec:i:1
// L	11	+	13	+	3M	ec:i:1
// L	12	+	11	-	4M	ec:i:1
// L	13	-	12	+	5M	ec:i:1
// L	13	-	11	-	3M	ec:i:1

// Expected output from gfatk linear
// There are two possible paths that could
// be output with the same probability.

// PATH: 11+ -> 12- -> 13+
// 11 ACCTT
// 12  CCTTGA
// 13   CTTGATT
// P ACCTTGATT

// PATH: 13- -> 12+ -> 11-
// 13 AATCAAG
// 12   TCAAGG
// 11     AAGGT
// P  AATCAAGGT

// TODO: test `-i` flag

#[test]
fn test_gfa_linear_stdout() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfatk")?;

    cmd.arg("linear").arg("./tests/test_linear.gfa");
    cmd.assert()
        .stdout(predicate::str::contains("ACCTTGATT").or(predicate::str::contains("AATCAAGGT")));

    Ok(())
}

// test `gfatk overlap`

// same test GFA as `gfatk linear`

// H	VN:Z:1.0
// S	11	ACCTT	ll:f:30.0
// S	12	TCAAGG	ll:f:60.0
// S	13	CTTGATT	ll:f:30.0
// L	11	+	12	-	4M	ec:i:1
// L	12	-	13	+	5M	ec:i:1
// L	11	+	13	+	3M	ec:i:1
// L	12	+	11	-	4M	ec:i:1
// L	13	-	12	+	5M	ec:i:1
// L	13	-	11	-	3M	ec:i:1

// walking through the links above
// the overlaps with one base pair either side should be

// 11 + -> 12 - : Overlap = CCTT; -1 = A; +1 = G --> ACCTTG
// 12 - -> 13 + : Overlap = CTTGA; -1 = C; +1 = T --> CCTTGAT
// 11 + -> 13 + : Overlap = CTT; -1 = C; +1 = G --> CCTTG
// 12 + -> 11 - : Overlap = AAGG; -1 = C; +1 = T --> CAAGGT
// 13 - -> 12 + : Overlap = TCAAG; -1 = A; +1 = G --> ATCAAGG
// 13 - -> 11 - : Overlap = AAG; -1 = C; +1 = G --> CAAGG

#[test]
fn test_gfa_overlap_stdout() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfatk")?;

    // just one base either side
    cmd.arg("overlap")
        .arg("./tests/test_linear.gfa")
        .arg("-s")
        .arg("1");

    // test all the overlaps are present.
    cmd.assert().stdout(predicate::str::contains("ACCTTG").and(
        predicate::str::contains("CCTTGAT").and(
            predicate::str::contains("CCTTG").and(
                predicate::str::contains("CAAGGT").and(
                    predicate::str::contains("ATCAAGG").and(predicate::str::contains("CAAGG")),
                ),
            ),
        ),
    ));

    Ok(())
}

// test `gfatk extract`

// # Duplicate 11,12,13
// # To 14,15,16
// H	VN:Z:1.0
// S	11	ACCTT	ll:f:30.0
// S	12	TCAAGG	ll:f:60.0
// S	13	CTTGATT	ll:f:30.0
// L	11	+	12	-	4M	ec:i:1
// L	12	-	13	+	5M	ec:i:1
// L	11	+	13	+	3M	ec:i:1
// L	12	+	11	-	4M	ec:i:1
// L	13	-	12	+	5M	ec:i:1
// L	13	-	11	-	3M	ec:i:1
// # duplication here
// S	14	ACCTT	ll:f:30.0
// S	15	TCAAGG	ll:f:60.0
// S	16	CTTGATT	ll:f:30.0
// L	14	+	15	-	4M	ec:i:1
// L	15	-	16	+	5M	ec:i:1
// L	14	+	16	+	3M	ec:i:1
// L	15	+	14	-	4M	ec:i:1
// L	16	-	15	+	5M	ec:i:1
// L	16	-	14	-	3M	ec:i:1

// two subgraphs present above, segment ID 11
// will extract a graph equivalent to
// ./test_linear.gfa

#[test]
fn test_subgraph_extraction() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfatk")?;

    // just one base either side
    cmd.arg("extract")
        .arg("./tests/test_subgraphs.gfa")
        .arg("-s")
        .arg("11");

    // should be the same output as ./tests/test_linear.gfa
    cmd.assert().stdout(predicate::str::contains(
        "H	VN:Z:1.0
S	11	ACCTT	ll:f:30
S	12	TCAAGG	ll:f:60
S	13	CTTGATT	ll:f:30
L	11	+	12	-	4M	ec:i:1
L	12	-	13	+	5M	ec:i:1
L	11	+	13	+	3M	ec:i:1
L	12	+	11	-	4M	ec:i:1
L	13	-	12	+	5M	ec:i:1
L	13	-	11	-	3M	ec:i:1
",
    ));

    Ok(())
}

// test `gfatk trim`

// H	VN:Z:1.0
// S	11	ACCTT	ll:f:30.0
// S	12	TCAAGG	ll:f:60.0
// S	13	CTTGATT	ll:f:30.0
// S	14	TTGGGG	ll:f:30.0
// L	11	+	12	-	4M	ec:i:1
// L	12	-	13	+	5M	ec:i:1
// L	11	+	13	+	3M	ec:i:1
// L	12	+	11	-	4M	ec:i:1
// L	13	-	12	+	5M	ec:i:1
// L	13	-	11	-	3M	ec:i:1
// L	14	+	11	+	2M	ec:i:1
// L	11	-	14	-	2M	ec:i:1

// added segment 14, with only two links
// therefore this should be removed

#[test]
fn test_gfa_trim_sterr() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfatk")?;

    cmd.arg("trim").arg("./tests/test_trim.gfa");
    cmd.assert().stderr(predicate::str::contains(
        "[+]	Removed segment 14 from GFA.
",
    ));

    Ok(())
}

// test for no edge coverage tags
// if user wants to use:
// `gfatk linear`, `gfatk dot`, `gfatk trim` or `gfatk stats`
// Edge coverage must be present.

#[test]
fn test_gfa_edge_coverage_failure_linear() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfatk")?;

    cmd.arg("linear").arg("./tests/test_no_ec.gfa");
    cmd.assert().failure();

    Ok(())
}

#[test]
fn test_gfa_edge_coverage_failure_dot() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfatk")?;

    cmd.arg("dot").arg("./tests/test_no_ec.gfa");
    cmd.assert().failure();

    Ok(())
}

#[test]
fn test_gfa_edge_coverage_failure_trim() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfatk")?;

    cmd.arg("trim").arg("./tests/test_no_ec.gfa");
    cmd.assert().failure();

    Ok(())
}

#[test]
fn test_gfa_edge_coverage_failure_stats() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfatk")?;

    cmd.arg("stats").arg("./tests/test_no_ec.gfa");
    cmd.assert().failure();

    Ok(())
}

// test segment coverage tag presence
// only relevant for:
// `gfatk linear -i <in.gfa>`

#[test]
fn test_gfa_node_coverage_failure_linear() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfatk")?;

    cmd.arg("linear").arg("./tests/test_no_ll.gfa").arg("-i");
    cmd.assert().failure();

    Ok(())
}
