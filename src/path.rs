use crate::gfa::gfa::GFAtk;
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;

use anyhow::{bail, ensure, Context, Result};
use gfa::gfa::Orientation;
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

/// Which option is used on the CLI, either a string
/// or a file path.
pub enum CLIOpt {
    /// From the command line.
    String,
    /// From a file.
    File,
}

/// A function to generate a path through a GFA file
/// given an input path.
///
/// This application is sufficiently different from
/// the auto-generated paths of `gfatk linear`, so it is
/// not included there.
///
/// For example:
/// ```bash
/// gfatk path ./input.gfa "12+, 11+, 2-, 2+"
/// ```
pub fn path(matches: &clap::ArgMatches) -> Result<()> {
    // read in path and parse gfa
    let gfa_file = matches.get_one::<PathBuf>("GFA");
    let path_cli = matches.get_one::<String>("path_cli");
    let path_file = matches.get_one::<PathBuf>("path_file");
    let all_p_lines = matches.get_flag("all_paths");

    let gfa: GFAtk = match gfa_file {
        Some(f) => {
            let ext = f.extension();
            match ext {
                Some(e) => {
                    if e == "gfa" {
                        GFAtk(load_gfa(f)?)
                    } else {
                        bail!("Input is not a GFA.")
                    }
                }
                None => bail!("Could not read file."),
            }
        }
        None => match utils::is_stdin() {
            true => GFAtk(load_gfa_stdin(std::io::stdin().lock())?),
            false => bail!("No input from STDIN. Run `gfatk path -h` for help."),
        },
    };

    if all_p_lines {
        let paths = gfa.get_path_lines()?;
        for (id, path) in paths {
            let (parsed_path, link_map) = parse_path(&path, CLIOpt::String, &gfa)?;
            gfa.from_path_cli(parsed_path, link_map, "path_all", Some(id))?;
        }
    } else {
        // we need some path specified
        if path_cli.is_none() && path_file.is_none() {
            bail!(
            "Please specify either a path as a positional argument string, or as a file `--path`."
        )
        }
        // but not both!
        if path_cli.is_some() && path_file.is_some() {
            bail!("Specify either <path>, or `--path`, not both.")
        }

        let (path, link_map) = match path_cli {
            Some(p) => parse_path(p, CLIOpt::String, &gfa),
            None => match path_file {
                Some(f) => match f.to_str() {
                    Some(string) => parse_path(string, CLIOpt::File, &gfa),
                    None => bail!("Could not convert 'path' file path to string."),
                },
                None => bail!("Should never reach here."),
            },
        }?;

        gfa.from_path_cli(path, link_map, "path", None)?;
    }

    Ok(())
}

/// Parse either a string, or a file, containing
/// the path.
pub fn parse_path(
    path: &str,
    is_cli: CLIOpt,
    gfa: &GFAtk,
) -> Result<(GFAPath, HashMap<String, usize>)> {
    match is_cli {
        CLIOpt::String => parse_path_string(path, gfa),
        CLIOpt::File => {
            let mut first_line = String::new();

            let file = fs::File::open(path)?;
            let mut buffer = BufReader::new(file);

            buffer.read_line(&mut first_line)?;

            parse_path_string(&first_line, gfa)
        }
    }
}

/// A GFA path element. Of the form `<segment ID><+/->`
#[derive(Debug, Clone)]
pub struct GFAPathElement {
    /// The ID of the GFA path element.
    /// It must match that of the segments in the GFA.
    pub segment_id: Vec<u8>,
    /// The orientation of the segment.
    pub orientation: Orientation,
    /// The index of the struct.
    pub index: usize,
}

#[derive(Debug, Clone)]
/// A series of GFA path elements.
pub struct GFAPath {
    pub inner: Vec<GFAPathElement>,
}

impl GFAPath {
    /// Constructor for [`GFAPath`].
    fn new() -> Self {
        Self { inner: vec![] }
    }
    /// Push an element to the end of the [`GFAPath`] buffer.
    fn push(&mut self, other: GFAPathElement) {
        self.inner.push(other);
    }
    /// Convert to a string for inclusion in the fasta header.
    pub fn to_fasta_header(&self) -> String {
        let mut output = String::new();

        for el in &self.inner {
            let el_d = std::str::from_utf8(&el.segment_id).unwrap();
            output += &format!("{}{},", el_d, el.orientation);
        }
        output.pop();
        output
    }
}

/// Parses a path string to a [`GFAPath`] object.
fn parse_path_string(path_string: &str, gfa: &GFAtk) -> Result<(GFAPath, HashMap<String, usize>)> {
    let gfa = &gfa.0;
    // make a map of the links
    let mut link_map = HashMap::new();
    for link in &gfa.links {
        let from_d = std::str::from_utf8(&link.from_segment)?;
        let to_d = std::str::from_utf8(&link.to_segment)?;
        let path_pair = format!("{}{}|{}{}", from_d, link.from_orient, to_d, link.to_orient);
        let cigar = utils::parse_cigar(&link.overlap)?;

        link_map.insert(path_pair, cigar);
    }

    // path_string consists of e.g.:
    // 1+, 2-, 3+, 4-, ...
    // split
    let mut split_path: Vec<&str> = path_string.split(',').collect();
    // trim whitespace
    for token in split_path.iter_mut() {
        *token = token.trim();
    }
    // now allocate to GFAPath
    let mut gfa_path = GFAPath::new();

    for (index, token) in split_path.into_iter().enumerate() {
        let mut token_string = token.to_owned();
        let orientation = token_string
            .pop()
            .context("Each path element should contain a character.")?;

        ensure!(
            orientation == '-' || orientation == '+',
            "The last char in the split_path token was {}, not \'+\' or \'-\'. Check path is specified correctly.",
            orientation
        );

        let o_enum = match orientation {
            '+' => Orientation::Forward,
            '-' => Orientation::Backward,
            _ => bail!("Orientation can only be the chars \'+\' or \'-\'."),
        };

        gfa_path.push(GFAPathElement {
            segment_id: token_string.into_bytes(),
            orientation: o_enum,
            index,
        });
    }

    Ok((gfa_path, link_map))
}
