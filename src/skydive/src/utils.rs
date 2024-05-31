use std::collections::HashSet;

use anyhow::Result;

pub fn parse_loci(loci_list: &Vec<String>) -> HashSet<(String, u64, u64)> {
    // Initialize a HashSet to store unique loci after parsing
    let mut loci = HashSet::new();
    
    // Iterate over each locus in the provided list
    for locus in loci_list {
        // Attempt to parse the locus using a function from the skydive module
        match parse_locus(locus.to_owned()) {
            Ok(l_fmt) => {
                // If parsing is successful, insert the formatted locus into the HashSet
                loci.insert(l_fmt);
            }
            Err(_) => {
                // If parsing fails, panic and terminate the program, providing an error message
                panic!("Could not parse locus '{}'.", locus);
            }
        }
    }

    loci
}

pub fn parse_locus(locus: String) -> Result<(String, u64, u64)> {
    let l_fmt = locus.replace(",", "");
    let parts: Vec<&str> = l_fmt.split(|c| (c == ':' || c == '-')).collect();

    let chr = parts[0].to_string();

    if parts.len() == 2 {
        let start = match parts[1].parse::<u64>() {
            Ok(val) => val,
            Err(e) => {
                return Err(anyhow::Error::new(e));
            }
        };

        Ok((chr, start - 1000, start + 1000))
    } else if parts.len() == 3 {
        let start = match parts[1].parse::<u64>() {
            Ok(val) => val,
            Err(e) => {
                return Err(anyhow::Error::new(e));
            }
        };

        let stop = match parts[2].parse::<u64>() {
            Ok(val) => val,
            Err(e) => {
                return Err(anyhow::Error::new(e));
            }
        };

        Ok((chr, start, stop))
    } else {
        anyhow::bail!("Locus format for '{}' is incorrect. It should be 'chr:start[-stop]'.", locus);
    }
}