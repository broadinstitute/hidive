use anyhow::Result;

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