//! This module implements the `WhatsHap` phasing algorithm as described in:
//! Murray Patterson, Tobias Marschall, Nadia Pisanti, Leo van Iersel, Leen Stougie, et al.
//! `WhatsHap`: Weighted Haplotype Assembly for Future-Generation Sequencing Reads.
//! Journal of Computational Biology, 2015, 22 (6), pp.498-509.
//! DOI: 10.1089/cmb.2014.0157
//! HAL ID: hal-01225988

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::Write;

use itertools::Itertools;

use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

#[derive(Debug)]
pub struct WMECData {
    pub reads: Vec<Vec<Option<u8>>>, // Reads matrix where None represents missing data
    pub confidences: Vec<Vec<Option<u32>>>, // Confidence degrees matrix
    pub num_snps: usize,             // Number of SNP positions
}

impl WMECData {
    // Initialize the data structure with given reads and confidences
    #[must_use]
    pub fn new(reads: Vec<Vec<Option<u8>>>, confidences: Vec<Vec<Option<u32>>>) -> Self {
        let num_snps = reads[0].len();
        WMECData {
            reads,
            confidences,
            num_snps,
        }
    }

    // Function to compute W^0(j, R) and W^1(j, R)
    // Cost to set all fragments in set R to 0 or 1 at SNP j
    #[must_use]
    pub fn compute_costs(&self, snp: usize, set_r: &BTreeSet<usize>) -> (u32, u32) {
        let mut w0 = 0; // Cost for setting to 0
        let mut w1 = 0; // Cost for setting to 1

        for &read_index in set_r {
            if let Some(allele) = self.reads[read_index][snp] {
                if let Some(confidence) = self.confidences[read_index][snp] {
                    if allele == 0 {
                        w1 += confidence; // Cost of flipping 0 -> 1
                    } else {
                        w0 += confidence; // Cost of flipping 1 -> 0
                    }
                }
            }
        }

        (w0, w1)
    }

    // Calculate minimum correction cost Delta C(j, (R, S))
    #[must_use]
    pub fn delta_c(&self, snp: usize, r: &BTreeSet<usize>, s: &BTreeSet<usize>) -> u32 {
        let (w0_r, w1_r) = self.compute_costs(snp, r);
        let (w0_s, w1_s) = self.compute_costs(snp, s);

        // Given a bipartition (R, S) of F(j), the minimum cost to make position j conflict-free is:
        std::cmp::min(w0_r, w1_r) + std::cmp::min(w0_s, w1_s)

        // Alternatively, under the all heterozygous assumption, where one wants to enforce all SNPs
        // to be heterozygous, the equation becomes:
        // std::cmp::min(w0_r + w1_r, w0_s + w1_s)
    }

    // Write reads matrix to a file in tab-separated format
    pub fn write_reads_matrix(&self, path: &str) -> std::io::Result<()> {
        let mut file = File::create(path)?;

        // Write header row with SNP positions
        for j in 0..self.num_snps {
            write!(file, "\tSNP_{}", j)?;
        }
        writeln!(file)?;

        // Write each read's data
        for (i, (read, conf)) in self.reads.iter().zip(self.confidences.iter()).enumerate() {
            write!(file, "Read_{}", i)?;
            for (allele, qual) in read.iter().zip(conf.iter()) {
                if let Some(allele) = allele {
                    if *allele == 0 {
                        write!(file, "\t0,{}", qual.unwrap())?
                    };
                    if *allele == 1 {
                        write!(file, "\t1,{}", qual.unwrap())?
                    };
                } else {
                    write!(file, "\t-,1")?;
                }
            }
            writeln!(file)?;
        }

        Ok(())
    }
}

// Function to generate all bipartitions of a set
fn generate_bipartitions(set: &BTreeSet<usize>) -> Vec<(BTreeSet<usize>, BTreeSet<usize>)> {
    let mut partitions = vec![];
    let set_vec: Vec<_> = set.iter().collect();

    let num_partitions = 1 << set_vec.len(); // 2^|set|
    for i in 0..num_partitions {
        let mut r = BTreeSet::new();
        let mut s = BTreeSet::new();

        for (j, &elem) in set_vec.iter().enumerate() {
            if i & (1 << j) == 0 {
                r.insert(*elem);
            } else {
                s.insert(*elem);
            }
        }

        partitions.push((r, s));
    }

    partitions
}

// Function to initialize the DP table for SNP 0
fn initialize_dp(
    data: &WMECData,
) -> (
    HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), u32>,
    HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), Option<(BTreeSet<usize>, BTreeSet<usize>)>>,
) {
    let mut dp: HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), u32> = HashMap::new();
    let mut backtrack: HashMap<
        (usize, BTreeSet<usize>, BTreeSet<usize>),
        Option<(BTreeSet<usize>, BTreeSet<usize>)>,
    > = HashMap::new();

    let active_fragments: BTreeSet<usize> = data
        .reads
        .iter()
        .enumerate()
        .filter(|(_, read)| read[0].is_some()) // Only consider fragments covering SNP 0
        .map(|(index, _)| index)
        .collect();

    let partitions = generate_bipartitions(&active_fragments);
    for (r, s) in partitions {
        let cost = data.delta_c(0, &r, &s);
        dp.insert((0, r.clone(), s.clone()), cost);
        backtrack.insert((0, r.clone(), s.clone()), None);
    }

    (dp, backtrack)
}

// Function to update the DP table for each SNP position
fn update_dp_old(
    data: &WMECData,
    dp: &mut HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), u32>,
    backtrack: &mut HashMap<
        (usize, BTreeSet<usize>, BTreeSet<usize>),
        Option<(BTreeSet<usize>, BTreeSet<usize>)>,
    >,
    snp: usize,
) {
    let active_fragments: BTreeSet<usize> = data
        .reads
        .iter()
        .enumerate()
        .filter(|(_, read)| read[snp].is_some()) // Only consider fragments covering SNP
        .map(|(index, _)| index)
        .collect();
    let partitions = generate_bipartitions(&active_fragments);

    for (r, s) in &partitions {
        let delta_cost = data.delta_c(snp, r, s);
        let mut min_cost = u32::MAX;
        let mut best_bipartition = None;

        let prev_active_fragments: BTreeSet<usize> = data
            .reads
            .iter()
            .enumerate()
            .filter(|(_, read)| read[snp - 1].is_some()) // Fragments covering the previous SNP
            .map(|(index, _)| index)
            .collect();

        for (prev_r, prev_s) in generate_bipartitions(&prev_active_fragments) {
            let r_compatible = r
                .intersection(&prev_active_fragments)
                .all(|&x| prev_r.contains(&x));
            let s_compatible = s
                .intersection(&prev_active_fragments)
                .all(|&x| prev_s.contains(&x));

            if r_compatible && s_compatible {
                if let Some(&prev_cost) = dp.get(&(snp - 1, prev_r.clone(), prev_s.clone())) {
                    let current_cost = delta_cost + prev_cost;

                    if current_cost < min_cost {
                        min_cost = current_cost;
                        best_bipartition = Some((prev_r.clone(), prev_s.clone()));
                    }
                }
            }
        }

        dp.insert((snp, r.clone(), s.clone()), min_cost);
        backtrack.insert((snp, r.clone(), s.clone()), best_bipartition);
    }
}

fn update_dp(
    data: &WMECData,
    dp: &mut HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), u32>,
    backtrack: &mut HashMap<
        (usize, BTreeSet<usize>, BTreeSet<usize>),
        Option<(BTreeSet<usize>, BTreeSet<usize>)>,
    >,
    snp: usize,
) {
    let active_fragments: BTreeSet<usize> = data
        .reads
        .iter()
        .enumerate()
        .filter(|(_, read)| read[snp].is_some())
        .map(|(index, _)| index)
        .collect();
    let partitions = generate_bipartitions(&active_fragments);

    // Pre-compute prev_active_fragments since it's used by all iterations
    let prev_active_fragments: BTreeSet<usize> = data
        .reads
        .iter()
        .enumerate()
        .filter(|(_, read)| read[snp - 1].is_some())
        .map(|(index, _)| index)
        .collect();

    // Collect results in parallel
    let results: Vec<_> = partitions
        .par_iter()
        .map(|(r, s)| {
            let delta_cost = data.delta_c(snp, r, s);
            let mut min_cost = u32::MAX;
            let mut best_bipartition = None;

            for (prev_r, prev_s) in generate_bipartitions(&prev_active_fragments) {
                let r_compatible = r
                    .intersection(&prev_active_fragments)
                    .all(|&x| prev_r.contains(&x));
                let s_compatible = s
                    .intersection(&prev_active_fragments)
                    .all(|&x| prev_s.contains(&x));

                if r_compatible && s_compatible {
                    if let Some(&prev_cost) = dp.get(&(snp - 1, prev_r.clone(), prev_s.clone())) {
                        let current_cost = delta_cost + prev_cost;

                        if current_cost < min_cost {
                            min_cost = current_cost;
                            best_bipartition = Some((prev_r.clone(), prev_s.clone()));
                        }
                    }
                }
            }

            ((r.clone(), s.clone()), (min_cost, best_bipartition))
        })
        .collect();

    // Update the hashmaps with the results
    for ((r, s), (min_cost, best_bipartition)) in results {
        dp.insert((snp, r.clone(), s.clone()), min_cost);
        backtrack.insert((snp, r, s), best_bipartition);
    }
}

fn backtrack_haplotypes(
    data: &WMECData,
    dp: &HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), u32>,
    backtrack: &HashMap<
        (usize, BTreeSet<usize>, BTreeSet<usize>),
        Option<(BTreeSet<usize>, BTreeSet<usize>)>,
    >,
) -> (Vec<u8>, Vec<u8>, BTreeSet<usize>, BTreeSet<usize>) {
    let mut best_cost = u32::MAX;
    let mut best_bipartition = None;

    // Restrict processing to reads that span variants within the window.
    let final_active_fragments: BTreeSet<usize> = data
        .reads
        .iter()
        .enumerate()
        .filter(|(_, read)| read[data.num_snps - 1].is_some())
        .map(|(index, _)| index)
        .collect();

    for (r, s) in generate_bipartitions(&final_active_fragments) {
        if let Some(&cost) = dp.get(&(data.num_snps - 1, r.clone(), s.clone())) {
            if cost < best_cost {
                best_cost = cost;
                best_bipartition = Some((r, s));
            }
        }
    }

    let mut haplotype1 = vec![0; data.num_snps];
    let mut haplotype2 = vec![0; data.num_snps];
    let mut current_bipartition = best_bipartition.unwrap();

    for snp in (0..data.num_snps).rev() {
        let (r, s) = &current_bipartition;
        let (w0_r, w1_r) = data.compute_costs(snp, r);
        let (w0_s, w1_s) = data.compute_costs(snp, s);

        if w0_r + w1_s <= w1_r + w0_s {
            haplotype1[snp] = 0;
            haplotype2[snp] = 1;
        } else {
            haplotype1[snp] = 1;
            haplotype2[snp] = 0;
        }

        if snp > 0 {
            if let Some(prev_bipartition) = backtrack
                .get(&(snp, r.clone(), s.clone()))
                .and_then(|x| x.as_ref())
            {
                current_bipartition = prev_bipartition.clone();
            } else {
                // This should not happen if the DP table is correctly filled
                panic!("No valid previous bipartition found for SNP {snp}");
            }
        }
    }

    crate::elog!("wmec score: {}", best_cost);

    (
        haplotype1,
        haplotype2,
        current_bipartition.0,
        current_bipartition.1,
    )
}

// Main function to perform WMEC using dynamic programming
#[must_use]
pub fn phase(data: &WMECData) -> (Vec<u8>, Vec<u8>, BTreeSet<usize>, BTreeSet<usize>) {
    let (mut dp, mut backtrack) = initialize_dp(data);

    for snp in 1..data.num_snps {
        update_dp(data, &mut dp, &mut backtrack, snp);
    }

    backtrack_haplotypes(data, &dp, &backtrack)
}

#[must_use]
pub fn phase_all(
    data: &WMECData,
    window: usize,
    stride: usize,
) -> (Vec<u8>, Vec<u8>, BTreeSet<usize>, BTreeSet<usize>) {
    data.write_reads_matrix("mat.tsv");

    // First, collect all window ranges
    let mut windows: Vec<_> = (0..data.num_snps)
        .step_by(stride)
        .map(|start| (start, (start + window).min(data.num_snps)))
        .collect();

    // Filter out last window if it completely overlaps with previous window
    if windows.len() > 1 {
        if let Some(&(_, prev_end)) = windows.get(windows.len() - 2) {
            if let Some(&(_, last_end)) = windows.last() {
                if last_end == prev_end {
                    windows.pop();
                }
            }
        }
    }

    let pb = crate::utils::default_bounded_progress_bar("Processing windows", windows.len() as u64);

    // Process windows in parallel
    let haplotype_pairs: Vec<_> = windows
        // .iter()
        .par_iter()
        .progress_with(pb)
        .map(|&(start, end)| {
            let (window_indices, window_reads): (Vec<usize>, Vec<Vec<Option<u8>>>) = data
                .reads
                .iter()
                .zip(data.confidences.iter())
                .enumerate()
                .map(|(read_idx, (read, confidence))| {
                    let window_read = read[start..end].to_vec();
                    let window_confidence = confidence[start..end].to_vec();

                    let none_count = window_read.iter().filter(|x| x.is_none()).count();
                    let lowqual_count = window_confidence
                        .iter()
                        .filter(|&&x| x.is_some() && x.unwrap() == 0)
                        .count();

                    (
                        none_count + lowqual_count,
                        read_idx,
                        window_read,
                        window_confidence,
                    )
                })
                .collect::<Vec<_>>()
                .into_iter()
                .sorted_by_key(|(bad_count, _, _, _)| *bad_count)
                .take(12)
                .map(|(_, read_idx, window_read, _)| (read_idx, window_read))
                .unzip();

            let (_, window_confidences): (Vec<usize>, Vec<Vec<Option<u32>>>) = data
                .reads
                .iter()
                .zip(data.confidences.iter())
                .enumerate()
                .map(|(read_idx, (read, confidence))| {
                    let window_read = read[start..end].to_vec();
                    let window_confidence = confidence[start..end].to_vec();

                    let none_count = window_read.iter().filter(|x| x.is_none()).count();
                    let lowqual_count = window_confidence
                        .iter()
                        .filter(|&&x| x.is_some() && x.unwrap() == 0)
                        .count();

                    (
                        none_count + lowqual_count,
                        read_idx,
                        window_read,
                        window_confidence,
                    )
                })
                .collect::<Vec<_>>()
                .into_iter()
                .sorted_by_key(|(bad_count, _, _, _)| *bad_count)
                .take(12)
                .map(|(_, read_idx, _, window_confidences)| (read_idx, window_confidences))
                .unzip();

            crate::elog!("window_reads: {} {} {}", start, end, window_reads.len());
            crate::elog!("window_indices: {} {} {:?}", start, end, window_indices);

            let window_data = WMECData::new(window_reads, window_confidences);

            let (mut dp, mut backtrack) = initialize_dp(&window_data);
            for snp in 1..window_data.num_snps {
                update_dp(&window_data, &mut dp, &mut backtrack, snp);
            }

            let (hap1, hap2, part1, part2) = backtrack_haplotypes(&window_data, &dp, &backtrack);

            let re1: BTreeSet<_> = part1.iter().map(|&i| window_indices[i as usize]).collect();
            let re2: BTreeSet<_> = part2.iter().map(|&i| window_indices[i as usize]).collect();

            crate::elog!("window_re1: {} {} {:?} {:?}", start, end, part1, re1);
            crate::elog!("window_re2: {} {} {:?} {:?}", start, end, part2, re2);

            (hap1, hap2, re1, re2)
        })
        .collect();

    let mut haplotype1 = Vec::new();
    let mut haplotype2 = Vec::new();
    let mut part1 = BTreeSet::new();
    let mut part2 = BTreeSet::new();

    let overlap = window - stride;

    let mut all_read_assignments = BTreeMap::new();

    let mut i = 0;
    for (hap1, hap2, reads1, reads2) in haplotype_pairs {
        // crate::elog!("{}", i);
        // crate::elog!("{:?}", reads1);
        // crate::elog!("{:?}", reads2);

        i += 1;

        if haplotype1.len() == 0 {
            for read_num in reads1.clone() {
                all_read_assignments
                    .entry(read_num)
                    .or_insert(vec![])
                    .push(1);
            }

            for read_num in reads2.clone() {
                all_read_assignments
                    .entry(read_num)
                    .or_insert(vec![])
                    .push(2);
            }

            haplotype1 = hap1;
            haplotype2 = hap2;
            part1 = reads1;
            part2 = reads2;
        } else {
            // Compare overlap regions to determine orientation
            crate::elog!(
                "haplotypes: {} {} {}",
                haplotype1.len(),
                haplotype2.len(),
                overlap
            );

            let new_overlap_len = std::cmp::min(overlap, hap1.len());
            let h1_overlap = &haplotype1[haplotype1.len() - new_overlap_len..];
            let h2_overlap = &haplotype2[haplotype2.len() - new_overlap_len..];
            let new_overlap = &hap1[..new_overlap_len];

            // Count matches between overlapping regions
            let h1_matches = h1_overlap
                .iter()
                .zip(new_overlap.iter())
                .filter(|(a, b)| a == b)
                .count();
            let h2_matches = h2_overlap
                .iter()
                .zip(new_overlap.iter())
                .filter(|(a, b)| a == b)
                .count();

            // Append new haplotypes in correct orientation based on best overlap match
            if h1_matches >= h2_matches {
                haplotype1.extend_from_slice(&hap1[std::cmp::min(overlap, hap1.len())..]);
                haplotype2.extend_from_slice(&hap2[std::cmp::min(overlap, hap2.len())..]);

                for read_num in reads1.clone() {
                    all_read_assignments
                        .entry(read_num)
                        .or_insert(vec![])
                        .push(1);
                }

                for read_num in reads2.clone() {
                    all_read_assignments
                        .entry(read_num)
                        .or_insert(vec![])
                        .push(2);
                }
            } else {
                haplotype1.extend_from_slice(&hap2[std::cmp::min(overlap, hap2.len())..]);
                haplotype2.extend_from_slice(&hap1[std::cmp::min(overlap, hap1.len())..]);

                for read_num in reads1.clone() {
                    all_read_assignments
                        .entry(read_num)
                        .or_insert(vec![])
                        .push(2);
                }

                for read_num in reads2.clone() {
                    all_read_assignments
                        .entry(read_num)
                        .or_insert(vec![])
                        .push(1);
                }
            }
        }
    }

    let mut final_part1 = BTreeSet::new();
    let mut final_part2 = BTreeSet::new();

    crate::elog!("all_read_assignments: {:?}", &all_read_assignments);

    for (read_num, assignments) in all_read_assignments {
        // Count frequency of assignments (1 or 2) for this read
        let most_common = assignments
            .iter()
            .fold(HashMap::new(), |mut counts, &value| {
                *counts.entry(value).or_insert(0) += 1;
                counts
            })
            .into_iter()
            .max_by_key(|&(_, count)| count)
            .map(|(value, _)| value)
            .unwrap_or(0); // Default to 0 if no assignments

        if most_common == 1 {
            final_part1.insert(read_num);
        } else if most_common == 2 {
            final_part2.insert(read_num);
        }
    }

    (haplotype1, haplotype2, final_part1, final_part2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_initialization_manuscript_example_1() {
        // Define the dataset as a matrix of reads
        let reads = vec![
            vec![Some(0), None],    // f0
            vec![Some(1), Some(0)], // f1
            vec![Some(1), Some(1)], // f2
            vec![None, Some(0)],    // f3
        ];

        // Define the confidence degrees
        let confidences = vec![
            vec![Some(5), None],    // f0
            vec![Some(3), Some(2)], // f1
            vec![Some(6), Some(1)], // f2
            vec![None, Some(2)],    // f3
        ];

        // Initialize the dataset
        let data = WMECData::new(reads, confidences);

        // Define the expected values for C(1, ·)
        let expected_values = vec![
            (vec![0, 1, 2], vec![], 5), // C(1, ({f0, f1, f2}, ∅)) = 5
            (vec![0, 1], vec![2], 3),   // C(1, ({f0, f1}, {f2})) = 3
            (vec![0, 2], vec![1], 5),   // C(1, ({f0, f2}, {f1})) = 5
            (vec![1, 2], vec![0], 0),   // C(1, ({f1, f2}, {f0})) = 0
        ];

        // Check that the values in the dynamic programming table match the expected values
        for (r, s, expected_cost) in expected_values {
            let r_set: BTreeSet<usize> = r.into_iter().collect();
            let s_set: BTreeSet<usize> = s.into_iter().collect();
            let cost = data.delta_c(0, &r_set, &s_set);
            assert_eq!(
                cost, expected_cost,
                "Cost for partition ({:?}, {:?}) does not match expected",
                r_set, s_set
            );
        }
    }

    #[test]
    fn test_recurrence_manuscript_example_2() {
        // Define the dataset as a matrix of reads
        let reads = vec![
            vec![Some(0), None],    // f0
            vec![Some(1), Some(0)], // f1
            vec![Some(1), Some(1)], // f2
            vec![None, Some(0)],    // f3
        ];

        // Define the confidence degrees
        let confidences = vec![
            vec![Some(5), None],    // f0
            vec![Some(3), Some(2)], // f1
            vec![Some(6), Some(1)], // f2
            vec![None, Some(2)],    // f3
        ];

        // Initialize the dataset
        let data = WMECData::new(reads, confidences);

        // Define the expected values for C(2, ·)
        let expected_values = vec![
            (vec![1, 2, 3], vec![], 1), // C(2, ({f1, f2, f3}, ∅)) = 1
            (vec![1, 2], vec![3], 1),   // C(2, ({f1, f2}, {f3})) = 1
            (vec![1, 3], vec![2], 3),   // C(2, ({f1, f3}, {f2})) = 3
            (vec![2, 3], vec![1], 4),   // C(2, ({f2, f3}, {f1})) = 4
        ];

        let (mut dp, mut backtrack) = initialize_dp(&data);

        for snp in 1..data.num_snps {
            update_dp(&data, &mut dp, &mut backtrack, snp);
        }

        // Verify that the results comport with expected_values
        for (r, s, expected_cost) in expected_values {
            let r_set: BTreeSet<usize> = r.into_iter().collect();
            let s_set: BTreeSet<usize> = s.into_iter().collect();
            let actual_cost = dp
                .get(&(1, r_set.clone(), s_set.clone()))
                .unwrap_or(&u32::MAX);
            assert_eq!(
                *actual_cost, expected_cost,
                "Cost for partition ({:?}, {:?}) does not match expected",
                r_set, s_set
            );
        }
    }

    /// This test case is based on Figure 1 from the paper:
    /// "WhatsHap: fast and accurate read-based phasing"
    /// by Marcel Martin et al.
    /// DOI: https://doi.org/10.1101/085050
    ///
    /// The figure illustrates a small example of the weighted minimum error correction problem,
    /// which is the core algorithmic component of WhatsHap.
    #[test]
    fn test_whatshap_manuscript_figure_1() {
        // Define the dataset as a matrix of reads
        let reads = vec![
            vec![Some(0), None, None, Some(0), Some(1)], // f0
            vec![Some(1), Some(0), Some(0), None, None], // f1
            vec![Some(1), Some(1), Some(0), None, None], // f2
            vec![None, Some(0), Some(0), Some(1), None], // f3
            vec![None, None, Some(1), Some(0), Some(1)], // f4
            vec![None, None, None, Some(0), Some(1)],    // f5
        ];

        // Define the confidence degrees
        let confidences = vec![
            vec![Some(32), None, None, Some(34), Some(17)], // f0
            vec![Some(15), Some(25), Some(13), None, None], // f1
            vec![Some(7), Some(3), Some(15), None, None],   // f2
            vec![None, Some(12), Some(23), Some(29), Some(31)], // f3
            vec![None, None, Some(25), Some(17), Some(19)], // f4
            vec![None, None, None, Some(20), Some(10)],     // f5
        ];

        // Initialize the dataset
        let data = WMECData::new(reads, confidences);

        // Perform the WMEC algorithm using dynamic programming
        let (haplotype1, haplotype2, part1, part2) = phase(&data);

        let expected_haplotype1 = vec![0, 1, 1, 0, 1];
        let expected_haplotype2 = vec![1, 0, 0, 1, 0];

        assert_eq!(
            haplotype1, expected_haplotype1,
            "Haplotype 1 does not match expected"
        );
        assert_eq!(
            haplotype2, expected_haplotype2,
            "Haplotype 2 does not match expected"
        );
    }
}
