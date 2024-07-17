use std::collections::HashMap;
use std::collections::HashSet;
use serde_json::Value;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;
use std::sync::Mutex;
use std::fs::File;
use flate2::read::GzDecoder;

#[derive(Debug)]
pub struct GraphicalGenome {
    pub anchor: HashMap<String, Value>,
    pub edges: HashMap<String, Value>,
    pub outgoing: HashMap<String, Vec<String>>,
    pub incoming: HashMap<String, Vec<String>>,
}

impl GraphicalGenome {
    pub fn load_graph(filename: &str) -> io::Result<GraphicalGenome> {
        let file = File::open(filename)?;
        let reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        let mut anchor_dict = HashMap::new();
        let mut edge_dict = HashMap::new();
        let mut outgoing: HashMap<String, Vec<String>> = HashMap::new();
        let mut incoming: HashMap<String, Vec<String>> = HashMap::new();
        for line in reader.lines() {
            let line = line?.trim_end().to_string();
            if line.starts_with('S') {
                let itemlist: Vec<&str> = line.split('\t').collect();
                let name = itemlist[1];
                let seq = itemlist[2];
                let annotation = &itemlist[3][5..];

                let json_value: Value = serde_json::from_str(annotation).unwrap();

                let mut value = json_value;

                value["seq"] = Value::String(seq.to_string());

                if name.starts_with('A') {
                    anchor_dict.insert(name.to_string(), value);
                } else if name.starts_with('E') {
                    edge_dict.insert(name.to_string(), value);
                }
            } else if line.starts_with('L') {
                let itemlist: Vec<&str> = line.split('\t').collect();
                let src = itemlist[1];
                let dst = itemlist[3];
                outgoing
                    .entry(src.to_string())
                    .or_default()
                    .push(dst.to_string());
                incoming
                    .entry(dst.to_string())
                    .or_default()
                    .push(src.to_string());
            }
        }
        Ok(GraphicalGenome {
            anchor: anchor_dict,
            edges: edge_dict,
            outgoing: outgoing,
            incoming: incoming,
        })
    }
    // Method to extract a single sample graph
    pub fn extract_single_sample_graph(
        &self,
        df_single_sample: &HashMap<String, HashMap<String, usize>>, // Assuming df_single_sample is a HashMap with anchor as key and another HashMap of readset and edge index as value
        sample: &str,
    ) -> (HashMap<String, Value>, HashMap<String, Value>, HashMap<String, Vec<String>>, HashMap<String, Vec<String>>) {
        let mut new_edges = HashMap::new();
        let mut new_incoming = HashMap::new();
        let mut new_outgoing = HashMap::new();

        for (anchor, readsets) in df_single_sample {
            if let Some(outgoinglist) = self.outgoing.get(anchor) {
                for (read, &edge_index) in readsets {
                    if let Some(edgename) = outgoinglist.get(edge_index) {
                        let default_edge_data = serde_json::json!({}); // Create a longer-lived value
                        let edge_data = self.edges.get(edgename).unwrap_or(&default_edge_data);
                        let mut edge_value = edge_data.clone();

                        // Update sequence data
                        edge_value["seq"] = edge_data["seq"].clone();

                        // Update reads and strain
                        let edge_object = edge_value.as_object_mut().unwrap();
                        {
                            let reads = edge_object.entry("reads").or_insert_with(|| serde_json::json!([]));
                            reads.as_array_mut().unwrap().push(serde_json::json!(read));
                        }
                        {
                            let strain = edge_object.entry("strain").or_insert_with(|| serde_json::json!(HashSet::<String>::new()));
                            strain.as_array_mut().unwrap().push(serde_json::json!(sample));
                        }
                        
                        new_edges.insert(edgename.to_string(), edge_value);

                        // Update incoming and outgoing
                        if let Some(dst) = outgoinglist.get(0) { // Assuming the first element is the destination
                            new_incoming.entry(edgename.to_string()).or_insert_with(Vec::new).push(anchor.to_string());
                            new_incoming.entry(dst.to_string()).or_insert_with(Vec::new).push(edgename.to_string());
                            new_outgoing.entry(anchor.to_string()).or_insert_with(Vec::new).push(edgename.to_string());
                            new_outgoing.entry(edgename.to_string()).or_insert_with(Vec::new).push(dst.to_string());
                        }
                    }
                }
            }
        }

        (self.anchor.clone(), new_edges, new_outgoing, new_incoming)
    }
}
