use clap::{Args, Parser, Subcommand};

mod prepare;
mod assemble;
mod join;
mod impute;

use skydive;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
#[clap(propagate_version = true)]
struct CLI {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Prepare a single long-read dataset for inclusion into a hidive graph.
    Prepare(Prepare),

    /// Assemble a hidive graph for a single long-read or short-read dataset.
    Assemble(Assemble),

    /// Join multiple hidive graphs.
    Join(Join),

    /// Impute the presence of missing k-mers in a joined hidive graph.
    Impute(Impute),
}

#[derive(Args)]
struct Prepare {
    /// first number
    number_one: i32,

    /// second number
    number_two: i32
}

#[derive(Args)]
struct Assemble {
    /// first number
    number_one: i32,

    /// second number
    number_two: i32
}

#[derive(Args)]
struct Join {
    /// first number
    number_one: i32,

    /// second number
    number_two: i32
}

#[derive(Args)]
struct Impute {
    /// first number
    number_one: i32,

    /// second number
    number_two: i32
}

fn main() {
    let num = 10;
    println!("Hello, world! {} plus one is {}!",
        num,
        skydive::add_one(num)
    );

    let cli = CLI::parse();

    match cli.command {
        Commands::Prepare(prepare) => prepare::start(&num, &2),
        Commands::Assemble(assemble) => assemble::start(&num, &3),
        Commands::Join(join) => join::start(&num, &4),
        Commands::Impute(impute) => impute::start(&num, &5),
    };
}
