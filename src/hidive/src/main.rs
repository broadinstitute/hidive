use skydive;

fn main() {
    let num = 10;
    println!("Hello, world! {} plus one is {}!",
        num,
        skydive::add_one(num)
    );
}
