fn build_constants() -> miette::Result<()>  {
    let out_dir: std::path::PathBuf = String::from("src/")
    .into();
    let mut code = Vec::new();
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=K");
    let k: usize = std::env::var("K")
        .unwrap_or_else(|_| "31".into())
        .parse()
        .expect("Failed to parse K");
    assert!(k >= 1, "K must be ≥ 1");
    assert!(k <= 61, "K must be ≤ 59");
    assert!(k % 2 == 1, "K must be odd");
    code.push(format!("pub const K: usize = {k};"));

    let kmer_bits = 2 * k;
    let kt = select_type(kmer_bits);
    code.push(format!("pub type KT = {kt};"));

    std::fs::write(out_dir.join("constants.rs"), code.join("\n"))
        .expect("Failed to write const file");
    Ok(())
}

fn select_type(n_bits: usize) -> &'static str {
    match n_bits.next_power_of_two() {
        1 | 2 | 4 | 8 => "u8",
        16 => "u16",
        32 => "u32",
        64 => "u64",
        128 => "u128",
        _ => panic!("Cannot fit {n_bits} bits in a primitive type"),
    }
}

fn main() -> miette::Result<()> {
    println!("cargo:rerun-if-changed=build.rs");
    build_constants()
}